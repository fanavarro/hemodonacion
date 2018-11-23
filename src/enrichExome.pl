use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use myUtils::CsvManager;
use myUtils::SeqUtils;
use myUtils::EnsemblUtils;
use IO::Handle;
use Data::Dumper;
use v5.10.0;


#Get input file from parameters
my $input = "";
if (scalar @ARGV == 1){
	$input = $ARGV[0];
} else {
	print "Usage: perl enrichExome.pl inputFile";
	exit;
}

# Needed to write to nohup.out in real time
STDOUT->autoflush(1);

my $ensemblUtils = myUtils::EnsemblUtils->instance();
#my $matchpwm     = myUtils::MatchPWM->instance();

# CSV file params
my %params = (
	csv_separator   => "\t",
	in_field_separator    => '-',
	file_name => $input,
	mode => '<'
);
my $csvManager = myUtils::CsvManager->new (%params);

# Read csv rows
my $rows = $csvManager->myUtils::CsvManager::readCsv();
my $rowCount = 0;
my $totalRows = scalar @{$rows};
foreach my $row (@{$rows}){
	$rowCount++;
	print "$rowCount/$totalRows\n";
	my $chrom = $row->{'CHROM'};
	if(isAlternative($chrom)){
		next;
	}
	my $variationGenomicPosition = $row->{'POS'};
	print "Variant position -> " . $chrom . ":" . $variationGenomicPosition . "\n";
	
	# Get transcripts that overlaps with the variation position.
	my $transcripts = getOverlappingTranscripts($chrom, $variationGenomicPosition);
	
	# Check if that position is affecting the initiation codon of these transcripts.
	foreach my $transcript ( @{$transcripts} ){
		if($transcript->translation()){
			print $transcript->stable_id() . "\n";
			if (overlapsInitiationCodon($chrom, $variationGenomicPosition, $transcript)){
				my $ref = $row->{'REF'};
				my $posRef = $row->{'POS'};
				my @alleles = split (/,/, $row->{'ALT'});
				foreach my $allele (@alleles){
					my $cdna = $transcript->seq->seq;
					my $cds = $transcript->translateable_seq;
					my $mutatedCdna = getMutatedCdna($transcript, $ref, $allele, $posRef);
					print "$ref -> $allele\n";
					#print "$cdna\n-\n$mutatedCdna\n";
				}
			}
		}
	}
}
$csvManager->myUtils::CsvManager::close();

sub getMutatedCdna {
	my $transcript = shift;
	my $ref = shift;
	my $alt = shift;
	my $pos = shift;
	my $strand = $transcript->strand();
	my $cdna = $transcript->seq->seq;
	my $cds = $transcript->translateable_seq;
	my $mapper = Bio::EnsEMBL::TranscriptMapper->new($transcript);
	my @coords = $mapper->genomic2cdna($pos, $pos, $transcript->strand());
	my $refCdnaPos = $coords[0]->start - 1;
	my $refLength = length($ref);
	if(defined ($strand) && $strand == -1){
		$ref = myUtils::SeqUtils::getReverseComplement($ref);
		$alt = myUtils::SeqUtils::getReverseComplement($alt);
	}
	my $mutatedCdna = $cdna;
	substr($mutatedCdna, $refCdnaPos, $refLength) = $alt;
	return $mutatedCdna;
}

# Check if a chromosome name references an alternative chromosome (as 'chr7_KI270803v1_alt')
sub isAlternative{
	my $chrName = shift;
	# If name contains the '_' character.
	if($chrName =~ m/.*_.*/){
		return 1;
	}
	return 0;
}
# Get a list of transcripts that overlap the position passed as argument.
# param 1 -> The chromosome.
# param 2 -> The position in the chromosome.
# return a list transcripts
sub getOverlappingTranscripts{
	my $chromosome = shift;
	my $position = shift;
	my $location = $chromosome . ':' . $position . '-' . $position;
	my $slice = $ensemblUtils->{slice_adaptor}->fetch_by_location($location, 'chromosome');
	my $transcripts = $ensemblUtils->{transcript_adaptor}->fetch_all_by_Slice($slice);
	return $transcripts;
	
}

# Checks if the position given by a chromosome and a position overlaps the initiation codon of the transcript.
# param 1 -> The chromosome.
# param 2 -> The position in the chromosome.
# param 3 -> The transcript.
# return true if the position is on the initiation codon of the transcript and false in other case.
sub overlapsInitiationCodon {
	my $chromosome = shift;
	my $position = shift;
	my $transcript = shift;
	my $initCodonCoordinateList = getInitiationCodonGenomicCoords($transcript);
	if(defined($initCodonCoordinateList)){
		foreach my $hash ( @{$initCodonCoordinateList} ){
			#my $str = $hash->{'chr'} . ":" . $hash->{'start'} . "-" . $hash->{'end'} . "\tstrand = " . $hash->{'strand'};
			#print $str . "\n";
			if(overlapsPosition($chromosome, $position, $hash)){
				#print "YEAH!";
				return 1;
			}
		}
	}
	return 0;
}

sub getInitiationCodonGenomicCoords {
	my $transcript = shift;
	my $mapper = Bio::EnsEMBL::TranscriptMapper->new($transcript);
	my @coords = $mapper->pep2genomic(1, 1);
	my $initiationCodonCoords = ();
	foreach my $coord (@coords) {
		if ($coord->isa("Bio::EnsEMBL::Mapper::Coordinate")){
			#print "genomic start is " . $coord->start . " and ends at " . $coord->end."\n";
			my $initiationCodonCoord = {};
			$initiationCodonCoord->{'start'} = $coord->start;
			$initiationCodonCoord->{'end'} = $coord->end;
			#$initiationCodonCoord->{'chr'} = $coord->name(); # No se encuentra el metodo name... lo mismo es de otra version.
			$initiationCodonCoord->{'chr'} = 'chr' . $transcript->seq_region_name();
			$initiationCodonCoord->{'strand'} = $coord->strand;
			push @{$initiationCodonCoords}, $initiationCodonCoord;
		}
		else{
			print "GAP from" . $coord->start . " to " . $coord->end . "\t(" . $transcript->stable_id() . ")\n";
		}
 	}
 	if(scalar @{ $initiationCodonCoords } == 0){
 		$initiationCodonCoords = undef;
 	}
 	return $initiationCodonCoords;
}

# Check if a coord object contains a single position given by a chromosome and a genomic position.
# param 1 -> The chromosome.
# param 2 -> The position in the chromosome.
# param 3 -> Coord hash, that includes start, end, chr and strand.
# return a list transcripts
sub overlapsPosition{
	my $chromosome = shift;
	my $position = shift;
	my $coord = shift;
	# If single position and coord are in the same chromosome
	if($chromosome eq $coord->{'chr'}){
		# If the position is between the start and the end
		if($position >= $coord->{'start'} && $position <= $coord->{'end'}){
			return 1;
		}
	}
	
	# The position is not in the coord.
	return 0;
}


