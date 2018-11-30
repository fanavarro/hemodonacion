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
my $output = "";
if (scalar @ARGV == 2){
	$input = $ARGV[0];
	$output = $ARGV[1];
} else {
	print "Usage: perl enrichExome.pl inputFile outputFile\n";
	exit;
}

# Needed to write to nohup.out in real time
STDOUT->autoflush(1);

my $ensemblUtils = myUtils::EnsemblUtils->instance();

# CSV file params
my %inputParams = (
	csv_separator   => "\t",
	in_field_separator    => '-',
	file_name => $input,
	mode => '<'
);
my $inputCSV = myUtils::CsvManager->new (%inputParams);

my @outputFields = qw(CHROMOSOME TRANSCRIPT_ID GENE_ID GENE_NAME PROTEIN_ID VARIATION_ID REF ALT SELECTED_ALT SIGNAL_PEPTIDE_START SIGNAL_PEPTIDE_END METS_IN_5_UTR APPROACH1_MET_POSITION APPROACH1_STOP_CODON_POSITION APPROACH1_MUTATED_SEQUENCE_LENGTH APPROACH1_READING_FRAME_STATUS APPROACH1_SIGNAL_PEPTIDE_CONSERVATION APPROACH2_MET_POSITION APPROACH2_INIT_CODON APPROACH2_STOP_CODON_POSITION APPROACH2_MUTATED_SEQUENCE_LENGTH APPROACH2_SCORE APPROACH2_READING_FRAME_STATUS APPROACH2_SIGNAL_PEPTIDE_CONSERVATION APPROACH3_MET_POSITION APPROACH3_INIT_CODON APPROACH3_STOP_CODON_POSITION APPROACH3_MUTATED_SEQUENCE_LENGTH APPROACH3_SCORE APPROACH3_READING_FRAME_STATUS APPROACH3_SIGNAL_PEPTIDE_CONSERVATION);
my %outputParams = (
	fields    => \@outputFields,
	csv_separator   => "\t",
	in_field_separator    => '-',
	file_name => $output,
	mode => '>'
);
my $outputCSV = myUtils::CsvManager->new (%outputParams);
$outputCSV->myUtils::CsvManager::writeHeader();
my @outputEntries;

# Read input csv rows
my $rows = $inputCSV->myUtils::CsvManager::readCsv();
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
	
	# Get transcripts that overlaps with the variation position.
	my $transcripts = getOverlappingTranscripts($chrom, $variationGenomicPosition);
	
	# For each transcript check
	foreach my $transcript ( @{$transcripts} ){
		# Skip transcript if there are error annotations in CDS.
		my $cds_errors = $ensemblUtils->get_cds_errors($transcript);
		if($cds_errors ne ""){
			next;
		}
		# Check if is protein coding
		if($transcript->translation()){
			# Check if variation is on the start codon
			if (overlapsInitiationCodon($chrom, $variationGenomicPosition, $transcript)){
				# Get transcript info
				my $signal_peptide_info = $ensemblUtils->get_signal_peptide_info($transcript);
				my $five_prime_utr = $transcript->five_prime_utr;
				if(defined($five_prime_utr)){
					my $five_prime_utr_seq = $five_prime_utr->seq;
					my $ref = $row->{'REF'};
					my $posRef = $row->{'POS'};
					my $transcript_info = generateTranscriptInfo($transcript, $row, $signal_peptide_info);
					my $cdna = $transcript->seq->seq;
					my $cds = $transcript->translateable_seq;
					my @alleles = split (/,/, $row->{'ALT'});
					# For each allele, get variation allele info
					foreach my $allele (@alleles){
						my $mutated_cdna = getMutatedCdna($transcript, $ref, $allele, $posRef);
						#print "\n> CDNA\n$cdna\n> CDS\n$cds\n> Mutated CDNA\n$mutated_cdna\n";
						my $alleleInfo = generateAlleleInfo($cdna, $mutated_cdna, $cds, $five_prime_utr_seq, $signal_peptide_info, isFivePrimeAffected());
						my %entry = (%{$alleleInfo}, %{$transcript_info});
						$entry{'SELECTED_ALT'} = $allele;
						push(@outputEntries, \%entry);
					}
				}
			}
		}
	}
}
$outputCSV->myUtils::CsvManager::writeEntries(\@outputEntries);
$outputCSV->myUtils::CsvManager::close();
$inputCSV->myUtils::CsvManager::close();


sub generateTranscriptInfo {
	my $transcript = shift;
	my $row = shift;
	my $signal_peptide_info = shift;
	my %entry;
	$entry{'CHROMOSOME'} = $row->{'CHROM'};
	$entry{'TRANSCRIPT_ID'} = $transcript->display_id;
	$entry{'GENE_ID'} = $transcript->get_Gene->stable_id;
	$entry{'GENE_NAME'} = $transcript->get_Gene->external_name;
	$entry{'PROTEIN_ID'} = $transcript->translation->display_id;
	$entry{'VARIATION_ID'} = $row->{'ID'};
	$entry{'REF'} = $row->{'REF'};
	$entry{'ALT'} = $row->{'ALT'};
	$entry{'SIGNAL_PEPTIDE_START'} = $signal_peptide_info->{'START'} if (defined($signal_peptide_info));
    $entry{'SIGNAL_PEPTIDE_END'} = $signal_peptide_info->{'END'} if (defined($signal_peptide_info));
    return \%entry;
}
sub generateAlleleInfo{
	my $cdna = shift;
	my $mutated_cdna = shift;
	my $cds = shift;
	my $five_prime_utr_seq = shift;
	my $signal_peptide_info = shift;
	my $five_prime_affected = shift;
	
	if(!defined($five_prime_affected)){
		$five_prime_affected = 0;
	}
	
	my $seq_info = myUtils::SeqUtils::get_met_mutation_info($cdna, $cds, $mutated_cdna, $five_prime_affected);
	my $pwm_info = $ensemblUtils->get_match_pwm_info_with_seqs($cdna, $mutated_cdna, $cds, $five_prime_affected);
	my $noderer_info = $ensemblUtils->get_noderer_info_with_seqs($cdna, $mutated_cdna, $cds, $five_prime_affected);
	my $five_prime_utr_info = $ensemblUtils->get_five_prime_utr_info_with_seqs($cdna, $mutated_cdna, $five_prime_utr_seq);
	
	my %entry;
	$entry{'METS_IN_5_UTR'} = $five_prime_utr_info;
	
	$entry{'APPROACH1_MET_POSITION'} = $seq_info->{'met_position'};
	$entry{'APPROACH1_STOP_CODON_POSITION'} = $seq_info->{'stop_codon_position'};
	$entry{'APPROACH1_MUTATED_SEQUENCE_LENGTH'} = $seq_info->{'seq_length'};
	$entry{'APPROACH1_READING_FRAME_STATUS'} = $seq_info->{'reading_frame'};
	$entry{'APPROACH1_SIGNAL_PEPTIDE_CONSERVATION'} = $ensemblUtils->get_signal_peptide_conservarion($signal_peptide_info->{'START'}, $signal_peptide_info->{'END'}, $seq_info->{'met_position'});
	
	
    $entry{'APPROACH2_MET_POSITION'} = $noderer_info->{'met_position'};
    $entry{'APPROACH2_INIT_CODON'} = $noderer_info->{'init_codon'};
    $entry{'APPROACH2_STOP_CODON_POSITION'} = $noderer_info->{'stop_codon_position'};
    $entry{'APPROACH2_MUTATED_SEQUENCE_LENGTH'} = $noderer_info->{'seq_length'};
    $entry{'APPROACH2_SCORE'} = $noderer_info->{'score'};
    $entry{'APPROACH2_READING_FRAME_STATUS'} = $noderer_info->{'reading_frame'};
    $entry{'APPROACH2_SIGNAL_PEPTIDE_CONSERVATION'} = $ensemblUtils->get_signal_peptide_conservarion($signal_peptide_info->{'START'}, $signal_peptide_info->{'END'}, $noderer_info->{'met_position'});
    
    
    $entry{'APPROACH3_MET_POSITION'} = $pwm_info->{'met_position'};
    $entry{'APPROACH3_INIT_CODON'} = $pwm_info->{'init_codon'};
    $entry{'APPROACH3_STOP_CODON_POSITION'} = $pwm_info->{'stop_codon_position'};
    $entry{'APPROACH3_MUTATED_SEQUENCE_LENGTH'} = $pwm_info->{'seq_length'};
    $entry{'APPROACH3_SCORE'} = $pwm_info->{'score'};
    $entry{'APPROACH3_READING_FRAME_STATUS'} = $pwm_info->{'reading_frame'};
    $entry{'APPROACH3_SIGNAL_PEPTIDE_CONSERVATION'} = $ensemblUtils->get_signal_peptide_conservarion($signal_peptide_info->{'START'}, $signal_peptide_info->{'END'}, $pwm_info->{'met_position'});
	return \%entry;
}

sub isFivePrimeAffected{
	# todo
	return 0;
}

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


