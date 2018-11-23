use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use myUtils::CsvManager;

#Get input file from parameters
my $input = "";
if (scalar @ARGV == 1){
	$input = $ARGV[0];
} else {
	print "Usage: perl enrichExome.pl inputFile";
	exit;
}

# Needed to write to nohup.out in real time
#STDOUT->autoflush(1);

# Registry configuration
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous',
);
# Set the flag to make sure that the database connection is dropped if not being used on each database.
# $registry->set_disconnect_when_inactive();
# Set the flag to make sure that the database connection is not lost before it's used. This is useful for long running jobs (over 8hrs).
$registry->set_reconnect_when_lost();


# Get the adaptor to get the Transcript, slices and transcript variation from the database
my $transcript_adaptor = $registry->get_adaptor( 'homo_sapiens', 'core', 'transcript' );
my $slice_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Slice' );
my $trv_adaptor = $registry->get_adaptor( 'homo_sapiens', 'variation', 'transcriptvariation' );

getOverlappingTranscripts('chr1', 201896455);


# Get a list of transcripts that overlap the position passed as argument.
# param 1 -> The chromosome.
# param 2 -> The position in the chromosome.
# return a list transcripts
sub getOverlappingTranscripts{
	my $chromosome = shift;
	my $position = shift;
	my $location = $chromosome . ':' . $position . '-' . $position;
	my $slice = $slice_adaptor->fetch_by_location($location, 'chromosome');
	print 'inicio slice = ' . $slice->start() . "\n";
	print 'fin slice = ' . $slice->end() . "\n";
	my $transcripts = $transcript_adaptor->fetch_all_by_Slice($slice);
	print "Transcritos: \n";
	foreach my $transcript ( @{$transcripts} ) {
		my $id = $transcript->stable_id();
		if($transcript->translation()){
			$transcript = $transcript->transform('chromosome');
			print  $id . "\n";
			my $initCodonCoordinateList = getInitiationCodonGenomicCoords($transcript);
			
			if(defined($initCodonCoordinateList)){
				foreach my $hash ( @{$initCodonCoordinateList} ){
					my $str = "chr" . $hash->{'chr'} . ":" . $hash->{'start'} . "-" . $hash->{'end'} . "\tstrand = " . $hash->{'strand'};
					print $str . "\n";
				}
			} else {
				print $id . "No se encontro inicio.\n";
			}
		} else{
			print $id . " no codificante.\n";
		}
	}
	
}

sub isInitiationCodon {
	my $chromosome = shift;
	my $position = shift;
	my $transcript = shift;
	my $codingRegionStart = $transcript->cdna_coding_start();
	if(!defined($codingRegionStart)){
		return 0;
	}
}

sub getInitiationCodonGenomicCoords {
	my $transcript = shift;
	my $mapper = Bio::EnsEMBL::TranscriptMapper->new($transcript);
	my @coords = $mapper->pep2genomic(1, 1);
	my $initiationCodonCoords = ();
	foreach my $coord (@coords) {
		if ($coord->isa("Bio::EnsEMBL::Mapper::Coordinate")){
			print "genomic start is " . $coord->start . " and ends at " . $coord->end."\n";
			my $initiationCodonCoord = {};
			$initiationCodonCoord->{'start'} = $coord->start;
			$initiationCodonCoord->{'end'} = $coord->end;
			#$initiationCodonCoord->{'chr'} = $coord->name(); # No se encuentra el metodo name... lo mismo es de otra version.
			$initiationCodonCoord->{'chr'} = $transcript->seq_region_name();
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


