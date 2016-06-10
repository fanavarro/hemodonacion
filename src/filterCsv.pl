use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use myUtils::CsvManager;
#Get input and output file from input parameters
my $output;
my $input;
if (scalar @ARGV == 2){
	$input = $ARGV[0];
	$output = $ARGV[1];
} else {
	print "Usage: perl filterCsv.pl inputFile outputFile";
	exit;
}

if (! -e $input){
	print "File $input does not exist.\n";
	exit;
}

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

my @fields = qw(CHROMOSOME GENE_ID GENE_NAME TRANSCRIPT_ID TRANSCRIPT_NAME TRANSCRIPT_BIOTYPE CDS_ERRORS PROTEIN_ID VARIATION_NAME SOURCE TRANSCRIPT_VARIATION_ALLELE_DBID MINOR_ALLELE_FREQUENCY CODON_CHANGE AMINOACID_CHANGE FIRST_MET_POSITION STOP_CODON_POSITION MUTATED_SEQUENCE_LENGTH READING_FRAME_STATUS CONSEQUENCE PHENOTYPE SO_TERM SIFT POLYPHEN PUBLICATIONS);
my $in_csv = myUtils::CsvManager->new (
	fields    => \@fields,
	csv_separator   => "\t",
	in_field_separator    => '-',
	file_name => $input,
	mode => '<'
);
my $out_csv = myUtils::CsvManager->new (
	fields    => \@fields,
	csv_separator   => "\t",
	in_field_separator    => '-',
	file_name => $output,
	mode => '>'
);

my $total_entries = $in_csv->myUtils::CsvManager::countEntries();
my $count = 1;
while ((my %entry = $in_csv->myUtils::CsvManager::readEntry())){
	print $count . "/" . $total_entries . "\n";
	$count = $count + 1;
	my $transcript_stable_id = $entry{'TRANSCRIPT_ID'};
	if (apply_filter($transcript_stable_id)){
		$out_csv->myUtils::CsvManager::writeEntry(%entry);
		print "$transcript_stable_id pass the filter\n";
	} else {
		print "$transcript_stable_id do not pass the filter\n";
	}
}

$in_csv->myUtils::CsvManager::close();
$out_csv->myUtils::CsvManager::close();

# Functions that receives a transcript identifier
# and returns if it pass the filter or not.
sub apply_filter{
	my $transcript_stable_id = $_[0];
	my $transcript = $transcript_adaptor->fetch_by_stable_id($transcript_stable_id);
	return (is_canonical($transcript) &&
		is_known($transcript) &&
		has_translation($transcript));
}

sub is_canonical{
	my $transcript = $_[0];
	return $transcript->is_canonical;
}
sub is_known{
	my $transcript = $_[0];
	return $transcript->is_known;
}

sub has_translation{
	my $transcript = $_[0];
	return $transcript->translation();
}
