use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use myUtils::CsvManager;

STDOUT->autoflush(1);
#Get input and output file from input parameters
my $output;
my $input;
if (scalar @ARGV == 2){
	$input = $ARGV[0];
	$output = $ARGV[1];
} else {
	print "Usage: perl add_column.pl inputFile outputFile";
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

my @fields1 = qw(CHROMOSOME GENE_ID GENE_NAME TRANSCRIPT_ID TRANSCRIPT_NAME TRANSCRIPT_BIOTYPE CDS_ERRORS PROTEIN_ID VARIATION_NAME SOURCE TRANSCRIPT_VARIATION_ALLELE_DBID MINOR_ALLELE_FREQUENCY CODON_CHANGE AMINOACID_CHANGE FIRST_MET_POSITION STOP_CODON_POSITION MUTATED_SEQUENCE_LENGTH READING_FRAME_STATUS CONSEQUENCE PHENOTYPE SO_TERM SIFT POLYPHEN PUBLICATIONS);
my @fields2 = qw(CHROMOSOME GENE_ID GENE_NAME TRANSCRIPT_ID TRANSCRIPT_REFSEQ_ID TRANSCRIPT_NAME TRANSCRIPT_BIOTYPE CDS_ERRORS PROTEIN_ID VARIATION_NAME SOURCE TRANSCRIPT_VARIATION_ALLELE_DBID MINOR_ALLELE_FREQUENCY CODON_CHANGE AMINOACID_CHANGE FIRST_MET_POSITION STOP_CODON_POSITION MUTATED_SEQUENCE_LENGTH READING_FRAME_STATUS CONSEQUENCE PHENOTYPE SO_TERM SIFT POLYPHEN PUBLICATIONS);
my $in_csv = myUtils::CsvManager->new (
	fields    => \@fields1,
	csv_separator   => "\t",
	in_field_separator    => '-',
	file_name => $input,
	mode => '<'
);
my $out_csv = myUtils::CsvManager->new (
	fields    => \@fields2,
	csv_separator   => "\t",
	in_field_separator    => '-',
	file_name => $output,
	mode => '>'
);


my @old_entries = @{$in_csv->myUtils::CsvManager::readCsv()};
my $total_entries = scalar(@old_entries);
my $count = 1;
$out_csv->myUtils::CsvManager::writeHeader();
foreach my $old_entry (@old_entries){
	print $count . "/" . $total_entries . "\n";
	$count = $count + 1;
	my $transcript_stable_id = $old_entry->{'TRANSCRIPT_ID'};
	my %new_entry = get_new_entry($old_entry);
	$out_csv->myUtils::CsvManager::writeEntry(%new_entry);
}

$in_csv->myUtils::CsvManager::close();
$out_csv->myUtils::CsvManager::close();


# Receives the old entry and returns
# a new entry with new fields
sub get_new_entry{
	my %old_entry = %{$_[0]};
	my %new_entry = %old_entry;
	my $transcript = $transcript_adaptor->fetch_by_stable_id($old_entry{'TRANSCRIPT_ID'});
	my @db_entries = @{$transcript->get_all_DBEntries('RefSeq_mRNA')};
	my $transcript_refseq_id = '';
	foreach my $dbentry ( @db_entries ){
		if ($dbentry->display_id){
			$transcript_refseq_id = $transcript_refseq_id . $dbentry->display_id . '-';
		}
	}
	chop($transcript_refseq_id);
	$new_entry{'TRANSCRIPT_REFSEQ_ID'} = $transcript_refseq_id;
	return %new_entry;
}


