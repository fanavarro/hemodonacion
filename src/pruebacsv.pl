use myUtils::CsvManager;
use strict;
use warnings;

##############
#my @fields = ('CHROMOSOME', 'GENE_ID', 'GENE_NAME', 'TRANSCRIPT_ID', 'VARIATION_NAME', 'MINOR_ALLELE_FREQUENCY', 'CODON_CHANGE', 'AMINOACID_CHANGE', 'CONSEQUENCE', 'SO_TERM', 'SIFT', 'POLYPHEN');
my @fields = qw(CHROMOSOME GENE_ID GENE_NAME TRANSCRIPT_ID VARIATION_NAME MINOR_ALLELE_FREQUENCY CODON_CHANGE AMINOACID_CHANGE CONSEQUENCE SO_TERM SIFT POLYPHEN);
my %hash1 = (
fields    => \@fields,
csv_separator   => ',',
in_field_separator    => '-',
file_name => "prueba.csv",
mode => '<'
);

my %hash2= (
fields    => \@fields,
csv_separator   => ',',
in_field_separator    => '-',
file_name => "prueba2.csv",
mode => '>'
);

my $in = myUtils::CsvManager->new (%hash1);
my $out= myUtils::CsvManager->new (%hash2);
$out->myUtils::CsvManager::writeHeader();
while ((my %entry = $in->myUtils::CsvManager::readEntry())){
	$out->myUtils::CsvManager::writeEntry(%entry);
}
$in->myUtils::CsvManager::close();
$out->myUtils::CsvManager::close();