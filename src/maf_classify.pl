use myUtils::CsvManager;
use strict;
use warnings;

#Get intput file from input parameters
my $input = "";
if (scalar @ARGV == 1){
        $input = $ARGV[0];
} else {
        print "Usage: perl maf_classify.pl file.csv";
        exit;
}
my $maf_cut = 0.01;
my @fields = qw(CHROMOSOME GENE_ID GENE_NAME TRANSCRIPT_ID TRANSCRIPT_REFSEQ_ID TRANSCRIPT_BIOTYPE CDS_ERRORS PROTEIN_ID VARIATION_NAME SOURCE TRANSCRIPT_VARIATION_ALLELE_DBID MINOR_ALLELE_FREQUENCY CODON_CHANGE AMINOACID_CHANGE FIRST_MET_POSITION STOP_CODON_POSITION MUTATED_SEQUENCE_LENGTH READING_FRAME_STATUS CONSEQUENCE PHENOTYPE SO_TERM SIFT POLYPHEN PUBLICATIONS);

my %params = (
fields    => \@fields,
csv_separator   => "\t",
in_field_separator    => '-',
file_name => $input,
mode => '<'
);
my %params_out1 = (
fields    => \@fields,
csv_separator   => "\t",
in_field_separator    => '-',
file_name => 'out1.csv',
mode => '>'
);
my %params_out2 = (
fields    => \@fields,
csv_separator   => "\t",
in_field_separator    => '-',
file_name => 'out2.csv',
mode => '>'
);
my $csvManager = myUtils::CsvManager->new (%params);
my $out1 = myUtils::CsvManager->new (%params_out1);
my $out2 = myUtils::CsvManager->new (%params_out2);
$out1->myUtils::CsvManager::writeHeader();
$out2->myUtils::CsvManager::writeHeader();
my $rows = $csvManager->myUtils::CsvManager::readCsv();
$csvManager->myUtils::CsvManager::close();

my @common =[];
my @no_common = [];
my @undefined = [];

my $n = scalar(@{$rows});
my $i = 1;
foreach my $reg (@{$rows}){
	print "$i" . "/" . "$n" . "\n";
	$i = $i + 1;
	if ($reg->{'MINOR_ALLELE_FREQUENCY'} && $reg->{'MINOR_ALLELE_FREQUENCY'} ne '-'){
		if($reg->{'MINOR_ALLELE_FREQUENCY'} < $maf_cut){
			push @no_common, $reg;
			$out1->myUtils::CsvManager::writeEntry(%{$reg});
		}
		else {
			push @common, $reg;
			$out2->myUtils::CsvManager::writeEntry(%{$reg});
		}
	} else{
		push @undefined, $reg;
	}
}
$out1->myUtils::CsvManager::close();
$out2->myUtils::CsvManager::close();
