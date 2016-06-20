use myUtils::CsvManager;
use strict;
use warnings;
use Scalar::Util qw(looks_like_number);

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

my @common;
my @no_common;
my @undefined;

my $n = scalar(@{$rows});
my $i = 1;
foreach my $reg (@{$rows}){
	# print "$i" . "/" . "$n" . "\n";
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

print "MAF >= $maf_cut \n";
print_statistics(get_statistics(\@common));
print "\nMAF < $maf_cut \n";
print_statistics(get_statistics(\@no_common));

# Receives a listref of entries and returns
# average first met position, the count of
# proteins that have a conserved frameshift
# or not, the count of proteins that dont have
# start codon
sub get_statistics{
    my $rows = $_[0];
	my $n = scalar(@{$rows});
	my %statistics;
	$statistics{'average_first_met'} = 0;
	$statistics{'count_frameshift_conserved'} = 0;
	$statistics{'count_frameshift_no_conserved'} = 0;
	$statistics{'count_no_met'} = 0;
        $statistics{'total'} = $n;
	foreach my $reg (@{$rows}){
	    if ($reg->{'FIRST_MET_POSITION'} && looks_like_number($reg->{'FIRST_MET_POSITION'})){
	        $statistics{'average_first_met'} += $reg->{'FIRST_MET_POSITION'};
	    } elsif ($reg->{'FIRST_MET_POSITION'} && $reg->{'FIRST_MET_POSITION'} =~ /^\s*$/ ){
		    $statistics{'count_no_met'} ++;
	    }
            if ($reg->{'READING_FRAME_STATUS'} eq 'Conserved'){
		    $statistics{'count_frameshift_conserved'} ++;
            } elsif($reg->{'READING_FRAME_STATUS'} eq 'Lost'){
		    $statistics{'count_frameshift_no_conserved'} ++;
            }
	}
	$statistics{'average_first_met'} = $statistics{'average_first_met'} / ($n - $statistics{'count_no_met'});
	return \%statistics;
}

sub print_statistics{
    my $statistics = $_[0];
    my $n = $statistics->{'total'};
    print "Average first met position: " . $statistics->{'average_first_met'} . "\n";
	print "Total conserved frameshift: " . $statistics->{'count_frameshift_conserved'} . "/" . "$n (" . $statistics->{'count_frameshift_conserved'}/$n .  ")\n";
	print "Total lost frameshift: " . $statistics->{'count_frameshift_no_conserved'} .  "/" . "$n (" . $statistics->{'count_frameshift_no_conserved'}/$n . ")\n";
	print "Total Peptide without start codon: " . $statistics->{'count_no_met'} .  "/" . "$n (" . $statistics->{'count_no_met'}/$n . ")\n";
}
