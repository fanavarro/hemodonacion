use myUtils::CsvManager;
use strict;
use warnings;
use Scalar::Util qw(looks_like_number);

my $INF = 1e9999;

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

generate_hist_file(\@common, 'hist_common.dat');
generate_hist_file(\@no_common, 'hist_no_common.dat');

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
	$statistics{'frameshift_conserved_average_first_met'} = 0;
	$statistics{'frameshift_no_conserved_average_first_met'} = 0;
	$statistics{'count_no_met'} = 0;
	$statistics{'total'} = $n;
	$statistics{'range'} = '';

	my $min_first_met_pos = $INF;
	my $max_first_met_pos = -$INF;

	my %total_statistics=get_list_statistics(get_first_met_pos_list($rows));
	my %conserved_frameshift_statistics=get_list_statistics(get_first_met_pos_list($rows, 'Conserved'));
	my %lost_frameshift_statistics=get_list_statistics(get_first_met_pos_list($rows, 'Lost'));

	$statistics{'average_first_met'} = $total_statistics{'mean'};
	$statistics{'range'} = $total_statistics{'min'} . " - " . $total_statistics{'max'};
	$statistics{'frameshift_conserved_average_first_met'} = $conserved_frameshift_statistics{'mean'};
	$statistics{'frameshift_conserved_range_first_met'} = $conserved_frameshift_statistics{'min'} . " - " . $conserved_frameshift_statistics{'max'};
	$statistics{'frameshift_no_conserved_average_first_met'} = $lost_frameshift_statistics{'mean'};
	$statistics{'frameshift_no_conserved_range_first_met'} = $lost_frameshift_statistics{'min'} . " - " . $lost_frameshift_statistics{'max'};
	foreach my $reg (@{$rows}){
	    my $first_met_pos = valid_met_pos_value($reg->{'FIRST_MET_POSITION'});
	   
	    if (!defined($first_met_pos)){
		$statistics{'count_no_met'} ++;
	    }
	    if ($reg->{'READING_FRAME_STATUS'} eq 'Conserved'){
		    $statistics{'count_frameshift_conserved'} ++;
	    } elsif($reg->{'READING_FRAME_STATUS'} eq 'Lost'){
		    $statistics{'count_frameshift_no_conserved'} ++;
	    }
	}

	return \%statistics;
}

sub print_statistics{
	my $statistics = $_[0];
	my $n = $statistics->{'total'};
	print "Average first met position: " . $statistics->{'average_first_met'} . "\trange: " . $statistics->{'range'} . "\n";
	print "Total conserved frameshift: " . $statistics->{'count_frameshift_conserved'} . "/" . "$n (" . $statistics->{'count_frameshift_conserved'}/$n .  ")\n";
	print "Average first met position on conserved frameshift: " . $statistics->{'frameshift_conserved_average_first_met'} . "\trange: " . $statistics->{'frameshift_conserved_range_first_met'} . "\n";
	print "Total lost frameshift: " . $statistics->{'count_frameshift_no_conserved'} .  "/" . "$n (" . $statistics->{'count_frameshift_no_conserved'}/$n . ")\n";
	print "Average first met position on no conserved frameshift: " . $statistics->{'frameshift_no_conserved_average_first_met'} . "\trange: " . $statistics->{'frameshift_no_conserved_range_first_met'} . "\n";
	print "Total Peptide without start codon: " . $statistics->{'count_no_met'} .  "/" . "$n (" . $statistics->{'count_no_met'}/$n . ")\n";
}

# Returns an array where each position
# have the count of sequences that have
# the first met in that position.
# param 0 -> Rows from csv file.
sub get_hist_array{
	my $rows = $_[0];
	my @hist;
	foreach my $reg (@{$rows}){
		my $first_met_pos = valid_met_pos_value($reg->{'FIRST_MET_POSITION'});
		if (defined($first_met_pos)){
			if ($hist[$first_met_pos]){
				$hist[$first_met_pos] += 1;
			} else {
				$hist[$first_met_pos] = 1;
			}
		}
	}
	return @hist;
}

# Generate a file that contains a first column
# with first met position, and a second column
# with the count of sequences that have the met
# in that position.
sub generate_hist_file{
	my $rows = $_[0];
	my $filename = $_[1];
	my @hist_array = get_hist_array($rows);
	open(my $fh, ">", $filename);
	for (my $i = 0; $i < scalar(@hist_array); $i++){
		my $count;
		if ($hist_array[$i]){
			$count = $hist_array[$i];
		} else {
			$count = 0;
		}
		print $fh "$i\t$count\n";
	}
	close($fh);
	
}

# Receives a value from FIRST_MET_POS
# column from csv, which could have a
# number or a number and '*'. This function
# returns the number without '*'.
sub valid_met_pos_value{
	my $value = $_[0];
	# If value is a number, return it.
	if(looks_like_number($value)){
		return $value;
	}
	# If value is a number with '*' at the end
	# return only the number.
	if($value =~  '\d+\*'){
		return substr($value, 0, length($value) - 1);
	}
	# if value does not match, return undef
	return undef;
}

# Receives a mapped csv and a frameshift mode.
# Return a list with the column first_met_pos
# filtered by READING_FRAME_STATUS using
# frameshift_mode.
sub get_first_met_pos_list{
	my $rows = $_[0];
	my $frameshift_mode = $_[1];
	my @list;
	foreach my $reg (@{$rows}){
		my $first_met_pos = valid_met_pos_value($reg->{'FIRST_MET_POSITION'});
		my $frameshift = $reg->{'READING_FRAME_STATUS'};
		if (defined($first_met_pos)){
			if(defined($frameshift_mode)){
				if($frameshift eq $frameshift_mode){
					push(@list, $first_met_pos);
				}
			} else{
				push(@list, $first_met_pos);
			}
		}
	}
	return \@list;
}

# Receives a numeric list and returns
# a hash with the following keys: 'max',
# 'min', 'mean' and 'n', accordingly to
# list values.
sub get_list_statistics{
	my @list = @{$_[0]};
	my $max = -$INF;
	my $min = $INF;
	my $mean = 0;
	my $n = scalar(@list);
	my %statistics;

	#print join(',', @list) . "\n\n\n";

	foreach my $item (@list){
		$mean += $item;
		
		if($item > $max){
			$max = $item;
		}

		if($item < $min){
			$min = $item;
		}
	}
	$mean = $mean / $n;
	
	$statistics{'max'} = $max;
	$statistics{'min'} = $min;
	$statistics{'mean'} = $mean;
	$statistics{'n'} = $n;
	return %statistics;
}

