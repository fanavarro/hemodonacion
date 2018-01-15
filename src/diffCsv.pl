use myUtils::CsvManager;
use strict;
use warnings;
use FileHandle qw( );

# Needed to write to nohup.out in real time
STDOUT->autoflush(1);

##############
my $file_name1;
my $file_name2;
if (scalar @ARGV == 2){
	$file_name1 = $ARGV[0];
	$file_name2 = $ARGV[1];
} else {
	print "Usage: perl diffCsv.pl file1.csv file2.csv\n";
	exit;
}
my $out_file_name1 = substr($file_name1, 0, index ($file_name1, '.')) . '_uniq.csv';
my $out_file_name2 = substr($file_name2, 0, index ($file_name2, '.')) . '_uniq.csv';

my %params1 = (
#fields    => \@fields,
csv_separator   => "\t",
in_field_separator    => '-',
file_name => $file_name1,
mode => '<'
);

my %params2= (
#fields    => \@fields,
csv_separator   => "\t",
in_field_separator    => '-',
file_name => $file_name2,
mode => '<'
);
my $csv1 = myUtils::CsvManager->new (%params1);
my $csv2= myUtils::CsvManager->new (%params2);

my %out_params1 = (
fields    => $csv1->myUtils::CsvManager::fields(),
csv_separator   => "\t",
in_field_separator    => '-',
file_name => $out_file_name1,
mode => '>'
);

my %out_params2= (
fields    => $csv1->myUtils::CsvManager::fields(),
csv_separator   => "\t",
in_field_separator    => '-',
file_name => $out_file_name2,
mode => '>'
);


my $out_csv1 = myUtils::CsvManager->new (%out_params1);
my $out_csv2= myUtils::CsvManager->new (%out_params2);

my $hashref = $csv1->myUtils::CsvManager::diffCsv($csv2);

my @entries1 = @{$hashref->{'this'}};
my @entries2 = @{$hashref->{'other'}};

$out_csv1->myUtils::CsvManager::writeHeader();
while((my $entry = shift @entries1)){
	$out_csv1->myUtils::CsvManager::writeEntry(%{$entry});
}

$out_csv2->myUtils::CsvManager::writeHeader();
while((my $entry = shift @entries2)){
	$out_csv2->myUtils::CsvManager::writeEntry(%{$entry});
}

$csv1->myUtils::CsvManager::close();
$csv2->myUtils::CsvManager::close();
$out_csv1->myUtils::CsvManager::close();
$out_csv2->myUtils::CsvManager::close();


#print "Entradas unicas de posY.csv.csv\n";
#foreach my $reg (@{$hashref->{'this'}}){
#	foreach my $key (@fields){
#		print $reg->{$key} if ($reg->{$key});
#		print "\t";
#	}
#	print "\n";
#}

#print "\nEntradas unicas de soY.csv\n";
#foreach my $reg (@{$hashref->{'other'}}){
#	foreach my $key (@fields){
#		print $reg->{$key} if ($reg->{$key});
#		print "\t";
#	}
#	print "\n";
#}
