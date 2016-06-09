use myUtils::CsvManager;
use strict;
use warnings;

my %params = (
csv_separator   => "\t",
in_field_separator    => '-',
file_name => 'by_pos_final.csv',
mode => '<'
);
my $csvManager = myUtils::CsvManager->new (%params);
my $rows = $csvManager->myUtils::CsvManager::readCsv();
$csvManager->myUtils::CsvManager::close();
#foreach my $reg (@{$rows}){
#	foreach my $key (@{$csvManager->{fields}}){
#		print $reg->{$key} if ($reg->{$key});
#		print "\t";
#	}
#	print "\n";
#}
