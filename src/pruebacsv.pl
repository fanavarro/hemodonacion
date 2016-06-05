use myUtils::CsvManager;
use strict;
use warnings;

##############
my @fields = qw(a b c);
my %hash1 = (
fields    => \@fields,
csv_separator   => ',',
in_field_separator    => '-',
file_name => "csvtest1.csv",
mode => '<'
);

my %hash2= (
fields    => \@fields,
csv_separator   => ',',
in_field_separator    => '-',
file_name => "csvtest2.csv",
mode => '<'
);

my $csv1 = myUtils::CsvManager->new (%hash1);
my $csv2= myUtils::CsvManager->new (%hash2);
my $hashref = $csv1->myUtils::CsvManager::compareCsv($csv2);
$csv1->myUtils::CsvManager::close();
$csv2->myUtils::CsvManager::close();

#my @list1 = $hashref->{'this'};
#my @list2 = $hashref->{'other'};
print $hashref->{'this'}->[0]->{'a'};
