package myUtils::CsvManager;
use strict;
use warnings;
use LWP::Simple;
use Data::Dumper;
use JSON;

my $DOI_PREFIX='http://doi.org/api/handles/';

# Find a link to publication by its Document Object Id.
# param 0 -> Document Object Id
# returns link to document.
sub find_by_doi{
    my $doi = $_[0];
    my $url = $DOI_PREFIX . $doi;
    print $url . "\n";
    my $response = get $url;
    if (!$response){
        die("Error getting $url");
    }
    my $hash = decode_json($response);
    foreach my %value ($hash->{'values'}){
	print %value{'type'};
    }
    return Dumper(decode_json($response));
}

my $response = find_by_doi('10.1016/j.humpath.2008.02.014');
#print $response . "\n";
#1;
