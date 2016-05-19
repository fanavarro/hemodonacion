package myUtils::CsvManager;
use strict;
use warnings;
use LWP::Simple;
use JSON;

my $DOI_PREFIX='http://doi.org/api/handles/';

# Find a link to publication by its Document Object Id.
# param 0 -> Document Object Id
# returns link to document.
sub find_link_by_doi{
    my $doi = $_[0];
    my $url = $DOI_PREFIX . $doi;
    my $response = get $url;
    if (!$response){
        die("Error getting $url");
    }
    my $hash = decode_json($response);
    if ($hash->{'responseCode'} != 1){
	die("Error getting $url. Code " . $hash->{'responseCode'});
    }
    my @array= @{$hash->{'values'}};
    my $link = '';
    foreach my $value (@array){
	if ($value->{'type'} eq 'URL'){
	    $link = $value->{'data'}->{'value'};
	    last;
	}
    }
    return $link;
}

1;
