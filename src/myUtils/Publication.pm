package myUtils::Publication;
use strict;
use warnings;
use LWP::Simple;
use JSON;

my $DOI_PREFIX='http://doi.org/api/handles/';
my $PUBMED_PREFIX = 'http://www.ncbi.nlm.nih.gov/pubmed/';

# Find a link to publication by its Document Object Id.
# param 0 -> Document Object Id
# returns link to document or undef if error.
sub find_link_by_doi{
    my $doi = $_[0];
    my $url = $DOI_PREFIX . $doi;
    my $response = get $url;
    if (!$response){
        return undef;
    }
    my $hash = decode_json($response);
    if ($hash->{'responseCode'} != 1){
	return undef;
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

# Build a link to the publication in pubmed.
# param 0 -> pmid
# returns url to the publication.
sub find_link_by_pmid{
    my $pmid = $_[0];
    my $url = $PUBMED_PREFIX . '?term=' . $pmid . '[pmid]';
    return $url;
}

# Build a link to the publication in pubmed.
# param 0 -> pmcid
# returns url to the publication.
sub find_link_by_pmcid{
    my $pmcid = $_[0];
    my $url = $PUBMED_PREFIX . '?term=' . $pmcid;
    return $url;
}

1;

