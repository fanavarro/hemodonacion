package myUtils::EVAService;
use strict;
use warnings;
use JSON;
use HTTP::Tiny;
use base 'Class::Singleton';

my $eva_server = "http://www.ebi.ac.uk/eva/webservices/rest/";
my $eva_version = "v1";

sub _new_instance {
	my $class = shift;
	my $this  = bless { }, $class;
	
	my $httpClient = HTTP::Tiny->new();

	$this->{http_client} = $httpClient;
	return $this;
}

# Use EVA service at http://www.ebi.ac.uk/eva/webservices/rest/v1/segments/
# in order to get all variants in a region.
# param 0 -> This.
# param 1 -> Taxonomic code for species.
# param 2 -> Assembly code for species.
# param 3 -> Chromosome interval; e.g. "X:73095000-73095100".
# param 4 -> Sequence Ontology term list to filter variations; e.g. "SO:0002012", which is start lost.
sub get_variants {
	my $this = shift;
	my $taxonomic_code = shift;
	my $assembly_code = shift;
	my $interval = shift;
	my $so_terms = shift;
	my $so_termsParam = _convert_list_to_get_param("annot-ct", $so_terms);
	my $get_url = $eva_server . $eva_version . "/segments/${interval}/variants?species=${taxonomic_code}_${assembly_code}&${so_termsParam}";
	
	return $this->_get_request($get_url);
}

# Use EVA service at http://www.ebi.ac.uk/eva/webservices/rest/v1/variants/
# in order to get all the information about a concrete variant.
# param 0 -> This.
# param 1 -> Taxonomic code for species.
# param 2 -> Assembly code for species.
# param 3 -> Variant ID.
sub get_variant_info {
	my $this = shift;
	my $taxonomic_code = shift;
	my $assembly_code = shift;
	my $variant_id = shift;
	
	my $get_url = $eva_server . $eva_version . "/variants/${variant_id}/info?species=${taxonomic_code}_${assembly_code}";
	return $this->_get_request($get_url);
}

# Use EVA service at http://www.ebi.ac.uk/eva/webservices/rest/v1/variants/
# in order to get the highest population MAF from a concrete variant.
# param 0 -> This.
# param 1 -> Taxonomic code for species.
# param 2 -> Assembly code for species.
# param 3 -> Variant ID.
sub get_highest_population_maf {
	my $this = shift;
	my $taxonomic_code = shift;
	my $assembly_code = shift;
	my $variant_id = shift;
	
	my $highest_population_maf = -1;
	my $variant_info = $this->get_variant_info($taxonomic_code, $assembly_code, $variant_id);
	if (defined($variant_info) && 	defined($variant_info->{response}) && defined($variant_info->{response}->[0]->{result}) && defined($variant_info->{response}->[0]->{result}->[0]->{sourceEntries})){
		my $source_entries = $variant_info->{response}->[0]->{result}->[0]->{sourceEntries};
		for my $source_entry_key (keys (%$source_entries)){
			my $source_entry = $source_entries->{$source_entry_key};
			if(defined($source_entry->{cohortStats})){
				if(defined($source_entry->{cohortStats}->{ALL})){
					my $population_maf = $source_entry->{cohortStats}->{ALL}->{maf};
					if(defined($population_maf) && $population_maf ne '0' && $population_maf > $highest_population_maf){
						$highest_population_maf = $population_maf;
					}
				}
			}
		}
	}
	if ($highest_population_maf == -1) {
		$highest_population_maf = undef;
	}
	return $highest_population_maf;
}

sub _get_request {
	my $this = shift;
	my $url = shift;
	my $response = $this->{http_client}->get($url, {
	  headers => { 'Content-type' => 'application/json' }
	});
	die "Error retrieving information from \"${url}\":\n" . $response->{content} unless $response->{success};
	
	my $hashedResponse = undef;
	if(length $response->{content}) {
		$hashedResponse = decode_json($response->{content});
	}
	return $hashedResponse;
}
sub _convert_list_to_get_param {
	my $varName = shift;
	my $list = shift;
	my $param = "";
	for my $element (@$list){
		$param = "${param}${varName}=${element}&";
	}
	my $lastChar = substr( $param, length($param) -1 , 1 );
	if ($lastChar eq "&"){
		$param = substr $param, 0, length($param) - 1;
	}
	
	return $param;
}

1;
