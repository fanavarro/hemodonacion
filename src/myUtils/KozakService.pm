package myUtils::KozakService;
use strict;
use warnings;
use LWP::Simple;
use HTML::TableExtract;
use base 'Class::Singleton';

sub _new_instance {
	my $class = shift;
	my $this  = bless { }, $class;
	
	my $browser = LWP::UserAgent->new;

	$this->{BROWSER} = $browser;
	return $this;
}


# Use service at http://atgpr.dbcls.jp/cgi-bin/atgpr.cgi
# in order to get information about kozak sequence in 
# parameter sequence.
# param 0 -> this.
# param 1 -> sequence to get kozak info.
sub get_kozak_info{
	my $this = shift;
	my $sequence = shift;
	my $n_results = shift;

	my $url = 'http://atgpr.dbcls.jp/cgi-bin/atgpr.cgi';
	my $response = $this->{BROWSER}->post( $url, { 'seq' => $sequence, 'number' => $n_results } );
	#foreach my $key (keys %{$response}){
	#	print $key . "\t-> " . $response->{$key} . "\n";
	#}

	my $te = HTML::TableExtract->new( );
	$te->parse($response->decoded_content());

	my $hash_ref_list = [];
	my @header = qw(PREVIOUS_ATGS RELIABILITY FRAME KOZAK_IDENTITY START FINISH ORF_AMINOACID_LENGTH STOP_CODON PROTEIN_SEQUENCE);

	# Initialize empty hash
	my $empty_hash;
	foreach my $key (@header){
		$empty_hash->{$key} = '';
	}

	# Examine first matching table
	my $ts = ($te->tables)[0];

	# If there are results, construct result by retrieving table information
	if (defined ($ts) && defined($ts->rows) && scalar($ts->rows) > 0){
		my @rows = $ts->rows;
		foreach my $row (@rows[1 .. scalar(@rows) - 1]) {
			my $hash_ref;
			for (my $i = 0; $i<scalar(@header); $i++){
				# Delete \n from the value.
				my $trim_value = @{$row}[$i];
				$trim_value =~ tr/\n//d;
				$hash_ref->{$header[$i]} = $trim_value;
				# print $header[$i] . " -> " . $hash_ref->{$header[$i]} . "\n";
			}
			push(@{$hash_ref_list}, $hash_ref);
		}
	} 
	# If ts is not defined, return empty hash
	else {
		push(@{$hash_ref_list}, $empty_hash);
	}
	return $hash_ref_list;

}

1;
