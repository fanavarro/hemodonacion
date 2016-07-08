package myUtils::KozakUtils;
use strict;
use warnings;
use myUtils::KozakService;
use v5.10.0;


# Param 0 -> CDNA sequence
# Param 1 -> Number of results to retrieve.
# Return list with info about kozak sequences.
sub get_kozak_info{
	my $cdna = shift;
	my $n_results = shift;

	my $kozak_service = myUtils::KozakService->instance();
	my $hash_ref_list = $kozak_service->myUtils::KozakService::get_kozak_info($cdna,$n_results);

	# Sort by start codon found
	my @sorted = sort { $a->{START} <=> $b->{START} } @{$hash_ref_list};
	foreach my $hash_ref ( @sorted ){
		# Minus 1 to start in 0.
		$hash_ref->{START}--;
		$hash_ref->{FINISH}--;
		#say "---------------------------------------------------------------";
		#foreach my $key (keys %{$hash_ref}){
		#	say $key . "\t-> " . $hash_ref->{$key};
		#}
	}
	return \@sorted;
}

1;
