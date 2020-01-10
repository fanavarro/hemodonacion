use strict;
use warnings;
use Data::Dumper;
use myUtils::EVAService;

local $Data::Dumper::Terse = 1;
local $Data::Dumper::Indent = 1;
my $eva_service = myUtils::EVAService->instance();
my $taxonomic_code = "mus_musculus";
my $assembly_code = "grcm38";
my @so_terms = ("SO:0002012"); # start lost;
my $pos = "X:73095000-73095100";

my $variant_id1 = "rs578978361";
my $variant_id2 = "rs220009852";

#my $variants = $eva_service->get_variants($taxonomic_code, $assembly_code, $pos, \@so_terms);
#print Dumper($variants);

#my $variant_info = $eva_service->get_variant_info($taxonomic_code, $assembly_code, $variant_id);
#print Dumper($variant_info);

my $highest_population_maf = $eva_service->get_highest_population_maf($taxonomic_code, $assembly_code, $variant_id1);
if (defined $highest_population_maf){
	print "Highest population MAF for ${variant_id1} is ${highest_population_maf}\n";
} else {
	print "Highest population MAF for ${variant_id1} not found\n";
}

$highest_population_maf = $eva_service->get_highest_population_maf($taxonomic_code, $assembly_code, $variant_id2);
if (defined $highest_population_maf){
	print "Highest population MAF for ${variant_id2} is ${highest_population_maf}\n";
} else {
	print "Highest population MAF for ${variant_id2} not found\n";
}
