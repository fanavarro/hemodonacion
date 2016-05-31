use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $varname = '';
if (scalar @ARGV == 1){
	$varname = $ARGV[0];
} else {
	print "Usage: perl vartest.pl varName";
	exit;
}

# Registry configuration
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous',
);

# my $transcript_adaptor = $registry->get_adaptor( 'homo_sapiens', 'core', 'transcript' );
# my $slice_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Slice' );
# my $trv_adaptor = $registry->get_adaptor( 'homo_sapiens', 'variation', 'transcriptvariation' );
my $var_adaptor = $registry->get_adaptor( 'homo_sapiens', 'variation', 'variation' );
my $vfa = $registry->get_adaptor( 'homo_sapiens', 'variation', 'variationfeature' );

# $var = $var_adaptor->fetch_by_name('CM034685');
my $var = $var_adaptor->fetch_by_name($varname);

if ($var){
	my $name = $var->name;
	my $stable_id = $var->stable_id;
	print "Name: $name\n";
	print "Stable Id: $stable_id\n";
	my $sources = $var->get_all_synonym_sources();
	foreach my $src (@$sources) {
		my @synonyms = $var->get_all_synonyms($src);
		print "$src: @synonyms\n";
	}

	my $variation_features = $vfa->fetch_all_by_Variation($var);
	print scalar @{$variation_features} . " variation features \n";
	my @transcript_variations = ();
	foreach my $variation_feature (@{$variation_features}){
		print feature2string($variation_feature) . "\n";
		push(@transcript_variations, @{$variation_feature->get_all_TranscriptVariations()});
	}

	foreach my $transcript_variation (@transcript_variations){
		print $transcript_variation->transcript_stable_id . "\n";
		print $transcript_variation->variation_feature->variation->name . "\n";
	}


} else {
	print "$varname no encontrado";
}

sub feature2string
{
    my $feature = shift;

    my $display_id  = $feature->display_id();
    my $name = $feature->name();
    my $source_name = $feature->source_name();
    my $so_term = $feature->class_SO_term();

    return sprintf( "%s:%s:%s:%s",
        $display_id, $name, $source_name, $so_term );
}
