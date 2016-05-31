use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $transcript_id = '';
if (scalar @ARGV == 1){
	$transcript_id = $ARGV[0];
} else {
	print "Usage: perl vartest2.pl transcriptId";
	exit;
}

# Registry configuration
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous',
);

my $transcript_adaptor = $registry->get_adaptor( 'homo_sapiens', 'core', 'transcript' );
# my $slice_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Slice' );
my $trv_adaptor = $registry->get_adaptor( 'homo_sapiens', 'variation', 'transcriptvariation' );
# my $var_adaptor = $registry->get_adaptor( 'homo_sapiens', 'variation', 'variation' );
# my $vfa = $registry->get_adaptor( 'homo_sapiens', 'variation', 'variationfeature' );

# $var = $var_adaptor->fetch_by_name('CM034685');

my @so_terms = ('start_lost');
my $transcript = $transcript_adaptor->fetch_by_stable_id($transcript_id);

# my $trvs = $trv_adaptor->fetch_all_by_Transcripts_SO_terms([$transcript], \@so_terms);
my $trvs = $trv_adaptor->fetch_all_by_Transcripts([$transcript], \@so_terms);

foreach my $trv (@{$trvs}){
	#print $trv->transcript_stable_id() . "\n";
	print $trv->variation_feature->variation->name . "\n";
}

#ENST00000487183
#CM034685
#ENST00000367698
#CM034685
#ENST00000494024
#CM034685
#ENST00000516422
#CM034685
#LRG_577t1
#CM034685
#ENST00000617423
#CM034685

