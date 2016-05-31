use strict;
use warnings;
use Bio::EnsEMBL::Registry;

# Registry configuration
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous',
);

my $gene_adaptor = $registry->get_adaptor('homo_sapiens', 'core', 'gene');
my $transcript_adaptor = $registry->get_adaptor( 'homo_sapiens', 'core', 'transcript' );
# my $slice_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Slice' );
my $trv_adaptor = $registry->get_adaptor( 'homo_sapiens', 'variation', 'transcriptvariation' );
# my $var_adaptor = $registry->get_adaptor( 'homo_sapiens', 'variation', 'variation' );
# my $vfa = $registry->get_adaptor( 'homo_sapiens', 'variation', 'variationfeature' );

my @genes = @{$gene_adaptor->fetch_all_by_external_name('SERPINC1')};
my $gene = $genes[0];
my $transcripts = $transcript_adaptor->fetch_all_by_Gene($gene);
my $constraint = '(translation_start=1 or translation_end=1)';
my $transcript_variations = $trv_adaptor->fetch_all_by_Transcripts_with_constraint($transcripts, $constraint);

foreach my $transcript_variation (@{$transcript_variations}){
	my $cds_start = $transcript_variation->cds_start;
	my $cds_end = $transcript_variation->cds_end;
	my $pep_start = $transcript_variation->translation_start;
	my $pep_end = $transcript_variation->translation_end;
	if (!defined($cds_start)){
		$cds_start = 'N';
	}
	if (!defined($cds_end)){
		$cds_end = 'N';
	}
	if (!defined($pep_start)){
		$pep_start = 'N';
	}
	if (!defined($pep_end)){
		$pep_end = 'N';
	}
	print $transcript_variation->variation_feature->variation->name . "\t" . $cds_start . " - " . $cds_end . "\t" . $pep_start . " - " . $pep_end . "\n";
}

sub position_in_transcript{
	my $tv = $_[0];
	
}
