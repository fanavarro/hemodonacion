use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous',
);
my $trv_adaptor = $registry->get_adaptor( 'homo_sapiens', 'variation', 'transcriptvariation' );
my $transcript_adaptor  = $registry->get_adaptor('human', 'core', 'Transcript');

my $transcript = $transcript_adaptor->fetch_by_stable_id('ENST00000444521');
my @transcripts = ($transcript);
my @so_terms = ('start_lost');
my $trvs = $trv_adaptor->fetch_all_by_Transcripts_SO_terms(\@transcripts, \@so_terms);
foreach my $tv ( @{$trvs} ) {
    my $tvas = $tv->get_all_alternate_TranscriptVariationAlleles();
    foreach my $tva ( @{$tvas} ) {
	print $tva;
	get_variation_seq($tva);
    }
}

get_variation_seq();
# param 0 -> TransciptVariationAllele object
# return -> Sequence of the variation including 5' and 3' regions.
sub get_variation_seq{
    my $tva = $_[0];
    # translateable_seq returns the coding part of the transcript
    # (it removes introns and 5' and 3' utr)
    # my $seq = $tva->transcript->translateable_seq;
    # seq contains 5' and 3' regions.
    my $seq = $tva->transcript->seq->seq;
    my $variation_start = $tva->transcript_variation->cdna_start - 1;
    my $variation_end = $tva->transcript_variation->cdna_end - 1;
    # If is a deletion, feature_seq is '-', so we will use '' instead
    # to build the final sequence.
    my $feature_seq = $tva->feature_seq eq "-" ? "" : $tva->feature_seq;
    substr($seq, $variation_start, $variation_end - $variation_start + 1) = $feature_seq;
    
    return $seq;
}
