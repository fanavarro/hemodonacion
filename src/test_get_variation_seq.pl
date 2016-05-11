use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $CODON_LENGTH = 3;
my $START_CODON = 'ATG';

my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous',
);
my $trv_adaptor = $registry->get_adaptor( 'homo_sapiens', 'variation', 'transcriptvariation' );
my $transcript_adaptor  = $registry->get_adaptor('human', 'core', 'Transcript');

my $transcript = $transcript_adaptor->fetch_by_stable_id('ENST00000338981');
my @transcripts = ($transcript);
my @so_terms = ('start_lost');
my $trvs = $trv_adaptor->fetch_all_by_Transcripts_SO_terms(\@transcripts, \@so_terms);
foreach my $tv ( @{$trvs} ) {
    my $tvas = $tv->get_all_alternate_TranscriptVariationAlleles();
    foreach my $tva ( @{$tvas} ) {
	#get_variation_coding_seq($tva);
	my $next_met = get_next_met($tva);
	print $next_met . "\n";
    }
}

sub get_next_met{
    my $tva = $_[0];
    my $seq = get_variation_cds_seq($tva);
    my $read_shift = 0;
    for (my $i = 0; $i < length($seq) - ($CODON_LENGTH) + 1; $i++){
	my $codon = substr($seq, $i, $CODON_LENGTH);
        if ($codon eq $START_CODON){
            return $i . " (shift +" . $i % $CODON_LENGTH . ")";
        }
    }
    return "No MET found";
}

sub get_variation_cds_seq{
    my $tva = $_[0];
    # translateable_seq returns the coding part of the transcript
    # (it removes introns and 5' and 3' utr)
    my $seq = $tva->transcript->translateable_seq;
    # Variation position starting at the begining of coding sequence.
    my $variation_start = $tva->transcript_variation->cds_start - 1;
    my $variation_end = $tva->transcript_variation->cds_end - 1;
    # If is a deletion, feature_seq is '-', so we will use '' instead
    # to build the final sequence.
    my $feature_seq = $tva->feature_seq eq "-" ? "" : $tva->feature_seq;

    print $tva->display_codon_allele_string . "\n";
    print $tva->transcript_variation->variation_feature->variation_name . "\t$variation_start-$variation_end\n";
    print "$seq\n";
    substr($seq, $variation_start, $variation_end - $variation_start + 1) = $feature_seq;
    print $seq . "\n";
    return $seq;
}

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
