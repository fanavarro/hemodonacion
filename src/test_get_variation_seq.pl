use strict;
use warnings;
use myUtils::KozakUtils;
use Bio::EnsEMBL::Registry;
use v5.10.0;

my $CODON_LENGTH = 3;
my $START_CODON = 'ATG';
my $MAX_KOZAK_RESULTS = 50000;

my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous',
);
my $trv_adaptor = $registry->get_adaptor( 'homo_sapiens', 'variation', 'transcriptvariation' );
my $transcript_adaptor  = $registry->get_adaptor('human', 'core', 'Transcript');

#my $transcript = $transcript_adaptor->fetch_by_stable_id('ENST00000366603');
my $transcript = $transcript_adaptor->fetch_by_stable_id('ENST00000313695');
my @transcripts = ($transcript);
my $constraint = "(translation_start=1 or translation_end=1)";
my $trvs = $trv_adaptor->fetch_all_by_Transcripts_with_constraint(\@transcripts, $constraint);
foreach my $tv ( @{$trvs} ) {
    my $tvas = $tv->get_all_alternate_TranscriptVariationAlleles();
    foreach my $tva ( @{$tvas} ) {
	print $tv->variation_feature->variation->stable_id . "\n" . get_variation_cdna_seq($tva) . "\n\n";
	#get_variation_cdna_seq($tva);
        #get_variation_cds_seq_upstream($tva);
	#my $next_met = get_first_met_info($tva);
	#print $next_met . "\n";
	#get_variation_seq($tva);
        #get_kozak_info2($tva);
	print "---------------------------\n";
    }
}

sub get_variation_cdna_seq{
    my $tva = $_[0];
    # seq contains 5' and 3' regions.
    my $seq = $tva->transcript->seq->seq;
    if (!defined($tva->transcript_variation->cdna_start) || !defined($tva->transcript_variation->cdna_end)){
        print "ERROR CDNA variation without start or end " . $tva->transcript_variation->variation_feature->variation_name . " " . $tva->transcript->display_id . "\n";
        return undef;
    }
    # Variation position counting utr regions.
    my $variation_start = $tva->transcript_variation->cdna_start - 1;
    my $variation_end = $tva->transcript_variation->cdna_end - 1;
    # If is a deletion, feature_seq is '-', so we will use '' instead
    # to build the final sequence.
    my $feature_seq = $tva->feature_seq eq "-" ? "" : $tva->feature_seq;

    if ( $feature_seq =~ m/[^ATGCNatgcn]/ ){
        print "Feature Seq not available -> " . $feature_seq . "\n";
        return undef;
    }
    substr($seq, $variation_start, $variation_end - $variation_start + 1) = $feature_seq;
    
    
    return $seq;
}

