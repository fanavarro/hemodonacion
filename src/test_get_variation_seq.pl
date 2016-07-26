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
my $transcript = $transcript_adaptor->fetch_by_stable_id('ENST00000373116');
my @transcripts = ($transcript);
my $constraint = "(translation_start=1 or translation_end=1)";
my $trvs = $trv_adaptor->fetch_all_by_Transcripts_with_constraint(\@transcripts, $constraint);
foreach my $tv ( @{$trvs} ) {
    my $tvas = $tv->get_all_alternate_TranscriptVariationAlleles();
    foreach my $tva ( @{$tvas} ) {
	get_variation_cdna_seq($tva);
	#get_variation_cdna_seq($tva);
        #get_variation_cds_seq_upstream($tva);
	#my $next_met = get_first_met_info($tva);
	#print $next_met . "\n";
	#get_variation_seq($tva);
        #get_kozak_info2($tva);
	print "---------------------------\n";
    }
}

sub get_first_met_info{
    my $tva = $_[0];
    my $seq = get_variation_cds_seq($tva);
    my $first_met_pos = index($seq, 'ATG');
    if ($first_met_pos != -1){
        my $cut_seq = substr($seq, $first_met_pos);
        return $first_met_pos . " (shift +" . length($cut_seq) % $CODON_LENGTH . ")";
    } else{
        return "No MET found";
    }
}

# Obtain info about the first kozak
# sequence on a mutated sequence.
# param 0 -> tva
sub get_kozak_info{
    my $tva = $_[0];
    my $source = $tva->variation_feature->variation->source;
    my $hash_kozak_info = {};
    if (defined($source)){
        if ($source->name ne 'dbSNP'){
             $hash_kozak_info->{'FRAMESHIFT'} = '';
             $hash_kozak_info->{'PREVIOUS_ATGS'} = '';
             $hash_kozak_info->{'RELIABILITY'} = '';
             $hash_kozak_info->{'KOZAK_IDENTITY'} = '';
             $hash_kozak_info->{'START'} = '';
             $hash_kozak_info->{'FINISH'} = '';
             $hash_kozak_info->{'ORF_AMINOACID_LENGTH'} = '';
             $hash_kozak_info->{'STOP_CODON'} = '';
             $hash_kozak_info->{'PROTEIN_SEQUENCE'} = '';
             return $hash_kozak_info;
        }
    }
    my $original_seq = $tva->transcript->seq->seq;
    my $mutated_seq = get_variation_cdna_seq($tva);
    
    # Singleton instance for kozak_service to retrieve
    # kozak information at http://atgpr.dbcls.jp/
    my $kozak_service = myUtils::KozakService->instance();

    my $original_kozak =  $kozak_service->myUtils::KozakService::get_kozak_info($original_seq);
    my $mutated_kozak =  $kozak_service->myUtils::KozakService::get_kozak_info($mutated_seq);
    my $first_original_kozak = $original_kozak->[0];
    my $first_mutated_kozak = $mutated_kozak->[0];
    
    # I use ne operator to avoid errors when kozak service returns empty string.
    if ($first_original_kozak->{'FRAME'} ne $first_mutated_kozak->{'FRAME'}){
        $hash_kozak_info->{'FRAMESHIFT'} = 'Lost';
    } else {
        $hash_kozak_info->{'FRAMESHIFT'} = 'Conserved';
    }
    
     $hash_kozak_info->{'PREVIOUS_ATGS'} = $first_mutated_kozak->{'PREVIOUS_ATGS'};
     $hash_kozak_info->{'RELIABILITY'} = $first_mutated_kozak->{'RELIABILITY'};
     $hash_kozak_info->{'KOZAK_IDENTITY'} = $first_mutated_kozak->{'KOZAK_IDENTITY'};
     # We perform the substraction to get the position in cds sequence
     $hash_kozak_info->{'START'} = $first_mutated_kozak->{'START'} - $first_original_kozak->{'START'};
     # We add +1 to point at the first base of the stop codon
     $hash_kozak_info->{'FINISH'} = $first_mutated_kozak->{'FINISH'} - $first_original_kozak->{'START'} + 1;
     $hash_kozak_info->{'ORF_AMINOACID_LENGTH'} = $first_mutated_kozak->{'ORF_AMINOACID_LENGTH'};
     $hash_kozak_info->{'STOP_CODON'} = $first_mutated_kozak->{'STOP_CODON'};
     $hash_kozak_info->{'PROTEIN_SEQUENCE'} = $first_mutated_kozak->{'PROTEIN_SEQUENCE'};

     say "Inicio kozak original = " . $first_original_kozak->{'START'};
     say "Inicio kozak mutada = " . $first_mutated_kozak->{'START'};
     say "Inicio kozak mutada con respecto a cds = " . $hash_kozak_info->{'START'};
     #get_variation_cds_seq($tva);
     #print $first_mutated_kozak->{'START'} . " - " . $first_original_kozak->{'START'} . " = " . $hash_kozak_info->{'START'} . "\n";

     return $hash_kozak_info;
}
sub get_variation_cds_seq{
    my $tva = $_[0];
    # translateable_seq returns the coding part of the transcript
    # (it removes introns and 5' and 3' utr)
    my $seq = $tva->transcript->translateable_seq;
    # Variation position starting at the begining of coding sequence.
    my $variation_start = $tva->transcript_variation->cds_start - 1;
    my $variation_end = $tva->transcript_variation->cds_end - 1;
    if (!defined($tva->transcript_variation->cds_start)){
        return get_variation_cds_seq_upstream($tva);
    }
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

sub get_variation_cds_seq_upstream{
    my $tva = $_[0];
    my $variation_start = $tva->transcript_variation->cdna_start - $tva->transcript_variation->cds_end - 1;
    my $variation_end = $tva->transcript_variation->cdna_end - 1;

    my $cdna_mutated_seq = get_variation_cdna_seq($tva);
    my $three_prime_utr_seq = undef;
    my $cds_end = undef;
    my $cds_start = $variation_start;
    if (defined($tva->transcript->three_prime_utr)){
        $three_prime_utr_seq = $tva->transcript->three_prime_utr->seq;
        $cds_end = index($cdna_mutated_seq, $three_prime_utr_seq);
    } else {
        $cds_end = length($cdna_mutated_seq) - 1;
    }
    my $final_seq = substr($cdna_mutated_seq, $cds_start, $cds_end - $cds_start);
    say $final_seq;
    return $final_seq
}
# param 0 -> TransciptVariationAllele object
# return -> Sequence of the variation including 5' and 3' regions.
sub get_variation_cdna_seq{
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
    
    print $tva->display_codon_allele_string . "\n";
    print $tva->transcript_variation->variation_feature->variation_name . "\t$variation_start-$variation_end\n";
    print "$seq\n";
    substr($seq, $variation_start, $variation_end - $variation_start + 1) = $feature_seq;
    print $seq . "\n";
    
    return $seq;
}
# Get the first kozak sequence in the
# original transcript and in the mutated
# transcript. This information is used
# to determine if the mutated transcript
# has lost the frameshift.
# param 0 -> TranscriptVariationAllele object
# return -> Hash with kozak sequence info in
# the mutated transcript.
sub get_kozak_info2{
    my $tva = $_[0];
    my $source = $tva->variation_feature->variation->source;
    my $hash_kozak_info = {};
    if (defined($source)){
        if ($source->name ne 'dbSNP'){
             $hash_kozak_info->{'FRAMESHIFT'} = '';
             $hash_kozak_info->{'PREVIOUS_ATGS'} = '';
             $hash_kozak_info->{'RELIABILITY'} = '';
             $hash_kozak_info->{'KOZAK_IDENTITY'} = '';
             $hash_kozak_info->{'START'} = '';
             $hash_kozak_info->{'FINISH'} = '';
             $hash_kozak_info->{'ORF_AMINOACID_LENGTH'} = '';
             $hash_kozak_info->{'STOP_CODON'} = '';
             $hash_kozak_info->{'PROTEIN_SEQUENCE'} = '';
             return $hash_kozak_info;
        }
    }
    my $original_cdna_seq = $tva->transcript->seq->seq;
    my $original_cds_seq = $tva->transcript->translateable_seq;
    my $mutated_cdna_seq = get_variation_cdna_seq($tva);
    my $mutated_cds_seq = get_variation_cds_seq($tva);
    # Find the position in which coding region starts at cdna sequence.
    my $default_kozak_pos = index($original_cdna_seq, $original_cds_seq);
    say "Inicio kozak = " . $default_kozak_pos;
    my $mutated_kozaks =  myUtils::KozakUtils::get_kozak_info($mutated_cdna_seq, $MAX_KOZAK_RESULTS);
    my $original_kozaks =  myUtils::KozakUtils::get_kozak_info($original_cdna_seq, $MAX_KOZAK_RESULTS);
    


    # Look for the natural kozak
    my $natural_kozak = ();
    foreach my $original_kozak (@{$original_kozaks}){
        if ($original_kozak->{'START'} == $default_kozak_pos){
            $natural_kozak = $original_kozak;
            last;
        }
    }

    # Look for the first kozak after default kozak in mutated.
    my $first_mutated_kozak = ();
    foreach my $mutated_kozak (@{$mutated_kozaks}){
        if ($mutated_kozak->{'START'} >= $default_kozak_pos){
        #if ($mutated_kozak->{'START'} == index($mutated_cdna_seq, $mutated_cds_seq)){
            $first_mutated_kozak = $mutated_kozak;
            last;
        }
    }

    if (!defined($natural_kozak->{'START'}) || !defined($first_mutated_kozak->{'START'})){
        print "Error kozak" . $tva->transcript->display_id . " " . $tva->transcript_variation->variation_feature->variation_name . "\n";
        return $hash_kozak_info;
    }
    # I use ne operator to avoid errors when kozak service returns empty string.
    if ($natural_kozak->{'FRAME'} ne $first_mutated_kozak->{'FRAME'}){
        $hash_kozak_info->{'FRAMESHIFT'} = 'Lost';
    } else {
        $hash_kozak_info->{'FRAMESHIFT'} = 'Conserved';
    }
    
     $hash_kozak_info->{'PREVIOUS_ATGS'} = $first_mutated_kozak->{'PREVIOUS_ATGS'};
     $hash_kozak_info->{'RELIABILITY'} = $first_mutated_kozak->{'RELIABILITY'};
     $hash_kozak_info->{'KOZAK_IDENTITY'} = $first_mutated_kozak->{'KOZAK_IDENTITY'};
     # We perform the substraction to get the position in cds sequence
     $hash_kozak_info->{'START'} = $first_mutated_kozak->{'START'} - $natural_kozak->{'START'};
     # We add +1 to point at the first base of the stop codon
     $hash_kozak_info->{'FINISH'} = $first_mutated_kozak->{'FINISH'} - $natural_kozak->{'START'} + 1;
     $hash_kozak_info->{'ORF_AMINOACID_LENGTH'} = $first_mutated_kozak->{'ORF_AMINOACID_LENGTH'};
     $hash_kozak_info->{'STOP_CODON'} = $first_mutated_kozak->{'STOP_CODON'};
     $hash_kozak_info->{'PROTEIN_SEQUENCE'} = $first_mutated_kozak->{'PROTEIN_SEQUENCE'};
     say "Inicio kozak original = " . $natural_kozak->{'START'};
     say "Inicio kozak mutada = " . $first_mutated_kozak->{'START'};
     say "Inicio kozak mutada con respecto a cds = " . $hash_kozak_info->{'START'};
     return $hash_kozak_info;
}
