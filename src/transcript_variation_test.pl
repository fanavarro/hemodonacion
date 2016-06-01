use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $CODON_LENGTH = 3;
my $MET = 'ATG';
my @STOP_CODONS = qw(TAG TAA TGA);

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
	print "Variation: " . $transcript_variation->variation_feature->variation_name . "\n";
	my $tvas = $transcript_variation->get_all_alternate_TranscriptVariationAlleles();
	foreach my $tva ( @{$tvas} ) {
		print "Alelo: " . $tva->dbID . "\n";
		print "codon: " . $tva->display_codon_allele_string  . "\n";
		my $hash = get_sequence_info($tva);
		
		#printHash($hash);
	}
}

sub get_sequence_info{
    my $tva = $_[0];
    my $hash_seq_info = {};
    $hash_seq_info->{'first_met_position'} = -1;
    $hash_seq_info->{'reading_frame'} = ' ';
    $hash_seq_info->{'stop_codon_position'} = ' ';
    my $seq = get_variation_cds_seq($tva);
    #print $seq . "\n";
    $hash_seq_info->{'seq_length'} = length($seq);
    my $first_met_pos = index($seq, $MET);
    if ($first_met_pos != -1){
        my $stop_codon_pos = get_stop_codon_position($seq);
        my $cut_seq = substr($seq, $first_met_pos);
	my $reading_frame = $first_met_pos % $CODON_LENGTH == 0 ? 'Conserved' : 'Lost';
	$hash_seq_info->{'first_met_position'} = $first_met_pos;
        $hash_seq_info->{'stop_codon_position'} = $stop_codon_pos != -1 ? $stop_codon_pos : 'No Stop';
	$hash_seq_info->{'reading_frame'} = $reading_frame;
    }
    return $hash_seq_info;
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
    print "Feature seq: $feature_seq \n";
    substr($seq, $variation_start, $variation_end - $variation_start + 1) = $feature_seq;
    
    return $seq;
}

sub printHash{
	my %hash = %{$_[0]};
	while (my ($key, $value) = each(%hash)) {
		print "$key -> $value\n";
	}
}

sub get_stop_codon_position{
    my $seq = $_[0];
    # get the first met in the sequence
    my $init_pos = index($seq, $MET);
    if ($init_pos != -1){
	# Search stop codons keeping the reading frame
        for (my $i = $init_pos; $i < length($seq); $i = $i + $CODON_LENGTH){
	    my $codon = substr($seq, $i, $CODON_LENGTH);
	    if (exists_in_list($codon, \@STOP_CODONS)){
	        return $i;
            }
        }
    }
    return -1;
}

sub exists_in_list {
    my $check = $_[0];
    my $list  = $_[1];
    foreach my $element ( @{$list} ) {
        if ( $check eq $element ) {
            return 1;
        }
    }
    return 0;
}
