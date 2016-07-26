package myUtils::SeqUtils;
use strict;
use warnings;
use List::Util qw[min max];
my $CODON_LENGTH = 3;
my $MET = 'ATG';
my @STOP_CODONS = qw(TAG TAA TGA);
my %codon_translation = (
    "TTT" => "F",    "TTC" => "F",    "TTA" => "L",    "TTG" => "L",
    "CTT" => "L",    "CTC" => "L",    "CTA" => "L",    "CTG" => "L",
    "ATT" => "I",    "ATC" => "I",    "ATA" => "I",    "ATG" => "M",
    "GTT" => "V",    "GTC" => "V",    "GTA" => "V",    "GTG" => "V",
    "TCT" => "S",    "TCC" => "S",    "TCA" => "S",    "TCG" => "S",
    "CCT" => "P",    "CCC" => "P",    "CCA" => "P",    "CCG" => "P",
    "ACT" => "T",    "ACC" => "T",    "ACA" => "T",    "ACG" => "T",
    "GCT" => "A",    "GCC" => "A",    "GCA" => "A",    "GCG" => "A",
    "TAT" => "Y",    "TAC" => "Y",    "TAA" => "-",    "TAG" => "-",
    "CAT" => "H",    "CAC" => "H",    "CAA" => "Q",    "CAG" => "Q",
    "AAT" => "N",    "AAC" => "N",    "AAA" => "K",    "AAG" => "K",
    "GAT" => "D",    "GAC" => "D",    "GAA" => "E",    "GAG" => "E",
    "TGT" => "C",    "TGC" => "C",    "TGA" => "-",    "TGG" => "W",
    "CGT" => "R",    "CGC" => "R",    "CGA" => "R",    "CGG" => "R",
    "AGT" => "S",    "AGC" => "S",    "AGA" => "R",    "AGG" => "R",
    "GGT" => "G",    "GGC" => "G",    "GGA" => "G",    "GGG" => "G",
);

# Get the aminoacid sequence from
# dna orf sequence.
# param 0 -> dna orf sequence.
# return aminoacid seq. Stop codon is represented as '-'.
sub get_translation{
    my $orf = shift;
    $orf = uc($orf); # to uppercase
    
    my $aa_seq = "";
    for (my $i = 0; $i < length($orf); $i += $CODON_LENGTH){
        my $codon = substr($orf, $i, $CODON_LENGTH);
        # print "$codon \n";
        $aa_seq = $aa_seq . $codon_translation{$codon} if defined $codon_translation{$codon};
    }
    return $aa_seq;
}

# Return true if both sequences are in frame.
# To figure out this, both sequences are translated
# and checked if one contains the other.
# param0 -> orf of the original sequence.
# param1 -> orf of the mutated sequence.
sub is_in_frame{
    my $original_orf_seq = shift;
    my $mutated_orf_seq = shift;
    my $original_aa = get_translation($original_orf_seq);
    my $mutated_aa = get_translation($mutated_orf_seq);
    
    # if mutation is a deletion is possible that first met in the
    # mutation is not changed; for example AT(GGAGAGTAA)GGATGA...
    # will produce the same met than the original, but three aminoacids
    # after that met are deleted, so we have to check from GGATGA...
    if (length($mutated_aa) < length($original_aa)){
        my $difference = length($original_aa) - length($mutated_aa);
        $original_aa = substr($original_aa, $difference + 1);
    } elsif(length($mutated_aa) > length($original_aa)){
        my $difference = length($mutated_aa) - length($original_aa);
        $mutated_aa = substr($mutated_aa, $difference + 1);
    }
    return (index($original_aa, $mutated_aa) != -1) || (index($mutated_aa, $original_aa) != -1);
}

# Receives a complete transcript sequence
# and a position and returns the orf according
# the first met found after position. It is possible
# that found sequence do not contain stop codon.
# param0 -> Sequence
# param1 -> Position to start the orf search
# return undef if there is not a MET. Return an orf
# if MET and STOP are found in frame. Return a sequence
# starting from MET until the end of sequence if no Stop codon
# found.
sub get_orf{
    my $seq = shift;
    my $pos = shift;
    my $orf = undef;

    # check if position is inside sequence
    if ($pos >= length($seq)){
        return $orf;
    }
    # find met pos
    my $pos_met = index($seq, $MET, $pos);

    # complete de orf until stop codon or the end of seq
    if ($pos_met != -1){
        $orf = '';
        for(my $i = $pos_met; $i < length($seq); $i += $CODON_LENGTH){
            my $codon = substr($seq, $i, $CODON_LENGTH);
            $orf = $orf . $codon;
            if (exists_in_list($codon, \@STOP_CODONS)){
                last;
            }
        }
    }

    return $orf;
}

# Return the position in which the translation
# starts in the cdna sequence.
# param0 -> the cdna sequence (transcript including 5' and 3')
# param1 -> the cds sequence (transcript with only coding sequence)
sub get_translation_start_pos{
    my $cdna = shift;
    my $cds = shift;
    return index($cdna, $cds);
}


# Get the position of the first stop codon following
# the reading frame with start position.
# param 0: string with the sequence
# param 1: start search position
# return the position of the first stop codon found after
# first MET keeping the reading frame.
sub get_stop_codon_position{
    my $seq = $_[0];
    # get the first met in the sequence
    my $init_pos = $_[1];
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

# Checks if an element exists in the given list.
# param 0 -> Element to check.
# param 1 -> List reference.
# returns true if elements exists inside the list and false if not.
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

# Receives the original cdna and cdssequences, the mutated
# affecting initiation codon cdna sequence, and a value indicating
# if the mutation are affecting to the five prime utr. The function
# returns a hash with the following keys:
# 'first_met_position' shows the position of the first met found
# in the mutated sequence, but in coordinates of the original cds
# sequence.
# 'reading_frame' shows if the reading frame is conserved or lost.
# 'stop_codon_position' indicates the position of the first stop codon found
# following the reading frame with 'first_met_pos'-
# 'seq_length' Percentage of the protein that is conserved. For example,
# if a protein has 50 aminoacid and the mutation causes the loss of 25
# aminoacids, this value will be 50%.
sub get_met_mutation_info{
    my $hash_seq_info = {};
    $hash_seq_info->{'first_met_position'} = '';
    $hash_seq_info->{'reading_frame'} = '';
    $hash_seq_info->{'stop_codon_position'} = '';
    $hash_seq_info->{'seq_length'} = '';

    
    my $cdna = $_[0];
    my $cds = $_[1];
    my $mutated_cdna = $_[2];
    my $five_affected = $_[3];
    my $translation_start_pos = get_translation_start_pos($cdna, $cds);

    # We start to count positions from reference
    my $reference;

    # If deletion, we move reference
    if(length($cdna) > length($mutated_cdna) && $five_affected){
        $reference = max($translation_start_pos - (length($cdna) - length($mutated_cdna)), $five_affected);
    } else {
        $reference = $translation_start_pos;
    }
    
    my $first_met_pos = index($mutated_cdna, $MET, $reference);

    # Correction in order to point to the original sequence position.
    # If muation is a insertion, we have to substract values to reference.
    # If deletion, we have to add values to reference only if first met found
    # in mutated seq is different from the natural.
    my $position_correction = 0;
    if (length($mutated_cdna) > length($cdna)){
        $position_correction = -(length($mutated_cdna) - length($cdna));
     } 
     if (length($mutated_cdna) < length($cdna) && !$five_affected){
        $position_correction = (length($cdna) - length($mutated_cdna));
     } 

    
    
    # if a met is found and it is inside cds region, fill result hash
    if ($first_met_pos != -1 &&  max($first_met_pos - $reference + $position_correction, 0) < length($cds)){
        my $stop_codon_pos = get_stop_codon_position($mutated_cdna, $first_met_pos);
        my $mutated_orf = get_orf($mutated_cdna, $first_met_pos);
        #say $mutated_orf;
        #say $cds;
        my $reading_frame = is_in_frame($cds, $mutated_orf) ? 'Conserved' : 'Lost';
        
        # If mutated and original met are different, apply the pos correction
        if ($first_met_pos != $translation_start_pos){
            # Use of max to avoid errors in met duplication cases, where first pos indicated -3.
            $hash_seq_info->{'first_met_position'} = max($first_met_pos - $reference + $position_correction, 0);
        }
        else {
            # Use of max to avoid errors in met duplication cases, where first pos indicated -3.
            $hash_seq_info->{'first_met_position'} = max($first_met_pos - $reference , 0);
        }
        
        $hash_seq_info->{'stop_codon_position'} = $stop_codon_pos != -1 ? $stop_codon_pos - $reference +  $position_correction : 'No Stop';
        $hash_seq_info->{'reading_frame'} = $reading_frame;
        $hash_seq_info->{'seq_length'} = (length($mutated_orf) * 100 / length($cds)) . '%';
        # Print warnings if first met pos calculated is not a MET in the original cds sequence.
        if (substr($cds, $hash_seq_info->{'first_met_position'}, $CODON_LENGTH) ne $MET){
            print ("Error\nfirst met pos = " . $hash_seq_info->{'first_met_position'} . " in $cds\n");
        }
        
        # Print warnings if stop codon pos calculated is not a stop codon in the original cds sequence.
        if ($hash_seq_info->{'stop_codon_position'} ne 'No Stop'){
            if (!exists_in_list(substr($cds, $hash_seq_info->{'stop_codon_position'}, $CODON_LENGTH), \@STOP_CODONS)){
                print ("Error\nstop codon pos = " . $hash_seq_info->{'stop_codon_position'} . " in $cds\n");
            }
        }
    }
    return $hash_seq_info;

}
1;
