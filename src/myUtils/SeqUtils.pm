package myUtils::SeqUtils;
use strict;
use warnings;
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
            if (exists_in_list($codon, \@STOP_CODONS)){
                last;
            }
            $orf = $orf . $codon;
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
1;
