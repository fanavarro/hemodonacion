use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use myUtils::SeqUtils;
use v5.10.0;
use List::Util qw[min max];

my $MET = 'ATG';
my $seq1 = 'ATGCTCATTATGCTCCAAATCTAA';
my $seq2 = 'ATGCTCCAAATCTAA';
my $seq3 = 'ATGCGGCAGCATGAATAA';
my $seq4 = 'ATGAATAACTAA';



# say myUtils::SeqUtils::get_translation($seq1);
# say myUtils::SeqUtils::get_translation($seq2);

# say myUtils::SeqUtils::get_orf($seq1, 0);

if (myUtils::SeqUtils::is_in_frame($seq1, $seq2)){
    say "TEST PASS";
} else {
    say "TEST FAIL";
}

if (myUtils::SeqUtils::is_in_frame($seq2, $seq1)){
    say "TEST PASS";
} else {
    say "TEST FAIL";
}

if (myUtils::SeqUtils::is_in_frame($seq3, $seq4)){
    say "TEST FAIL";
} else {
    say "TEST PASS";
}

if (myUtils::SeqUtils::is_in_frame($seq4, $seq3)){
    say "TEST FAIL";
} else {
    say "TEST PASS";
}

if ($seq1 eq myUtils::SeqUtils::get_orf($seq1, 0)){
    say "TEST PASS";
} else{
    say "TEST FAIL";
}

if ($seq2 eq myUtils::SeqUtils::get_orf($seq1, 1)){
    say "TEST PASS";
} else{
    say "TEST FAIL";
}

if (!defined(myUtils::SeqUtils::get_orf($seq1, 90))){
    say "TEST PASS";
} else {
    say "TEST FAIL";
}

if (!defined(myUtils::SeqUtils::get_orf('AAAAAAAAAAAAA', 90))){
    say "TEST PASS";
} else {
    say "TEST FAIL";
}

if ($seq3 eq myUtils::SeqUtils::get_orf($seq3, 0)){
    say "TEST PASS";
} else{
    say "TEST FAIL";
}

if ('ATGAATAA' eq myUtils::SeqUtils::get_orf($seq3, 1)){
    say "TEST PASS";
} else{
    say "TEST FAIL";
}

print "\n\n\n\n\n\n";
my $mutated_cdna = 'GATGATTCCGGCCCCCTAAGATAA';
my $cdna = 'GATTATATGATTCCGGCCCCCTAAGATAA';
my $cds = 'ATGATTCCGGCCCCCTAA';
my $hashref = get_sequence_info_dbsnp2($cdna, $cds, $mutated_cdna);
if ($hashref->{first_met_position} == 0 && $hashref->{stop_codon_position} == 15 && $hashref->{reading_frame} eq 'Conserved'){
    say "TEST PASS";
} else {
    say "TEST FAIL";
}

$mutated_cdna = 'GATTATACGATTCCGGCCCCCTAAGATAA';
$cdna = 'GATTATATGATTCCGGCCCCCTAAGATAA';
$cds = 'ATGATTCCGGCCCCCTAA';
$hashref = get_sequence_info_dbsnp2($cdna, $cds, $mutated_cdna);
if ($hashref->{first_met_position} eq '' && $hashref->{stop_codon_position} eq '' && $hashref->{reading_frame} eq ''){
    say "TEST PASS";
} else {
    say "TEST FAIL";
}

$mutated_cdna = 'GATTATACGATTCCGGCCATGCCCTAAGATAA';
$cdna = 'GATTATATGATTCCGGCCATGCCCTAAGATAA';
$cds = 'ATGATTCCGGCCATGCCCTAA';
$hashref = get_sequence_info_dbsnp2($cdna, $cds, $mutated_cdna);
if ($hashref->{first_met_position} == 12 && $hashref->{stop_codon_position} == 18 && $hashref->{reading_frame} eq 'Conserved'){
    say "TEST PASS";
} else {
    say "TEST FAIL";
}

$mutated_cdna = 'GATTATATTATTCATGCCATCCCCTAAGATAA';
$cdna = 'GATTATATGATTCATGCCATCCCCTAAGATAA';
$cds = 'ATGATTCATGCCATCCCCTAA';
$hashref = get_sequence_info_dbsnp2($cdna, $cds, $mutated_cdna);
if ($hashref->{first_met_position} == 7 && $hashref->{stop_codon_position} eq 'No Stop' && $hashref->{reading_frame} eq 'Lost'){
    say "TEST PASS";
} else {
    say "TEST FAIL";
}

$mutated_cdna = 'TGATTCATGCCATCCCCTAAGATAA';
$cdna = 'GATTATATGATTCATGCCATCCCCTAAGATAA';
$cds = 'ATGATTCATGCCATCCCCTAA';
$hashref = get_sequence_info_dbsnp2($cdna, $cds, $mutated_cdna);
if ($hashref->{first_met_position} == 6 && $hashref->{stop_codon_position} eq 'No Stop' && $hashref->{reading_frame} eq 'Lost'){
    say "TEST PASS";
} else {
    say "TEST FAIL";
}

$mutated_cdna = 'GATTATAAAAATGATTCATGCCATCCCCTAAGATAA';
$cdna = 'GATTATATGATTCATGCCATCCCCTAAGATAA';
$cds = 'ATGATTCATGCCATCCCCTAA';
$hashref = get_sequence_info_dbsnp2($cdna, $cds, $mutated_cdna);
if ($hashref->{first_met_position} == 4 && $hashref->{stop_codon_position} == 22 && $hashref->{reading_frame} eq 'Conserved'){
    say "TEST PASS";
} else {
    say "TEST FAIL";
}

sub get_sequence_info_dbsnp2{
    my $hash_seq_info = {};
    $hash_seq_info->{'first_met_position'} = '';
    $hash_seq_info->{'reading_frame'} = '';
    $hash_seq_info->{'stop_codon_position'} = '';
    $hash_seq_info->{'seq_length'} = '';

    my $mutated_cdna = 'GATGATTCCGGCCCCCTAAGATAA';
    my $cdna = 'GATTATATGATTCCGGCCCCCTAAGATAA';
    my $cds = 'ATGATTCCGGCCCCCTAA';
    my $mutated_cdna = $_[2];
    my $cdna = $_[0];
    my $cds = $_[1];
    my $translation_start_pos = myUtils::SeqUtils::get_translation_start_pos($cdna, $cds);

    # We start to count positions from reference
    my $reference;

    # If deletion, we move reference
    if(length($cdna) > length($mutated_cdna)){
        $reference = max($translation_start_pos - (length($cdna) - length($mutated_cdna)), 0);
    } else {
        $reference = $translation_start_pos;
    }



    
    my $first_met_pos = index($mutated_cdna, $MET, $reference);
    if ($first_met_pos != -1){
        my $stop_codon_pos = myUtils::SeqUtils::get_stop_codon_position($mutated_cdna, $first_met_pos);
        my $mutated_orf = myUtils::SeqUtils::get_orf($mutated_cdna, $first_met_pos);
        my $reading_frame = myUtils::SeqUtils::is_in_frame($mutated_orf, $cds) ? 'Conserved' : 'Lost';
        $hash_seq_info->{'first_met_position'} = $first_met_pos - $reference;
        $hash_seq_info->{'stop_codon_position'} = $stop_codon_pos != -1 ? $stop_codon_pos - $reference : 'No Stop';
        $hash_seq_info->{'reading_frame'} = $reading_frame;
        $hash_seq_info->{'seq_length'} = (length($mutated_orf) * 100 / length($cds)) . '%';

        # if an ORF exists into the seq, seq file is generated.
        # if($stop_codon_pos != -1){
        #    my $orf = substr($seq, $first_met_pos, $stop_codon_pos - $first_met_pos + 3);
        #    generate_variation_seq_files($tva, $orf);
        # }
    #say "reference -> " . $reference;
    #say "first met pos -> " . $hash_seq_info->{'first_met_position'};
    #say "stop codon pos -> " . $hash_seq_info->{'stop_codon_position'};
    #say "reading frame -> " . $hash_seq_info->{'reading_frame'};
    #say "seq length -> " . $hash_seq_info->{'seq_length'};
    }
    return $hash_seq_info;

}


