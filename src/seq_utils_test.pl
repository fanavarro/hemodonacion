use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use myUtils::SeqUtils;
use v5.10.0;

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
