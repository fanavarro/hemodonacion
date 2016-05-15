use strict;
use warnings;
use List::Util qw[min max];
my $seq1 = 'AGTGAATTATATGAAATGAATT';
my $MET = 'ATG';
my @STOP_CODONS = qw(TAG TAA TGA);

print get_stop_codon_position($seq1);

sub get_stop_codon_position{
    my $seq = $_[0];
    my $init_pos = index($seq, $MET);

    for (my $i = $init_pos; $i < length($seq); $i = $i + 3){
	my $codon = substr($seq, $i, 3);
	if (exists_in_list($codon, \@STOP_CODONS)){
	    return $i;
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
