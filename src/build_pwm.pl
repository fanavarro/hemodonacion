use strict;
use warnings;
use Bio::Matrix::PSM::IO::masta;

my $input = "";
if (scalar @ARGV == 1){
	$input = $ARGV[0];
} else {
	print "Usage: perl build_pwm.pl inputFile format";
	exit;
}
my $fh = \*STDOUT;
my $psmIN =  Bio::Matrix::PSM::IO->new(-format=> 'masta',
                                                 -file  => $input);

my $psmOUT =  Bio::Matrix::PSM::IO->new(-format=> 'masta',
                                                 -fh     => $fh);


my $matrix = $psmIN->next_matrix; #Returns Bio::Matrix::PSM::SiteMatrix object

$psmOUT->write_psm($matrix);
