package myUtils::MatchPWM;
use Statistics::R;
use strict;
use warnings;
use base 'Class::Singleton';

sub _new_instance {
	my $class = shift;
	my $this  = bless { }, $class;
	
	my $R = Statistics::R->new();

	$this->{R} = $R;
	my $out = $this->{R}->run(q`library(PWMEnrich)`);
	$out = $this->{R}->run(q`pwm = data.matrix(read.table("myUtils/pwm/kozak_context.pwm"))`);
	return $this;
}

# Perform a matching between kozak PWM and a sequence.
# param 0 -> Sequence where kozak sequences are going to be found.
# return a hash of lists:
#	$hits->{START} = list of the beginings of matches.
#	$hits->{END} = list of the endings of matches.
#	$hits->{WIDTH} = list of the widths of matches.
#	$hits->{SCORE} = list of the scores of matches.
#	$hits->{INIT_CODON} = list of init codon positions
sub get_kozak_matches{
	my $this = shift;
	my $seq = shift;
	my $hits;
	
	$this->{R}->set('seq', $seq);
	$this->{R}->run(q`hits = matchPWM(pwm, seq, with.score = T)`);
	$this->{R}->run(q`start = start(hits) - 1`);
	$this->{R}->run(q`end = end(hits) - 1`);
	$this->{R}->run(q`width = width(hits)`);
	$this->{R}->run(q`score = mcols(hits)$score`);
	$this->{R}->run(q`init_codon = (start(hits) + ((width(hits)-3)/2))-1`);

	
	$hits->{START} = $this->{R}->get('start');
	$hits->{END} = $this->{R}->get('end');
	$hits->{WIDTH} = $this->{R}->get('width');
	$hits->{SCORE} = $this->{R}->get('score');
	$hits->{INIT_CODON} = $this->{R}->get('init_codon');
	return $hits;
}


1;
