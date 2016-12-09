package myUtils::MatchPWM;
use Statistics::R;
use strict;
use warnings;
use Try::Tiny;
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
#	$hits->{INIT_CODON_POS} = list of init codon positions
sub get_kozak_matches{
	my $this = shift;
	my $seq = shift;
	my $hits;
	try {
		$this->{R}->run('seq = "' . $seq . '"');
	} catch {
		# Por lo visto R se reinicia cuando falla... hay que reinicializar...
		$this->{R}->run(q`library(PWMEnrich)`);
		$this->{R}->run(q`pwm = data.matrix(read.table("myUtils/pwm/kozak_context.pwm"))`);
		try {
			$this->{R}->set('seq', $seq);
		} catch {
			$this->{R}->run(q`library(PWMEnrich)`);
			$this->{R}->run(q`pwm = data.matrix(read.table("myUtils/pwm/kozak_context.pwm"))`);
			return undef;
		};
	};


	
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
	$hits->{INIT_CODON_POS} = $this->{R}->get('init_codon');
	

	# Checkings needed if we only have 1 result or 0 results.
	# If we only have 1 result, this result is a scalar, not a list. We have to create the list with the scalar.
	# If there are not results, the values are integer(0) indicating empty list. In this case, we create the empty list.
	if (!ref($hits->{START})){
		$hits->{START} = $hits->{START} =~ m/\(0\)/ ? [] : [$hits->{START}];
		$hits->{END} =  $hits->{END} =~ m/\(0\)/ ? [] : [$hits->{END}];
		$hits->{WIDTH} =  $hits->{WIDTH} =~ m/\(0\)/ ? [] : [$hits->{WIDTH}];
		$hits->{SCORE} =  $hits->{SCORE} =~ m/\(0\)/ ? [] : [$hits->{SCORE}];
		$hits->{INIT_CODON_POS} =  $hits->{INIT_CODON_POS} =~ m/\(0\)/ ? [] : [$hits->{INIT_CODON_POS}];
	}
	return $hits;
}


1;
