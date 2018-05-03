package myUtils::Noderer;
use strict;
use warnings;
use myUtils::SeqUtils;
use base 'Class::Singleton';

sub _new_instance {
	my $class = shift;
	my $this  = bless { }, $class;
	
	my $kozak_scores = read_kozak_scores('myUtils/noderer/kozak_scores.txt');

	$this->{kozak_scores} = $kozak_scores;
	return $this;
}

sub read_kozak_scores {
  my $file_path = shift;
  open(my $file_handler, '<', $file_path) or die "Unable to open file, $!";
  
  my %kozak_scores;
  while (my $line = <$file_handler>){
    $line =~ s/^\s+|\s+$//g;  #Trim spaces.
    if (length($line) > 0 && substr($line, 0, 1) ne '#'){
      my @tokens = split(' ', $line);
      my $seq = $tokens[0];
      my $score = $tokens[1];
      $kozak_scores{$seq} = $score;
    }
  }
  close($file_handler) or warn "Unable to close the file handle: $!";
  
  return \%kozak_scores;
}

sub get_kozak_matches {
  my $this = shift;
  my $cdna = shift;
  my $min_efficiency = shift;
  my $hits;
	
	my @starts = ();
	my @ends = ();
	my @widths = ();
	my @scores = ();
	my @init_codon_positions = ();
	
	my $atg_positions = myUtils::SeqUtils::get_met_positions($cdna);
	my $index = 0;
	for (my $i = 0; $i <  scalar (@{$atg_positions}); $i++){
	  # The Kozak context must have 11 nucleotides including ATG (xxxxxxATGxx).
	  my $atg_position = $atg_positions->[$i];
	  if ($atg_position > 6 && $atg_position < length($cdna) - 4){
	    my $kozak_context = substr($cdna, $atg_position - 6, 11);
	    my $efficiency = $this->{kozak_scores}->{$kozak_context};
	    if($efficiency >= $min_efficiency){
	      $starts[$index] = $atg_position - 6;
	      $ends[$index] = $atg_position + 4;
	      $widths[$index] = 11;
	      $scores[$index] = $efficiency;
	      $init_codon_positions[$index] = $atg_position;
	      $index = $index + 1;
	    }
	  }
	}
	
  $hits->{START} = \@starts;
	$hits->{END} = \@ends;
	$hits->{WIDTH} = \@widths;
	$hits->{SCORE} = \@scores;
	$hits->{INIT_CODON_POS} = \@init_codon_positions;
	
	return $hits;
}

1;
