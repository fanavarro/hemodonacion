package myUtils::Noderer;
use strict;
use warnings;
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

sub get_kozak_info {
  my $cdna = shift;
  my $hits;
  
  $hits->{START} = '';
	$hits->{END} = '';
	$hits->{WIDTH} = '';
	$hits->{SCORE} = '';
	$hits->{INIT_CODON_POS} = '';
}

1;
