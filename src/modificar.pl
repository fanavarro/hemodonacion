use strict;
use warnings;

open my $fd, '<', "kozaks.txt" or die "Could not open file $!";

my $id=1;

while (my $line = <$fd>){
	chomp $line;
	print ">SEQUENCE-$id\n";
	print $line . "\n";
	$id = $id + 1;
}

close($fd);
