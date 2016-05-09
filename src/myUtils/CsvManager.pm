package myUtils::CsvManager;
use strict;
use warnings;
use Fcntl qw(SEEK_SET);

# Receives a map with parameters
sub new{
	# get parameters
	my $class = shift @_;
	my %hash = @_;
	my $this = \%hash;
	open(my $file_descriptor, $this->{mode}, $this->{file_name}) or die "Could not open file '$this->{file_name}' $!";
	$this->{file_descriptor} = $file_descriptor;

	bless($this, $class);

	return $this;
}

sub close{
	my $this = shift;
	close($this->{file_descriptor});
}

sub readEntry{
	my $this = shift;
	
	my @fieldsName=@{$this->{fields}};
	my $fd = $this->{file_descriptor};
	my $currentPos = tell $fd;
	my $line ;
	if ($currentPos == 0){
		$line = <$fd>
	}
	$line = <$fd>;
	
	my %entry;
	if ($line){
		chomp $line;
		my @splitted = split($this->{csv_separator}, $line);
		for (my $i = 0; $i < scalar(@splitted); $i++){
			$entry{$fieldsName[$i]} = $splitted[$i];
		}
	}
	return %entry;
}

sub writeHeader{
	my $this = shift;
	my $header = $this->{fields};
	my $fd = $this->{file_descriptor} ;
	
	print $fd join($this->{csv_separator}, @{$header}) . "\n";
}
sub writeEntry{
	my $this = shift;
	my %entry = @_;
	my @fieldNames=@{$this->{fields}};
	my $fd = $this->{file_descriptor};
	my $csv_separator = $this->{csv_separator};
	
	my $string = "";
	for (my $i = 0; $i < scalar  keys  %entry; $i++){
		my $fieldName = $fieldNames[$i];
		$string = $string . $entry{$fieldName} . $csv_separator;
	}
	chop ($string);
	$string = $string . "\n";
	print $fd $string;
}
sub in_field_separator {
	my $this = shift;
	my $in_field_separator = shift;
	$this->{in_field_separator} = $in_field_separator if defined $in_field_separator;
	return $this->{in_field_separator};
}

sub csv_separator {
	my $this = shift;
	my $csv_separator = shift;
	$this->{csv_separator} = $csv_separator if defined $csv_separator;
	return $this->{csv_separator};
}

sub countEntries{
	my $this = shift;
	my $fd = $this->{file_descriptor};
	my $current_pos = tell $fd;
	seek($fd, 0, SEEK_SET);
	my $lines = 0;
	while (<$fd>){
		$lines = $lines + 1;
	}
	seek($fd, $current_pos, SEEK_SET);
	return $lines - 1; # Remove header.
}
1;
