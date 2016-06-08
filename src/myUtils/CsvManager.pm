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

# Compare this csv with other csv. Returns hash reference
# with two keys: 'this' contains exclusive csv entries from
# this csv. 'other' contains exclusive csv entries from other
# csv.
sub diffCsv{
	my $this = shift;
	my $other = shift;
	my $fd1 = $this->{file_descriptor};
	my $current_pos1 = tell $fd1;
	my $fd2 = $this->{file_descriptor};
	my $current_pos2 = tell $fd2;

	my $refhash = {};
	$refhash->{'this'} = [];
	$refhash->{'other'} = [];

	seek($fd1, 0, SEEK_SET);
	while(my %entry = $this->readEntry()){
		if(! $other->existsEntry(\%entry)){
			push @{$refhash->{'this'}}, \%entry;
		}
	}

	seek($fd2, 0, SEEK_SET);
	while(my %entry = $other->readEntry()){
		if(! $this->existsEntry(\%entry)){
			push @{$refhash->{'other'}}, \%entry;
		}
	}

	seek($fd1, $current_pos1, SEEK_SET);
	seek($fd2, $current_pos2, SEEK_SET);
	return $refhash;
}

# Return true if entry exists
# into the csv.
sub existsEntry{
	my $this = shift;
	my $entry = shift;
	my $fd = $this->{file_descriptor};
	my $current_pos = tell $fd;
	seek($fd, 0, SEEK_SET);
	while (my %csvEntry = $this->readEntry()){
		if (equalsEntry(\%csvEntry, $entry)){
			seek($fd, $current_pos, SEEK_SET);
			return 1;
		}
	}
	seek($fd, $current_pos, SEEK_SET);
	return 0;
}

sub equalsEntry{
	my $entry1 = shift;
	my $entry2 = shift;
	if (keys %$entry1 ne keys %$entry2){
		return 0;
	}

	foreach my $key (keys %$entry1){
		if ($entry1->{$key} ne $entry2->{$key}){
			return 0;
		}
	}
	return 1;
}
1;





