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

	# If fields are not in params, we extract it from file.
	if (!$this->{fields}){
		$this->{fields} = _getFieldsFromFile($file_descriptor, $this->{csv_separator});
	}

	bless($this, $class);

	return $this;
}

# Extracts fields from file header.
sub _getFieldsFromFile{
	my $fd = $_[0];
	my $separator = $_[1];
	my $currentPos = tell $fd;
	seek($fd, 0, SEEK_SET);
	my $line = <$fd>;
	chomp($line);
	my @fields = split($separator, $line);
	seek($fd, $currentPos, SEEK_SET);
	return \@fields;
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
	for (my $i = 0; $i < scalar  @fieldNames; $i++){
		my $fieldName = $fieldNames[$i];
		my $value = defined($entry{$fieldName}) ? $entry{$fieldName} : '';
		$string = $string . $value . $csv_separator;
	}
	chop ($string);
	$string = $string . "\n";
	print $fd $string;
}
sub writeEntries{
	my $this = shift;
	my $entry_list = shift;

	foreach my $entry (@{$entry_list}){
		$this->writeEntry(%{$entry});
	}
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

sub fields {
	my $this = shift;
	my @fields = @_;
	$this->{fields} = \@fields if scalar(@fields) > 0;
	return $this->{fields};
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

# Reads a csv file and saves it into
# a list of hashes, where each item
# on the list is a row and each key
# of the hash is a column.
# param 0 -> this
# return list of hashes
sub readCsv{
	my $this = shift;
	my $fd = $this->{file_descriptor};
	my $rows = [];
	my $currentPos = tell $fd;
	seek($fd, 0, SEEK_SET);
	while(my %entry = $this->readEntry()){
		push @{$rows}, \%entry;
	}
	seek($fd, $currentPos, SEEK_SET);
	return $rows;
}
# Compare this csv with other csv. Returns hash reference
# with two keys: 'this' contains exclusive csv entries from
# this csv. 'other' contains exclusive csv entries from other
# csv.
sub diffCsv{
	my $this = shift;
	my $other = shift;

	my $refhash = {};
	$refhash->{'this'} = [];
	$refhash->{'other'} = [];
	
	my $entries1 = $this->readCsv();
	my $entries2 = $other->readCsv();

	my $count = 1;
	my $total = scalar @{$entries1} + scalar @{$entries2};

	foreach my $entry (@{$entries1}){
		print $count . "/" . $total . "\n";
		$count = $count + 1;
		if(! _existsEntry($entries2, $entry)){
			push @{$refhash->{'this'}}, $entry;
		}
	}

	foreach my $entry (@{$entries2}){
		print $count . "/" . $total . "\n";
		$count = $count + 1;
		if(! _existsEntry($entries1, $entry)){
			push @{$refhash->{'other'}}, $entry;
		}
	}


	return $refhash;
}

# Return true if entry exists
# into the csv.
sub _existsEntry{
	my $entries = shift;
	my $entry = shift;
	foreach my $csvEntry (@{$entries}){
		if (equalsEntry($csvEntry, $entry)){
			return 1;
		}
	}
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





