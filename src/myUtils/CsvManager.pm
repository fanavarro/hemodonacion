package myUtils::CsvManager;
use strict;
use warnings;

# Receives a map with parameters
sub new{
	# get parameters
	my $class = shift @_;
	my %hash = @_;
	my $this = \%hash;
	open(my $file_descriptor, $this->{mode}, $this->{file_name}) or die "Could not open file '$this->{file_name}' $!";
	$this->{file_descriptor} = $file_descriptor;
	
	
	print "csv separator " . $this->{csv_separator} . "\n";
	print "in_field_separator " . $this->{in_field_separator} . "\n";
	print "file_name " . $this->{file_name} . "\n";
	print "fields " . $this->{fields} . "\n";
	print "mode " . $this->{mode} . "\n";

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

	chomp $line;
	my @splitted = split($this->{csv_separator}, $line);
	for (my $i = 0; $i < scalar(@splitted); $i++){
		$entry{$fieldsName[$i]} = $splitted[$i];
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

1;
