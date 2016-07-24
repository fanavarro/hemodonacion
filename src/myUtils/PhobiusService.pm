package myUtils::PhobiusService;
use strict;
use warnings;
use LWP::Simple::REST qw/http_post http_get/;
use XML::LibXML;

# REST api constants
my $EMAIL = 'francisco.abad@um.es';
my $URL_BASE = 'http://www.ebi.ac.uk/Tools/services/rest/phobius/';
my $URL_RUN = $URL_BASE . 'run/';
my $URL_STATUS = $URL_BASE . 'status/';
my $URL_RESULT_TYPES = $URL_BASE . 'resulttypes/';
my $URL_RESULT = $URL_BASE . 'result/';

# Standalone application constants
my $PHOBIUS_LOCAL_DIR = 'myUtils/phobius/phobius.pl';

my $RESULT_SEP = '//';

# Receives a hash with seq_id -> seq
# and returns info about signal peptide
# in these sequences using phobius rest service.
sub get_info_signal_peptide_rest{
	my %seqs = %{shift()};
	my $fasta = hash_to_fasta(\%seqs);
	my $job_id = http_post($URL_RUN, {email => $EMAIL, format => 'long', stype => 'protein', sequence => $fasta});
	# print "Job id = $job_id\n";

	# Wait until job ends.
	my $job_status = http_get($URL_STATUS . $job_id);
	while ($job_status eq 'RUNNING'){
		sleep 2;
		$job_status = http_get($URL_STATUS . $job_id);
	}

	if ($job_status ne 'FINISHED'){
		# tratar error
		print "ERROR\nStatus = " . $job_status . "\n";
		return undef;
	}
	
	# my $result_types_xml = http_get($URL_RESULT_TYPES . $job_id);
	# my $result_types = XML::LibXML->load_xml(string => $result_types_xml);
	# foreach my $result_type ($result_types->findnodes('//type')){
	# 	my $result_type_id = $result_type->findvalue('./identifier');
	#	my $result = http_get($URL_RESULT . $job_id . '/' . $result_type_id);
	#	print "Result type = $result_type_id\n";
	#	print "Result = $result\n";
	# }
	my $result_type_id = 'out';
	my $result_text = http_get($URL_RESULT . $job_id . '/' . $result_type_id);
	chomp($result_text );

	return text_result_to_hash($result_text);
}

# Receives a hash with seq_id -> seq
# and returns info about signal peptide
# in these sequences using local phobius.
sub get_info_signal_peptide_local{
	my %seqs = %{shift()};
	my $fasta = hash_to_fasta(\%seqs);
	
	# use qx to call local phobius program.
	my $phobius_output = qx(echo "$fasta" |  $PHOBIUS_LOCAL_DIR -long);
	chomp($phobius_output );

	return text_result_to_hash($phobius_output);
}

# Convert a hash of type id_seq -> seq
# to a fasta string.
# param0 -> Hash whose key is a seq id and whose value is the seq.
# return a fasta formatted string .
sub hash_to_fasta{
	my %seqs = %{shift()};
	
	my $fasta = '';

	while (my($seq_id, $seq) = each %seqs){
		$fasta = $fasta . '>' . $seq_id . "\n";
		$fasta = $fasta . $seq . "\n";
	}
	return $fasta;
}

sub text_result_to_hash{
	my $result_text = shift;
	my %result_hash = ();
	if (defined($result_text)){
		# Iterate over each sequence results
		foreach my $seq_text (split($RESULT_SEP, $result_text)){
			# Remove newlines from the beginning and the ending
			chomp($seq_text);
			$seq_text =~ s/^\n//; 
			my @lines = split("\n", $seq_text);
			my $id = (split " ", $lines[0])[1];
			my @feature_list = ();
			# Iterate over each feature of the current sequence
			for(my $i = 1; $i < scalar(@lines); $i++){
				my %entry;
				# Separate line in words separated by two or more whitespaces.
				my @splitted_line = split (/\s{2,}/, $lines[$i]);
				$entry{'TYPE'} = $splitted_line[1];
				$entry{'START'} = $splitted_line[2];
				$entry{'END'} = $splitted_line[3];
				$entry{'LOCATION'} = $splitted_line[4];
				push (@feature_list, \%entry);
			}
			$result_hash{$id} = \@feature_list;
		}
	}

	if (%result_hash){
		return \%result_hash;
	}
	else {
		return undef;
	}
	
}
1;
