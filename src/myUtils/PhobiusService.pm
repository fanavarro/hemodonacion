package myUtils::PhobiusService;
use strict;
use warnings;
use LWP::Simple::REST qw/http_post http_get/;
use XML::LibXML;

my $EMAIL = 'francisco.abad@um.es';

my $URL_BASE = 'http://www.ebi.ac.uk/Tools/services/rest/phobius/';
my $URL_RUN = $URL_BASE . 'run/';
my $URL_STATUS = $URL_BASE . 'status/';
my $URL_RESULT_TYPES = $URL_BASE . 'resulttypes/';
my $URL_RESULT = $URL_BASE . 'result/';

# Receives a hash with seq_id -> seq
# and returns info about signal peptide
# in these sequences.
sub get_info_signal_peptide{
	my %seqs = %{shift()};
	my $fasta = hash_to_fasta(\%seqs);
	my $job_id = http_post($URL_RUN, {email => $EMAIL, format => 'long', stype => 'protein', sequence => $fasta});
	print "Job id = $job_id\n";

	# Wait until job ends.
	my $job_status = http_get($URL_STATUS . $job_id);
	while ($job_status eq 'RUNNING'){
		sleep 2;
		$job_status = http_get($URL_STATUS . $job_id);
	}

	if ($job_status ne 'FINISHED'){
		# tratar error
		print "ERROR\nStatus = " . $job_status . "\n";
		return;
	}
	
	my $result_types_xml = http_get($URL_RESULT_TYPES . $job_id);
	my $result_types = XML::LibXML->load_xml(string => $result_types_xml);
	foreach my $result_type ($result_types->findnodes('//type')){
		my $result_type_id = $result_type->findvalue('./identifier');
		my $result = http_get($URL_RESULT . $job_id . '/' . $result_type_id);
		print "Result type = $result_type_id\n";
		print "Result = $result\n";
	}
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

1;
