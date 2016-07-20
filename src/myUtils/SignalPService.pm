package myUtils::SignalPService;
use strict;
use warnings;
use SOAP::Lite;
use base 'Class::Singleton';

sub _new_instance {
	my $class = shift;
	my $this  = bless { }, $class;
	
	my $signalp = SOAP::Lite
	-> service('http://www.cbs.dtu.dk/ws/SignalP/SignalP_3_1.wsdl')
	-> on_fault(
			sub { my($self, $res) = @_;
				die
				 "faultcode:", ref $res ? $res->faultcode : $self->transport->status, "\n" ,
				 "faultstring:", ref $res ? $res->faultstring : $self->transport->status, "\n";
			}
		);

	$this->{entry_point} = $signalp;
	return $this;
}

# Receives a hash of id_sequence->AA_sequence
# and returns signalP info.
sub get_signal_peptide_info{
	my $this = shift;
	my $entry_point = $this->{entry_point};
	my %seqs = %{(shift)};
	my $input = _build_signalp_input(\%seqs);
	# WS returns a job, which we have to check in other call until it is finished.
	my $job = $entry_point->runService($input);
	my $jobid = SOAP::Data->name(
		'job' => \SOAP::Data->value(
			SOAP::Data->name('jobid' => $job->{jobid}
		)
	));
	$job = $this->_wait_until_job_done($job);

	my $results = $entry_point->fetchResult($jobid);

	my @predictions = ();
	if (ref ( $results->{prediction}->{gff_record} ) eq 'ARRAY') {
		@predictions = @{$results->{prediction}->{gff_record}};
	} else {
		@predictions = ($results->{prediction}->{gff_record});
	}

	# test print results
	foreach my $prediction (@predictions){
		my %pred = %{$prediction};
		while (my($key, $value) = each %pred){
			print "$key -> $value\n";
		}
		print "\n";
	}

}

# Receives a hash of id_sequence->AA_sequence
# and returns an array of seqs to be envelopped inside
# WS input request.
# param0 -> hash of id_sequence->AA_sequence.
# return array of entries of sequences to use as input for signalp ws.
sub _build_signalp_seqs_input{
	my %seqs = %{shift()};
	my @signalp_seqs_input;
	while (my($seq_id, $seq) = each %seqs){
		my $entry = ( SOAP::Data->name('entry' => \SOAP::Data->value(SOAP::Data->name('ident' => $seq_id) ,SOAP::Data->name('seq' => $seq)))) ;
		push @signalp_seqs_input, $entry;
	}
	return @signalp_seqs_input;
	
}

# Receives a hash of id_sequence->AA_sequence
# and returns an input object to use as signalP WS.
# param0 -> hash of id_sequence->AA_sequence.
# return array of entries of sequences to use as input for signalp ws.
sub _build_signalp_input{
	my %seqs = %{shift()};
	my @seqs_signalp_input = _build_signalp_seqs_input(\%seqs);
	my $input = SOAP::Data->name(
			'parameters' => \SOAP::Data->value(
				SOAP::Data->name('organism' => "euk"),
				SOAP::Data->name('sequences' => \SOAP::Data->value(@seqs_signalp_input))
				)
			);
	return $input;
}

sub _wait_until_job_done{
	my $this = shift;
	my $job = shift;
	my $entry_point = $this->{entry_point};
	my $jobid = SOAP::Data->name(
		'job' => \SOAP::Data->value(
			SOAP::Data->name('jobid' => $job->{jobid}
		)
	));

	$job = $entry_point -> pollQueue($jobid);

	while ($job->{status} !~ /FINISHED/) {
		print  "The job is: $job->{status}\n";
		$job = $entry_point->pollQueue($jobid);
		sleep 5;
	}
	return $job;
}

1;
