use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $transcript_id = "";
if (scalar @ARGV == 1){
	$transcript_id = $ARGV[0];
} else {
	print "Usage: perl test.pl transcript_id\n";
	exit;
}

# Registry configuration
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous',
);
# Set the flag to make sure that the database connection is dropped if not being used on each database.
# $registry->set_disconnect_when_inactive();
# Set the flag to make sure that the database connection is not lost before it's used. This is useful for long running jobs (over 8hrs).
$registry->set_reconnect_when_lost();

# Get the adaptor to get the Transcript, slices and transcript variation from the database
my $transcript_adaptor = $registry->get_adaptor( 'homo_sapiens', 'core', 'transcript' );

my $transcript = $transcript_adaptor->fetch_by_stable_id($transcript_id);

if (defined $transcript){
	print "status -> " . $transcript->status . "\n";
	if(defined ($transcript->havana_transcript)){
		print "havana transcript -> " . $transcript->havana_transcript->display_id . "\n";
	}
	print "source -> " . $transcript->source . "\n";
	print "tsl -> " . $transcript->tsl . "\n";
	print "biotype -> " . $transcript->biotype . "\n";
}
