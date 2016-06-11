use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $transcript_id = '';
if (scalar @ARGV == 1){
	$transcript_id = $ARGV[0];
} else {
	print "Usage: perl vartest2.pl transcriptId";
	exit;
}

# Registry configuration
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous',
);

my $transcript_adaptor = $registry->get_adaptor( 'homo_sapiens', 'core', 'transcript' );
# my $slice_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Slice' );
# my $trv_adaptor = $registry->get_adaptor( 'homo_sapiens', 'variation', 'transcriptvariation' );
# my $var_adaptor = $registry->get_adaptor( 'homo_sapiens', 'variation', 'variation' );
# my $vfa = $registry->get_adaptor( 'homo_sapiens', 'variation', 'variationfeature' );

# $var = $var_adaptor->fetch_by_name('CM034685');

my $transcript = $transcript_adaptor->fetch_by_stable_id($transcript_id);

print "Display Id -> " . $transcript->display_id() . "\n";
print "Description -> " . $transcript->description() . "\n";
print "External name -> " . $transcript->external_name() . "\n";


foreach my $dbentry ( @{$transcript->get_all_DBEntries('RefSeq_mRNA')} ){
	print "\tDatabase ->" . $dbentry->database . "\n";
	print "\tDB display name ->" . $dbentry->db_display_name . "\n";
	print "\tDB ID ->" . $dbentry->dbID . "\n";
	print "\tDB name ->" . $dbentry->dbname . "\n";
	print "\tDescription ->" . $dbentry->description . "\n";
	print "\tDisplay id ->" . $dbentry->display_id . "\n";
	print "\tPrimary id ->" . $dbentry->primary_id . "\n";
	print "\tOptional id ->" . $dbentry->optional_id . "\n";
	print "\n";
}
