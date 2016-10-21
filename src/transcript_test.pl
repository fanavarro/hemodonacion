use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use v5.10;

my @transcript_ids;

if (scalar @ARGV > 0){
	@transcript_ids = @ARGV;
} else {
	say "Usage: perl transcript_test id1 id2 id3 ...";
	exit;
}

# Registry configuration
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous',
);

my $transcript_adaptor = $registry->get_adaptor( 'homo_sapiens', 'core', 'transcript' );
# my $trv_adaptor = $registry->get_adaptor( 'homo_sapiens', 'variation', 'transcriptvariation' );

my $transcripts = $transcript_adaptor->fetch_all_by_stable_id_list(\@transcript_ids);

foreach my $transcript (@{$transcripts}){
	my $havana = $transcript->havana_transcript;
	my $havana_text = "";
	if (defined($havana)){
		$havana_text = $havana->display_id;
	}
	say $transcript->stable_id . "  " . $transcript->status . "  " . $havana_text . "  " . $transcript->gencode_basic;
	
}


sub printHash{
	my %hash = %{$_[0]};
	while (my ($key, $value) = each(%hash)) {
		print "$key -> $value\n";
	}
}

sub exists_in_list {
    my $check = $_[0];
    my $list  = $_[1];
    foreach my $element ( @{$list} ) {
        if ( $check eq $element ) {
            return 1;
        }
    }
    return 0;
}
