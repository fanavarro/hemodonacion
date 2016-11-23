use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use IO::File;
use myUtils::SeqUtils;

#Get output file from input parameters
my $output = "";
if (scalar @ARGV == 1){
	$output = $ARGV[0];
} else {
	print "Usage: perl get_kozak_seqs.pl outputFile\n";
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
my $slice_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Slice' );


# Transcript constraints. Source is ensembl_havana, which has experimentally observed transcripts.
# Biotype is protein coding.
# Status is known.
my $transcript_constraints = "source = 'ensembl_havana' and biotype = 'protein_coding' and t.status = 'KNOWN'";


my @chromosomes = qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y);

open my $fd, '>', $output or die "Could not open file $output $!";
my $positions_after = 10;
my $positions_before = 10;
foreach my $chromosome (@chromosomes){
    my $transcript_count = 1;
    my $transcripts = get_transcripts_by_chromosome_with_constraints($chromosome, $transcript_constraints);
    my $total_transcripts = scalar @{$transcripts};
    foreach my $transcript (@{$transcripts}){
        my $id = $transcript->stable_id;
        print "Chromosome $chromosome\tTranscript $transcript_count of $total_transcripts ($id)\n";
        my $cdna = $transcript->seq->seq;
        my $cds = $transcript->translateable_seq;
        my $kozak_context = myUtils::SeqUtils::get_kozak_context($cdna, $cds, $positions_before, $positions_after);
        if (defined $kozak_context){
            print $fd $kozak_context . "\n";
        }
        if (!$transcript->is_canonical){
            print "se ha colado un transcrito no canonico.";
        }
        $transcript_count = $transcript_count + 1;
    }
}
close($fd);


# Get all transcripts that belongs to chromosome.
# param 0 -> Chromosome name.
# param 1 -> Constraints to pass to DB.
# return list reference of transcripts.
sub get_transcripts_by_chromosome_with_constraints {
    my $chromosome = $_[0];
    my $constraints = $_[1];
    my $slice = $slice_adaptor->fetch_by_region( 'chromosome', $chromosome );
    my $transcripts = $transcript_adaptor->fetch_all_by_Slice_constraint($slice, $constraints);
    my $canonical_transcripts = filter_canonical_transcripts($transcripts);
    return $canonical_transcripts;
}


sub filter_canonical_transcripts{
    my $transcripts = shift;
    my @canonical_transcripts;
    foreach my $transcript (@{$transcripts}){
        if ($transcript->is_canonical){
            push @canonical_transcripts, $transcript;
        }
    }
    return \@canonical_transcripts;
}






