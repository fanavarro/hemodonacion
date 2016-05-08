use strict;
use warnings;
use Bio::EnsEMBL::Registry;

#Get output file from input parameters
my $output = "";
if (scalar @ARGV == 1){
	$output = $ARGV[0];
} else {
	print "Usage: perl ensemblMining.pl outputFile";
	exit;
}
open(my $fh, '>', $output) or die "Could not open file '$output' $!";
print "Results will be printed in $output";

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
my $trv_adaptor = $registry->get_adaptor( 'homo_sapiens', 'variation', 'transcriptvariation' );

# Chromosomes to be treated
my @chromosomes = qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y);

# Sequence Ontology terms
my @so_terms = ('start_lost');

# CSV separator character configuration
my $csv_separator = ',';
my $in_field_separator = '-';
# Print csv header
my @header = qw(CHROMOSOME GENE_ID GENE_NAME TRANSCRIPT_ID VARIATION_NAME MINOR_ALLELE_FREQUENCY CODON_CHANGE AMINOACID_CHANGE  CONSEQUENCE SO_TERM SIFT POLYPHEN);
print $fh join($csv_separator, @header) . "\n";

# For each chromosome, get its variations with specified so terms.
foreach my $chromosome (@chromosomes) {
	get_variations_by_chromosome_so_terms($chromosome, \@so_terms);
}
close($fh);

# Get variations from chromosome (param 1) with the so term (param 2).
# param 1 -> Chromosome name
# param 2 -> SO term list
sub get_variations_by_chromosome_so_terms {
	my $chromosome = $_[0];
	my $so_terms = $_[1];
	
	# Get all transcripts from chromosome
	my $transcript  = get_transcripts_by_chromosome( $chromosome );
	# Get variation in transcripts with so terms
	my $trvs = $trv_adaptor->fetch_all_by_Transcripts_SO_terms($transcript, $so_terms);

	my $count = 1;
	foreach my $tv ( @{$trvs} ) {
	    print "Chromosome $chromosome\t Transcript Variation: $count/" . scalar @{$trvs} . "\n";
	    $count = $count + 1;
	    my $tvas = $tv->get_all_alternate_TranscriptVariationAlleles();
	    my $minor_allele_frequency = $tv->variation_feature->minor_allele_frequency ? $tv->variation_feature->minor_allele_frequency : '-';
	    foreach my $tva ( @{$tvas} ) {
		my @ensembl_consequences;
		my @so_consequences;
		my $ocs = $tva->get_all_OverlapConsequences();

		foreach my $oc ( @{$ocs} ) {
		    push @ensembl_consequences, $oc->display_term;
		    push @so_consequences,      $oc->SO_term;
		}

		my $sift     = $tva->sift_prediction;
		my $polyphen = $tva->polyphen_prediction;

		print $fh $chromosome . $csv_separator . $tv->transcript->get_Gene->stable_id . $csv_separator . 
		  $tv->transcript->get_Gene->external_name . $csv_separator . 
		  $tv->transcript->display_id . $csv_separator .
		  $tv->variation_feature->variation_name . $csv_separator .
		  $minor_allele_frequency . $csv_separator .
		  $tva->transcript_variation->get_reference_TranscriptVariationAllele->codon . "/" . $tva->codon . $csv_separator . $tva->pep_allele_string .
		  $csv_separator . join( $in_field_separator, @ensembl_consequences ) .
		  $csv_separator . join( $in_field_separator, @so_consequences );

		print $fh $csv_separator;
		if ( defined($sift) ) {
		    print $fh "$sift";
		}
		print $fh $csv_separator;
		if ( defined($polyphen) ) {
		    print $fh "$polyphen";
		}

		print $fh "\n";
	    }
	}
}

# Checks if an element exists in the given list.
# param 0 -> Element to check.
# param 1 -> List reference.
# returns true if elements exists inside the list and false if not.
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

# Get all transcripts that belongs to chromosomes in list.
# param 0 -> List reference with chromosome names.
# return list reference of transcripts.
sub get_transcripts_by_chromosomes {
    my $chromosomes = $_[0];
    my @transcripts = ();
    while ( my $chromosome = shift @{$chromosomes} ) {
        my $aux = get_transcripts_by_chromosome($chromosome);
        push( @transcripts, @{$aux} );
    }
    return \@transcripts;
}

# Get all transcripts that belongs to chromosome.
# param 0 -> Chromosome name.
# return list reference of transcripts.
sub get_transcripts_by_chromosome {
    my $chromosome = $_[0];
    my $slice = $slice_adaptor->fetch_by_region( 'chromosome', $chromosome );
    my $transcripts = $transcript_adaptor->fetch_all_by_Slice($slice);
    return $transcripts;
}
