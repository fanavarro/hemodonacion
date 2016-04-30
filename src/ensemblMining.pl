use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous',
);


#get the adaptor to get the Transcript, slices and transcript variation from the database
my $transcript_adaptor = $registry->get_adaptor( 'homo_sapiens', 'core', 'transcript' );
my $slice_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Slice' );
my $trv_adaptor = $registry->get_adaptor( 'homo_sapiens', 'variation', 'transcriptvariation' );

my @chromosomes = ('Y');
my @so_terms = ('start_lost');



my @header = qw(GENE CODON_CHANGE AMINOACID_CHANGE  CONSEQUENCE SO_TERM SIFT POLYPHEN);
print join(",", @header) . "\n";
foreach my $chromosome (@chromosomes) {
	get_variations_by_chromosome_so_terms($chromosome, \@so_terms);
}

# Get variations from chromosome (param 1) with the so term (param 2).
# param 1 -> Chromosome name
# param 2 -> SO term list
sub get_variations_by_chromosome_so_terms {
	my $chromosome = $_[0];
	my $so_terms = $_[1];
	
	# Get all transcripts from chromosome
	my $transcript  = get_transcripts_by_chromosome( $chromosome );
	# Get variation in transcripts with so terms
	my $trvs = $trv_adaptor->fetch_all_by_Transcripts_SO_terms($transcript, \@so_terms);

	
	foreach my $tv ( @{$trvs} ) {
	    my $tvas = $tv->get_all_alternate_TranscriptVariationAlleles();

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
		
		print $tv->transcript->get_Gene->stable_id, ",",  $tva->transcript_variation->get_reference_TranscriptVariationAllele->codon,
		  "/",$tva->codon,
		  ",", $tva->pep_allele_string,
		  ",", join( "-", @ensembl_consequences ),
		  ",", join( "-", @so_consequences );

		if ( defined($sift) ) {
		    print ",$sift";
		}
		if ( defined($polyphen) ) {
		    print ",$polyphen";
		}

		print "\n";
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
