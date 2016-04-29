use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous',
);

my $stable_id = 'ENST00000367698';  #this is the stable_id of a human transcript

#get the adaptor to get the Transcript and slices from the database
my $transcript_adaptor =
  $registry->get_adaptor( 'homo_sapiens', 'core', 'transcript' );
my $slice_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Slice' );

#get the Transcript object
#my $transcript = $transcript_adaptor->fetch_by_stable_id($stable_id);
my @chromosomes = ('Y');
my $transcript  = get_transcripts_by_chromosomes( \@chromosomes );
print "Encontrados "
  . scalar @{$transcript}
  . " transcritos en los cromosomas "
  . join( ",", @chromosomes );
my @so_terms = ('start_lost');

#get the adaptor to get TranscriptVariation objects
my $trv_adaptor =
  $registry->get_adaptor( 'homo_sapiens', 'variation', 'transcriptvariation' );

#my $trvs = $trv_adaptor->fetch_all_by_VariationFeatures_SO_terms([$transcript], [\@so_terms]);
my $trvs = $trv_adaptor->fetch_all_by_Transcripts($transcript)
  ;    #get ALL effects of Variations in the Transcript
print "Encontrados " . scalar @{$trvs} . " variaciones\n";
$trvs = filter_by_so_terms( $trvs, \@so_terms );
print "Fltradas " . scalar @{$trvs} . " variaciones\n";

my @header = qw(GENE CODON_CHANGE AMINOACID_CHANGE  CONSEQUENCE SO_TERM SIFT POLYPHEN);
print join(",", @header) . "\n";
        
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

# Receive a transcript variation object and SO terms list.
# param 1 -> transcript variation object.
# param 2 -> SO terms list
# return a new transcript object variation filtering so_terms in param 1.
sub filter_by_so_terms {
    my $trvs         = $_[0];
    my $so_terms     = $_[1];
    my @filteredList = ();
    foreach my $tv ( @{$trvs} ) {
        my $tvas = $tv->get_all_alternate_TranscriptVariationAlleles();
        foreach my $tva ( @{$tvas} ) {
            my $ocs = $tva->get_all_OverlapConsequences();
            foreach my $oc ( @{$ocs} ) {
                if ( exists_in_list( $oc->SO_term, $so_terms ) ) {
                    push @filteredList, $tv;
                }
            }
        }
    }
    return \@filteredList;
}

# Checks if an element exists in the given list.
# param 0 -> Element to check.
# param 1 -> Reference to list.
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
# param 0 -> List with chromosome names.
# return list of transcripts.
sub get_transcripts_by_chromosomes {
    my $chromosomes = $_[0];
    my @transcripts = ();
    while ( my $chromosome = shift @{$chromosomes} ) {
        my $slice =
          $slice_adaptor->fetch_by_region( 'chromosome', $chromosome );
        my $aux = $transcript_adaptor->fetch_all_by_Slice($slice);
        push( @transcripts, @{$aux} );
    }
    return \@transcripts;
}
