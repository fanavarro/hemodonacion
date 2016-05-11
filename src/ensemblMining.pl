use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use myUtils::CsvManager;

#Get output file from input parameters
my $output = "";
if (scalar @ARGV == 1){
	$output = $ARGV[0];
} else {
	print "Usage: perl ensemblMining.pl outputFile";
	exit;
}
open(my $fh, '>', $output) or die "Could not open file '$output' $!";
print "Results will be printed in $output\n";

# CSV file configuration
my @fields = qw(CHROMOSOME GENE_ID GENE_NAME TRANSCRIPT_ID TRANSCRIPT_BIOTYPE CDS_ERRORS PROTEIN_ID VARIATION_NAME MINOR_ALLELE_FREQUENCY CODON_CHANGE AMINOACID_CHANGE NEXT_MET CONSEQUENCE SO_TERM SIFT POLYPHEN);
my $out_csv = myUtils::CsvManager->new (
	fields    => \@fields,
	csv_separator   => ',',
	in_field_separator    => '-',
	file_name => $output,
	mode => '>'
);
$out_csv->myUtils::CsvManager::writeHeader();

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
# my @chromosomes = qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y);
my @chromosomes = qw(Y);

# Sequence Ontology terms
# start_lost -> a codon variant that changes
# at least one base of the canonical start codon.
my @so_terms = ('start_lost');

# For each chromosome, get its variations with specified so terms.
foreach my $chromosome (@chromosomes) {
	get_variations_by_chromosome_so_terms($chromosome, \@so_terms, $out_csv);
}

$out_csv->myUtils::CsvManager::close();

# Get variations from chromosome (param 1) with the so term (param 2).
# param 1 -> Chromosome name
# param 2 -> SO term list
# param 3 -> CsvManager object to print results.
sub get_variations_by_chromosome_so_terms {
	my $chromosome = $_[0];
	my $so_terms = $_[1];
	my $out_csv = $_[2];
	# Get all transcripts from chromosome
	my $transcript  = get_transcripts_by_chromosome( $chromosome );
	# Get variation in transcripts with so terms
	my $trvs = $trv_adaptor->fetch_all_by_Transcripts_SO_terms($transcript, $so_terms);

	my $count = 1;
	foreach my $tv ( @{$trvs} ) {
	    print "Chromosome $chromosome\t Transcript Variation: $count/" . scalar @{$trvs} . "\n";
	    $count = $count + 1;
	    
	    my $minor_allele_frequency = $tv->variation_feature->minor_allele_frequency ? $tv->variation_feature->minor_allele_frequency : '-';
	    my $cds_errors = get_cds_errors($tv->transcript);

	    my $tvas = $tv->get_all_alternate_TranscriptVariationAlleles();
	    foreach my $tva ( @{$tvas} ) {
		my %entry;
	    	
		my @ensembl_consequences;
		my @so_consequences;
		my $ocs = $tva->get_all_OverlapConsequences();

		foreach my $oc ( @{$ocs} ) {
		    push @ensembl_consequences, $oc->display_term;
		    push @so_consequences,      $oc->SO_term;
		}

		my $sift     = $tva->sift_prediction;
		my $polyphen = $tva->polyphen_prediction;
		$entry{'CHROMOSOME'} = $chromosome;
		$entry{'GENE_ID'} = $tv->transcript->get_Gene->stable_id;
		$entry{'GENE_NAME'} = $tv->transcript->get_Gene->external_name;
		$entry{'TRANSCRIPT_ID'} = $tv->transcript->display_id;
		$entry{'TRANSCRIPT_BIOTYPE'} = $tv->transcript->biotype;
		$entry{'CDS_ERRORS'} = $cds_errors;
		$entry{'PROTEIN_ID'} = $tv->transcript->translation->display_id;
		$entry{'VARIATION_NAME'} = $tv->variation_feature->variation_name;
		$entry{'MINOR_ALLELE_FREQUENCY'} = $minor_allele_frequency;
		$entry{'CODON_CHANGE'} = $tva->display_codon_allele_string;
		$entry{'AMINOACID_CHANGE'} = $tva->pep_allele_string;
		#get_next_met($tva);
		get_variation_seq($tva);
		$entry{'NEXT_MET'} = '-';
		$entry{'CONSEQUENCE'} = join( $out_csv->myUtils::CsvManager::in_field_separator(), @ensembl_consequences );
		$entry{'SO_TERM'} = join( $out_csv->myUtils::CsvManager::in_field_separator(), @so_consequences );
		$entry{'SIFT'} = '';
		$entry{'POLYPHEN'} = '';
		

		if ( defined($sift) ) {
		    $entry{'SIFT'} = "$sift";
		}
		if ( defined($polyphen) ) {
		    $entry{'POLYPHEN'} = "$polyphen";
		}
		$out_csv->myUtils::CsvManager::writeEntry(%entry);

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

# param 0 -> TranscriptVariationAllele object.
# return -> Position of first Met found relative to peptide.
sub get_next_met{
    my $transcript_variation_allele = $_[0];
    my $seq = length $transcript_variation_allele->feature_seq;
    print $seq . "\n";
}

sub get_variation_seq{
    my $tva = $_[0];
    # translateable_seq returns the coding part of the transcript
    # (it removes introns and 5' and 3' utr)
    # my $seq = $tva->transcript->translateable_seq;
    # seq contains 5' and 3' regions.
    my $seq = $tva->transcript->seq->seq;
    my $variation_start = $tva->transcript_variation->cdna_start - 1;
    my $variation_end = $tva->transcript_variation->cdna_end - 1;
    # If is a deletion, feature_seq is '-', so we will use '' instead
    # to build the final sequence.
    my $feature_seq = $tva->feature_seq eq "-" ? "" : $tva->feature_seq;
    substr($seq, $variation_start, $variation_end - $variation_start + 1) = $feature_seq;
    
    return $seq;
}

# param 0 -> Transcript object
# param 1 -> Attribute code.
# return -> Attribute name.
sub get_attribute {
    my ($transcript, $code) = @_;
    my ($attr) = @{$transcript->get_all_Attributes($code)};
    my $attribute_name;
    if($attr) {
      $attribute_name = $attr->name();
    }
    return $attribute_name;
}

sub get_cds_errors{
    my $transcript = $_[0];
    my $cds_errors="";
    my $five_cds = get_attribute($transcript, 'cds_start_NF');
    my $three_cds = get_attribute($transcript, 'cds_end_NF');
    if ($five_cds){
	$cds_errors = $cds_errors . $five_cds . '-';
    }
    if ($three_cds){
	$cds_errors = $cds_errors . $three_cds . '-';
    }
    chop $cds_errors;
    return $cds_errors;
}


