use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use myUtils::CsvManager;
use myUtils::Publication;
use File::Path qw(make_path);

my $CODON_LENGTH = 3;
my $MET = 'ATG';
my @STOP_CODONS = qw(TAG TAA TGA);

# Needed to write to nohup.out in real time
STDOUT->autoflush(1);

#Get output file from input parameters
my $output = "";
if (scalar @ARGV == 1){
	$output = $ARGV[0];
} else {
	print "Usage: perl ensemblMining.pl outputFile";
	exit;
}
#open(my $fh, '>', $output) or die "Could not open file '$output' $!";
print "Results will be printed in $output\n";

# CSV file configuration
my @fields = qw(CHROMOSOME GENE_ID GENE_NAME TRANSCRIPT_ID TRANSCRIPT_BIOTYPE CDS_ERRORS PROTEIN_ID VARIATION_NAME SOURCE TRANSCRIPT_VARIATION_ALLELE_DBID MINOR_ALLELE_FREQUENCY CODON_CHANGE AMINOACID_CHANGE FIRST_MET_POSITION STOP_CODON_POSITION MUTATED_SEQUENCE_LENGTH READING_FRAME_STATUS CONSEQUENCE PHENOTYPE SO_TERM SIFT POLYPHEN PUBLICATIONS);
my $out_csv = myUtils::CsvManager->new (
	fields    => \@fields,
	csv_separator   => "\t",
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
#my @chromosomes = qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y);
my @chromosomes = qw(Y);

# Sequence Ontology terms
# start_lost -> a codon variant that changes
# at least one base of the canonical start codon.
my @so_terms = ('start_lost');
my $variation_pos_at_peptide = 1;
# For each chromosome, get its variations with specified so terms.
foreach my $chromosome (@chromosomes) {
	#get_variations_by_chromosome_so_terms($chromosome, \@so_terms, $out_csv);
	get_variations_by_chromosome_peptide_position($chromosome, $variation_pos_at_peptide, $out_csv);
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

	fill_csv($chromosome, $trvs, $out_csv);
}
# Get variations from chromosome (param 1) at the peptide position (param 2).
# param 1 -> Chromosome name
# param 2 -> Peptide position
# param 3 -> CsvManager object to print results.
sub get_variations_by_chromosome_peptide_position {
	my $chromosome = $_[0];
	my $peptide_position = $_[1];
	my $out_csv = $_[2];
	# Get all transcripts from chromosome
	my $transcript  = get_transcripts_by_chromosome( $chromosome );

	# Get variation in transcripts adding the position constraint
        my $constraint = "(translation_start=$peptide_position or translation_end=$peptide_position)";
	my $trvs = $trv_adaptor->fetch_all_by_Transcripts_with_constraint($transcript, $constraint);

	fill_csv($chromosome, $trvs, $out_csv);
}

# Write info about TranscriptVariation list into csv file.
# param 0 -> Chromosome
# param 1 -> TranscriptVariation ref list object.
# param 2 -> CsvManager object.
sub fill_csv{
    my $chromosome = $_[0];
    my $trvs = $_[1];
    my $out_csv = $_[2];

    my $count = 1;
    foreach my $tv ( @{$trvs} ) {
	print "Chromosome $chromosome\t Transcript Variation: $count/" . scalar @{$trvs} . "\n";
	$count = $count + 1;
	
	my $variation = $tv->variation_feature->variation;
	my $minor_allele_frequency = $tv->variation_feature->minor_allele_frequency ? $tv->variation_feature->minor_allele_frequency : '-';
	my $cds_errors = get_cds_errors($tv->transcript);
	my $phenotype_info = get_phenotype_info($variation);
	my $publications_info = get_publications_info($variation);
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
            my $seq_info = get_sequence_info($tva);
            $entry{'CHROMOSOME'} = $chromosome;
            $entry{'GENE_ID'} = $tv->transcript->get_Gene->stable_id;
            $entry{'GENE_NAME'} = $tv->transcript->get_Gene->external_name;
            $entry{'TRANSCRIPT_ID'} = $tv->transcript->display_id;
            $entry{'TRANSCRIPT_BIOTYPE'} = $tv->transcript->biotype;
            $entry{'CDS_ERRORS'} = $cds_errors;
            $entry{'PROTEIN_ID'} = $tv->transcript->translation->display_id;
            $entry{'VARIATION_NAME'} = $tv->variation_feature->variation_name;
            $entry{'SOURCE'} = $variation->source_name;
            $entry{'TRANSCRIPT_VARIATION_ALLELE_DBID'} = $tva->dbID;
            $entry{'MINOR_ALLELE_FREQUENCY'} = $minor_allele_frequency;
            $entry{'CODON_CHANGE'} = defined($tva->display_codon_allele_string) ? $tva->display_codon_allele_string : ' ';
            $entry{'AMINOACID_CHANGE'} = defined($tva->pep_allele_string) ? $tva->pep_allele_string : ' ';
            $entry{'FIRST_MET_POSITION'} = $seq_info->{'first_met_position'};
            $entry{'STOP_CODON_POSITION'} = $seq_info->{'stop_codon_position'};
            $entry{'MUTATED_SEQUENCE_LENGTH'} = $seq_info->{'seq_length'};
            $entry{'READING_FRAME_STATUS'} = $seq_info->{'reading_frame'};
            $entry{'CONSEQUENCE'} = join( '-', @ensembl_consequences );
            $entry{'PHENOTYPE'} = $phenotype_info;
            $entry{'SO_TERM'} = join( '-', @so_consequences );
            $entry{'SIFT'} = '';
            $entry{'POLYPHEN'} = '';
            $entry{'PUBLICATIONS'} = $publications_info;
	

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
# test borrar
#    $transcripts = $transcript_adaptor->fetch_all_by_stable_id_list(["ENST00000414219"]);
    return $transcripts;
}

# Extract information about variation sequence.
# param 0 -> TranscriptVariationAllele object.
# return -> Hash with the following keys:
# 'first_met_position': position of first Met found relative to peptide.
# 'reading_frame': Conserved or Lost.
# 'stop_codon_position': Position of the first stop codon (keeping reading frame).
# 'seq_length': Length of the mutated sequence.
sub get_sequence_info{
    my $tva = $_[0];
    my $hash_seq_info = {};
    $hash_seq_info->{'first_met_position'} = ' ';
    $hash_seq_info->{'reading_frame'} = ' ';
    $hash_seq_info->{'stop_codon_position'} = ' ';
    $hash_seq_info->{'seq_length'} = ' ';
    my $source = $tva->variation_feature->variation->source;
    if (defined($source) && $source->name eq 'dbSNP'){
        my $seq = get_variation_cds_seq($tva);
        $hash_seq_info->{'seq_length'} = length($seq);
        my $first_met_pos = index($seq, $MET);
        if ($first_met_pos != -1){
            my $stop_codon_pos = get_stop_codon_position($seq);
            my $cut_seq = substr($seq, $first_met_pos);
	    my $reading_frame = $first_met_pos % $CODON_LENGTH == 0 ? 'Conserved' : 'Lost';
	    $hash_seq_info->{'first_met_position'} = $first_met_pos;
            $hash_seq_info->{'stop_codon_position'} = $stop_codon_pos != -1 ? $stop_codon_pos : 'No Stop';
	    $hash_seq_info->{'reading_frame'} = $reading_frame;

            # if an ORF exists into the seq, seq file is generated.
            if($stop_codon_pos != -1){
                my $orf = substr($seq, $first_met_pos, $stop_codon_pos - $first_met_pos + 3);
                generate_variation_seq_files($tva, $orf);
            }
        }
    }
    return $hash_seq_info;
}

# Get the position of the first stop codon following
# the reading frame.
# param 1: string with the sequence
# return the position of the first stop codon found after
# first MET keeping the reading frame.
sub get_stop_codon_position{
    my $seq = $_[0];
    # get the first met in the sequence
    my $init_pos = index($seq, $MET);
    if ($init_pos != -1){
	# Search stop codons keeping the reading frame
        for (my $i = $init_pos; $i < length($seq); $i = $i + $CODON_LENGTH){
	    my $codon = substr($seq, $i, $CODON_LENGTH);
	    if (exists_in_list($codon, \@STOP_CODONS)){
	        return $i;
            }
        }
    }
    return -1;
}

# param 0 -> TrancscriptVariationAllele object.
# return the sequence of the transcript with the mutation, including 5' and 3'.
sub get_variation_cdna_seq{
    my $tva = $_[0];
    # seq contains 5' and 3' regions.
    my $seq = $tva->transcript->seq->seq;
    # Variation position counting utr regions.
    my $variation_start = $tva->transcript_variation->cdna_start - 1;
    my $variation_end = $tva->transcript_variation->cdna_end - 1;
    # If is a deletion, feature_seq is '-', so we will use '' instead
    # to build the final sequence.
    my $feature_seq = $tva->feature_seq eq "-" ? "" : $tva->feature_seq;
    substr($seq, $variation_start, $variation_end - $variation_start + 1) = $feature_seq;
    
    return $seq;
}

# param 0 -> TrancscriptVariationAllele object.
# return the coding sequence of the transcript with the mutation.
sub get_variation_cds_seq{
    my $tva = $_[0];
    # translateable_seq returns the coding part of the transcript
    # (it removes introns and 5' and 3' utr)
    my $seq = $tva->transcript->translateable_seq;
    # Variation position starting at the begining of coding sequence.
    my $variation_start = $tva->transcript_variation->cds_start - 1;
    my $variation_end = $tva->transcript_variation->cds_end - 1;
    # If is a deletion, feature_seq is '-', so we will use '' instead
    # to build the final sequence.
    my $feature_seq = $tva->feature_seq eq "-" ? "" : $tva->feature_seq;
    substr($seq, $variation_start, $variation_end - $variation_start + 1) = $feature_seq;
    
    return $seq;
}

# Generates the final ORF sequence of the mutated transcript.
# This sequence is stored by the following path:
# data/transcript id/variation id/alleleid.txt
# param 0 -> TranscriptVariationObject needed to build the path.
# param 1 -> The ORF final sequence.
sub generate_variation_seq_files{
    my $tva = $_[0];
    my $seq = $_[1];
    my $transcript_id = $tva->transcript->display_id;
    my $variation_id = $tva->transcript_variation->variation_feature->variation_name;
    my $transcript_variation_allele_id = $tva->dbID;
    my $path = '../data/sequences/' . $transcript_id . '/' . $variation_id . '/' ;
    make_path($path);
    my $file_path = $path . $transcript_variation_allele_id . '.txt';
    open(my $fd, '>', $file_path) or die "Could not open file '$file_path' $!";
    print $fd $seq;
    close($fd);
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

# Receives a transcript object and extract information about
# cds start and cds end errors.
# param 0 -> Transcript object
# returns a string with info about cds errors at 5' and 3'.
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

# Receives a Variation object and returns a string with the following format:
# phenotype description 1 (source 1) - phenotype description 2 (source 2) ...
# param 0 -> Variation object.
# returns a string with info about phenotype of the variation.
sub get_phenotype_info{
    my $variation = $_[0];
    my $pf_list = $variation->get_all_PhenotypeFeatures;
    my $pf_formatted = '';
    foreach my $pf (@{$pf_list}){
        $pf_formatted = $pf_formatted . $pf->phenotype->description . '(' . $pf->source_name . ')-';
    }
    chop($pf_formatted);
    return $pf_formatted;
}

# Build a string with publications about a variation and its links
# param 0-> Variation object.
# returns a string with format title1 -> link1 title2 -> link2...
sub get_publications_info{
    my $variation = $_[0];
    my $publication_list = $variation->get_all_Publications;
    my $publication_info = '';
    foreach my $publication (@{$publication_list}){
        $publication_info = $publication_info . get_publication_info($publication) . ' ';
    }
    return $publication_info;
}

# Receives a publication object and obtains its url.
# param 0 -> Publication object.
# returns a string with format title -> link.
sub get_publication_info{
    my $publication = $_[0];
    my $publication_info = '';
    my $publication_url;
    if ($publication->title){
        $publication_info = $publication->title;
    } else {
        $publication_info = "Title not provided"
    }
    $publication_info = $publication_info . " -> ";
    $publication_url = myUtils::Publication::find_link_by_doi($publication->doi) if ($publication->doi);
    $publication_url = myUtils::Publication::find_link_by_pmid($publication->pmid) if (!$publication_url && $publication->pmid);
    $publication_url = myUtils::Publication::find_link_by_pmcid($publication->pmcid) if (!$publication_url && $publication->pmcid);
    $publication_url = 'url not provided' if (!$publication_url);

    $publication_info = $publication_info . $publication_url;
    return $publication_info;
}











