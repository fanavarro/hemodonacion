package myUtils::EnsemblUtils;
use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use myUtils::CsvManager;
use myUtils::Publication;
use myUtils::KozakUtils;
use myUtils::SeqUtils;
use myUtils::PhobiusService;
use File::Path qw(make_path);
use List::Util qw[min max];
use DateTime;
use myUtils::MatchPWM;
use myUtils::Noderer;
use base 'Class::Singleton';

my $MAX_KOZAK_RESULTS = 50000;

sub _new_instance {
	my $class = shift;
	my $this  = bless { }, $class;
	$this->{registry} = getRegistry();
	$this->{transcript_adaptor} = $this->{registry}->get_adaptor( 'homo_sapiens', 'core', 'transcript' );
	$this->{slice_adaptor} = $this->{registry}->get_adaptor( 'Human', 'Core', 'Slice' );
	$this->{trv_adaptor} = $this->{registry}->get_adaptor( 'homo_sapiens', 'variation', 'transcriptvariation' );
	return $this;
}

sub getRegistry{
	my $this = shift;
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
	return $registry;
}

# Apply a filter to the TranscriptVariation list passed as argument.
# param 0 -> Reference list of transcripts.
sub filter_transcript_variations{
	my $this = shift;
    my @original_list = @{$_[1]};
    my @filtered_list;
    my $cont = 1;
    my $max = scalar (@original_list);
    my $start_time = DateTime->now;
    foreach my $transcript_variation (@original_list){
        print $cont . "/" . $max . "\n";
        $cont = $cont + 1;
        if (pass_filter($transcript_variation)){
            push @filtered_list, $transcript_variation;
        }
    }
    my $end_time = DateTime->now;
    my $elapse = $end_time - $start_time;
    print "Elapsed time : ".$elapse->in_units('minutes')."m\n";
    return \@filtered_list;
}
# Functions that receives a TranscriptVariation object
# and returns if it pass the filter or not.
sub pass_filter{
	my $this = shift;
    my $transcript_variation = shift;
    my $transcript = $transcript_variation->transcript;
    my $variation_feature = $transcript_variation->variation_feature;

    # Checkings at transcript level:

    # If it contains something in cds_error, transcript does not pass the filter.
    my $cds_errors = get_cds_errors($transcript);
    if ($cds_errors ne ""){
        print $transcript->stable_id . " with CDS errors.\n";
        return 0;
    }

    # If transcript does not have translation, return false.
    if (!defined($transcript->translation)){
        #print $transcript->stable_id . " with no translation.\n";
        return 0;
    }

    my @evidence_status = @{$variation_feature->get_all_evidence_values};
    if (scalar(@evidence_status) == 0){
        #print $variation_feature->display_id . " with no evidences.\n";
        return 0;
    }
    return 1;
}

# Get TranscriptVariation objects from transcript list and so term list.
# param 1 -> transcript reference list
# param 2 -> SO term reference list
# return a list ref of hashes with information about transcript variation objects in param 1.
sub get_transcript_variations_by_transcripts_SO_terms {
	my $this = shift;
    my $transcripts = shift;
    my $so_terms = shift;
    return $this->{trv_adaptor}->fetch_all_by_Transcripts_SO_terms($transcripts, $so_terms);

}

# Get TranscriptVariation objects from transcripts and variation position.
# param 0 -> Reference list of transcripts.
# param 1 -> Position of the variation in amino acid sequence.
# return Reference list of transcript variations.
sub get_transcript_variations_by_transcripts_peptide_position{
	my $this = shift;
    my $transcripts =shift;
    my $var_pos = shift;
    my $constraint = "(translation_start=$var_pos or translation_end=$var_pos)";
    return $this->{trv_adaptor}->fetch_all_by_Transcripts_with_constraint($transcripts, $constraint);
}


# Get info about TranscriptVariation list.
# param 0 -> Chromosome
# param 1 -> TranscriptVariation ref list object.
# return a list ref of hashes with information about transcript variation objects in param 1.
sub get_transcript_variation_info{
	my $this = shift;
    my $chromosome = $_[0];
    my $trvs = $_[1];

    # List where we are adding all the info about the transcript variation allele objects.
    my @result_list;

    my $count = 1;
    my $n = scalar @{$trvs};
    foreach my $tv ( @{$trvs} ) {
        print "Chromosome $chromosome\t Transcript Variation: $count/$n\t" . $tv->transcript->display_id . '-' . $tv->variation_feature->variation_name . "\n";
        $count = $count + 1;
	
        my $variation = $tv->variation_feature->variation;
        my $transcript = $tv->transcript;
        my $cds_errors = get_cds_errors($transcript);
        # If there exist errors in cds or variation does not have evidences, we skip the transcript_variation.
        if ($cds_errors ne ""){next};
        if (scalar(@{$variation->get_all_evidence_values}) == 0){next};
        my $minor_allele_frequency = $tv->variation_feature->minor_allele_frequency ? $tv->variation_feature->minor_allele_frequency : '';
        my $phenotype_info = get_phenotype_info($variation);
        my $publications_info = get_publications_info($variation);
        my $ref_seq_mrna_ids = get_ref_seq_mrna_ids($transcript);
        my $signal_peptide_info = get_signal_peptide_info($transcript);
        my $cds_coords = '[' . $tv->cds_start . ', ' . $tv->cds_end . ']' if (defined $tv->cds_start && defined $tv->cds_end) ;
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
            my $atgpr_info = get_atgpr_info($tva);
            my $pwm_info = get_match_pwm_info($tva);
            my $noderer_info = get_noderer_info($tva);
            my $mets_in_5_utr = get_five_prime_utr_info($tva);

            $entry{'CHROMOSOME'} = $chromosome;
            $entry{'GENE_ID'} = $transcript->get_Gene->stable_id;
            $entry{'GENE_NAME'} = $transcript->get_Gene->external_name;
            $entry{'TRANSCRIPT_ID'} = $transcript->display_id;
            $entry{'TRANSCRIPT_REFSEQ_ID'} = $ref_seq_mrna_ids;
            $entry{'TRANSCRIPT_BIOTYPE'} = $transcript->biotype;
            $entry{'METS_IN_5_UTR'} = $mets_in_5_utr;
            $entry{'CDS_ERRORS'} = $cds_errors;
            $entry{'PROTEIN_ID'} = $transcript->translation->display_id;
            $entry{'VARIATION_NAME'} = $tv->variation_feature->variation_name;
            $entry{'VARIATION_TYPE'} = $tv->variation_feature->var_class;
            $entry{'SOURCE'} = $variation->source_name;
            $entry{'TRANSCRIPT_VARIATION_ALLELE_DBID'} = $tva->dbID;
            $entry{'MINOR_ALLELE_FREQUENCY'} = $minor_allele_frequency;
            $entry{'CODON_CHANGE'} = defined($tva->display_codon_allele_string) ? $tva->display_codon_allele_string : '';
            $entry{'CDS_COORDS'} = $cds_coords;
            $entry{'AMINOACID_CHANGE'} = defined($tva->pep_allele_string) ? $tva->pep_allele_string : '';
            
            # method 1 columns (first atg)
            $entry{'APPROACH1_MET_POSITION'} = $seq_info->{'met_position'};
            $entry{'APPROACH1_STOP_CODON_POSITION'} = $seq_info->{'stop_codon_position'};
            $entry{'APPROACH1_MUTATED_SEQUENCE_LENGTH'} = $seq_info->{'seq_length'};
            $entry{'APPROACH1_READING_FRAME_STATUS'} = $seq_info->{'reading_frame'};
            $entry{'APPROACH1_SIGNAL_PEPTIDE_CONSERVATION'} = get_signal_peptide_conservarion($signal_peptide_info->{'START'}, $signal_peptide_info->{'END'}, $seq_info->{'met_position'});
            
            # method 2 columns (noderer scores)
            $entry{'APPROACH2_MET_POSITION'} = $noderer_info->{'met_position'};
            $entry{'APPROACH2_INIT_CODON'} = $noderer_info->{'init_codon'};
            $entry{'APPROACH2_STOP_CODON_POSITION'} = $noderer_info->{'stop_codon_position'};
            $entry{'APPROACH2_MUTATED_SEQUENCE_LENGTH'} = $noderer_info->{'seq_length'};
            $entry{'APPROACH2_SCORE'} = $noderer_info->{'score'};
            $entry{'APPROACH2_READING_FRAME_STATUS'} = $noderer_info->{'reading_frame'};
            $entry{'APPROACH2_SIGNAL_PEPTIDE_CONSERVATION'} = get_signal_peptide_conservarion($signal_peptide_info->{'START'}, $signal_peptide_info->{'END'}, $noderer_info->{'met_position'});

            # Old method 2 columns (atgPR)
            # $entry{'MET_POSITION_2'} = $atgpr_info->{'START'};
            # $entry{'STOP_CODON_POSITION_2'} = $atgpr_info->{'STOP_CODON'} eq 'Yes' ? $atgpr_info->{'FINISH'} : '';
            # $entry{'MUTATED_SEQUENCE_LENGTH_2'} = $atgpr_info->{'KOZAK_MUTATED_SEQUENCE_LENGTH'};
            # $entry{'SCORE_2'} = $atgpr_info->{'RELIABILITY'};
            # $entry{'READING_FRAME_STATUS_2'} = $atgpr_info->{'FRAMESHIFT'};
            # $entry{'SIGNAL_PEPTIDE_CONSERVATION_2'} = get_signal_peptide_conservarion($signal_peptide_info->{'START'}, $signal_peptide_info->{'END'}, $atgpr_info->{'START'});

            # method 3 columns (Kozak pwm matches)
            $entry{'APPROACH3_MET_POSITION'} = $pwm_info->{'met_position'};
            $entry{'APPROACH3_INIT_CODON'} = $pwm_info->{'init_codon'};
            $entry{'APPROACH3_STOP_CODON_POSITION'} = $pwm_info->{'stop_codon_position'};
            $entry{'APPROACH3_MUTATED_SEQUENCE_LENGTH'} = $pwm_info->{'seq_length'};
            $entry{'APPROACH3_SCORE'} = $pwm_info->{'score'};
            $entry{'APPROACH3_READING_FRAME_STATUS'} = $pwm_info->{'reading_frame'};
            $entry{'APPROACH3_SIGNAL_PEPTIDE_CONSERVATION'} = get_signal_peptide_conservarion($signal_peptide_info->{'START'}, $signal_peptide_info->{'END'}, $pwm_info->{'met_position'});

            $entry{'SIGNAL_PEPTIDE_START'} = $signal_peptide_info->{'START'} if (defined($signal_peptide_info));
            $entry{'SIGNAL_PEPTIDE_END'} = $signal_peptide_info->{'END'} if (defined($signal_peptide_info));
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
            push(@result_list, \%entry);
        }
    }
    return \@result_list;
}

# Checks if an element exists in the given list.
# param 0 -> Element to check.
# param 1 -> List reference.
# returns true if elements exists inside the list and false if not.
sub exists_in_list {
	my $this = shift;
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
# param 1 -> Constraints to pass to DB.
# return list reference of transcripts.
sub get_transcripts_by_chromosomes_with_constraint {
	my $this = shift;
    my $chromosomes = $_[0];
    my $constraints = $_[1];
    my @transcripts = ();
    while ( my $chromosome = shift @{$chromosomes} ) {
        my $aux = get_transcripts_by_chromosome_with_constraints($chromosome, $constraints);
        push( @transcripts, @{$aux} );
    }
    return \@transcripts;
}

# Get all transcripts that belongs to chromosome.
# param 0 -> Chromosome name.
# param 1 -> Constraints to pass to DB.
# return list reference of transcripts.
sub get_transcripts_by_chromosome_with_constraints {
	my $this = shift;
    my $chromosome = $_[0];
    my $constraints = $_[1];
    my $slice = $this->{slice_adaptor}->fetch_by_region( 'chromosome', $chromosome );
    my $transcripts = $this->{transcript_adaptor}->fetch_all_by_Slice_constraint($slice, $constraints);
    return $transcripts;
}


# Extract information about variation sequence
# when we have information about alleles, that it is
# to say, when the source from the tva is dbsnp.
# param 0 -> TranscriptVariationAllele object.
# return -> Hash with the following keys:
# 'first_met_position': position of first Met found relative to peptide.
# If there is a deletion, initial affected codon is filled with 5' intron
# for example AAAA|C(A)T|GATTAGCACC -> AA(A)AGATTAGCACC where () indicate
# reference to count positions.
# 'reading_frame': Conserved or Lost.
# 'stop_codon_position': Position of the first stop codon (keeping reading frame).
# 'seq_length': Percentage of the orf found in the mutated sequence size according to original orf.
sub get_sequence_info{
	my $this = shift;
    my $tva = $_[0];
    my $hash_seq_info = {};
    $hash_seq_info->{'met_position'} = '';
    $hash_seq_info->{'reading_frame'} = '';
    $hash_seq_info->{'stop_codon_position'} = '';
    $hash_seq_info->{'seq_length'} = '';

    my $mutated_cdna = get_variation_cdna_seq($tva);
    if (defined($mutated_cdna)){
        my $cdna = $tva->transcript->seq->seq;
        my $cds = $tva->transcript->translateable_seq;
        my $five_prime_affected = !defined($tva->transcript_variation->cds_start);
        $hash_seq_info = myUtils::SeqUtils::get_met_mutation_info($cdna, $cds, $mutated_cdna, $five_prime_affected);
    }

   
    return $hash_seq_info;

}

# param 0 -> TrancscriptVariationAllele object.
# return the sequence of the transcript with the mutation, including 5' and 3'.
sub get_variation_cdna_seq{
	my $this = shift;
    my $tva = $_[0];
    # seq contains 5' and 3' regions.
    my $seq = $tva->transcript->seq->seq;
    if (!defined($tva->transcript_variation->cdna_start) || !defined($tva->transcript_variation->cdna_end)){
        print "ERROR CDNA variation without start or end " . $tva->transcript_variation->variation_feature->variation_name . " " . $tva->transcript->display_id . "\n";
        return undef;
    }
    # Variation position counting utr regions.
    my $variation_start = $tva->transcript_variation->cdna_start - 1;
    my $variation_end = $tva->transcript_variation->cdna_end - 1;
    # If is a deletion, feature_seq is '-', so we will use '' instead
    # to build the final sequence.
    my $feature_seq = $tva->feature_seq eq "-" ? "" : $tva->feature_seq;

    # If feature seq not contains nucleotides or it is very large, possible error on data...
    if ( $feature_seq =~ m/[^ATGCNatgcn]/ || length($feature_seq) > length($seq)){
        print "Feature Seq not available -> " . $feature_seq . "\n";
        return undef;
    }
    substr($seq, $variation_start, $variation_end - $variation_start + 1) = $feature_seq;
    
    
    return $seq;
}

# param 0 -> TrancscriptVariationAllele object.
# return the coding sequence of the transcript with the mutation.
sub get_variation_cds_seq{
	my $this = shift;
    my $tva = $_[0];
    # translateable_seq returns the coding part of the transcript
    # (it removes introns and 5' and 3' utr)
    my $seq = $tva->transcript->translateable_seq;
    if (!defined($tva->transcript_variation->cds_end) || !defined($tva->transcript_variation->cds_start)){
        print "ERROR variation end or start not defined inside CDS sequence. " . $tva->transcript_variation->variation_feature->variation_name . " " . $tva->transcript->display_id . "\n";
        return undef;
    }
    
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
	my $this = shift;
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
	my $this = shift;
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
	my $this = shift;
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
	my $this = shift;
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
	my $this = shift;
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
	my $this = shift;
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


# Return the ref seq mrna identifiers
# from a transcript. If there are several ids,
# they will be separated by '-'.
# param 0 -> Transcript object.
# return string with ref seq mrna identifiers.
sub get_ref_seq_mrna_ids{
	my $this = shift;
    my $transcript = $_[0];
    my @db_entries = @{$transcript->get_all_DBEntries('RefSeq_mRNA')};
    my $transcript_refseq_id = '';
    foreach my $dbentry ( @db_entries ){
        if ($dbentry->display_id){
            $transcript_refseq_id = $transcript_refseq_id . $dbentry->display_id . '-';
        }
    }
    chop($transcript_refseq_id);
    return $transcript_refseq_id;
}

# Get the first sequence in the
# original transcript and in the mutated
# transcript by using atgpr prediction of init codon.
# This information is used
# to determine if the mutated transcript
# has lost the frameshift.
# param 0 -> TranscriptVariationAllele object
# return -> Hash with atgpr predicted sequence info in
# the mutated transcript.
sub get_atgpr_info{
	my $this = shift;
    my $tva = $_[0];
    my $hash_atgpr_info = {};
    $hash_atgpr_info->{'FRAMESHIFT'} = '';
    $hash_atgpr_info->{'PREVIOUS_ATGS'} = '';
    $hash_atgpr_info->{'RELIABILITY'} = '';
    $hash_atgpr_info->{'KOZAK_IDENTITY'} = '';
    $hash_atgpr_info->{'START'} = '';
    $hash_atgpr_info->{'FINISH'} = '';
    $hash_atgpr_info->{'ORF_AMINOACID_LENGTH'} = '';
    $hash_atgpr_info->{'STOP_CODON'} = '';
    $hash_atgpr_info->{'PROTEIN_SEQUENCE'} = '';
    $hash_atgpr_info->{'KOZAK_MUTATED_SEQUENCE_LENGTH'} = '';
   


    my $mutated_cdna_seq = get_variation_cdna_seq($tva);
    if (!defined($mutated_cdna_seq)){
        return $hash_atgpr_info;
    }
    my $original_cdna_seq = $tva->transcript->seq->seq;
    my $original_cds_seq = $tva->transcript->translateable_seq;

    # Find the position in which coding region starts at cdna sequence.
    my $default_kozak_pos = myUtils::SeqUtils::get_translation_start_pos($original_cdna_seq, $original_cds_seq);
    my $mutated_kozaks =  myUtils::KozakUtils::get_kozak_info($mutated_cdna_seq, $MAX_KOZAK_RESULTS);
    my $original_kozaks =  myUtils::KozakUtils::get_kozak_info($original_cdna_seq, $MAX_KOZAK_RESULTS);
    
    # We start counting positions from reference
    my $reference;

    # If deletion affecting 5' region, we move reference
    if(length($original_cdna_seq) > length($mutated_cdna_seq) && !defined($tva->transcript_variation->cds_start)){
        $reference = max($default_kozak_pos - (length($original_cdna_seq) - length($mutated_cdna_seq)), 0);
    } else {
        $reference = $default_kozak_pos;
    }

    # Look for the natural kozak
    my $natural_kozak = ();
    foreach my $original_kozak (@{$original_kozaks}){
        if ($original_kozak->{'START'} == $default_kozak_pos){
            $natural_kozak = $original_kozak;
            last;
        }
    }

    # Correction in order to point to the original sequence position.
    # If muation is a insertion, we have to add values to reference.
    my $position_correction = 0;
    if (length($mutated_cdna_seq) > length($original_cdna_seq)){
        $position_correction = -(length($mutated_cdna_seq) - length($original_cdna_seq));
     } 
     if (length($mutated_cdna_seq) < length($original_cdna_seq) && defined($tva->transcript_variation->cds_start)){
        $position_correction = (length($original_cdna_seq) - length($mutated_cdna_seq));
     } 

    # Look for the first kozak after default kozak in mutated.
    my $first_mutated_kozak = ();
    foreach my $mutated_kozak (@{$mutated_kozaks}){
        # If we are iterating outside coding region, break the loop. Kozak not found.
        if (max($mutated_kozak->{'START'} - $reference + $position_correction, 0) >= length($original_cds_seq)){
            last;
        }
        if ($mutated_kozak->{'START'} >= $reference && $mutated_kozak->{'RELIABILITY'} >= 0.25){
            $first_mutated_kozak = $mutated_kozak;
            last;
        }
    }

    # If natural kozak or first kozak in mutated sequence are not been found, return empty hash.
    if (!defined($natural_kozak->{'START'}) || !defined($first_mutated_kozak->{'START'})){
        return $hash_atgpr_info;
    }



    # Calculate frameshift with SeqUtils
    my $mutated_orf = myUtils::SeqUtils::get_orf($mutated_cdna_seq, $first_mutated_kozak->{'START'});   
    $hash_atgpr_info->{'FRAMESHIFT'} = myUtils::SeqUtils::is_in_frame($mutated_orf, $original_cds_seq) ? 'Conserved' : 'Lost';

    $hash_atgpr_info->{'PREVIOUS_ATGS'} = $first_mutated_kozak->{'PREVIOUS_ATGS'};
    $hash_atgpr_info->{'RELIABILITY'} = $first_mutated_kozak->{'RELIABILITY'};
    $hash_atgpr_info->{'KOZAK_IDENTITY'} = $first_mutated_kozak->{'KOZAK_IDENTITY'};
    # If mutated and original met are different, apply the pos correction
    if($first_mutated_kozak->{'START'} != $natural_kozak->{'START'}){
        # We perform the substraction to get the position in cds sequence
        $hash_atgpr_info->{'START'} = max($first_mutated_kozak->{'START'} - $reference + $position_correction, 0);
    } else {
        $hash_atgpr_info->{'START'} = max($first_mutated_kozak->{'START'} - $reference, 0);
    }
    # We add +1 to point at the first base of the stop codon
    $hash_atgpr_info->{'FINISH'} = $first_mutated_kozak->{'FINISH'} - $reference + 1 + $position_correction;
    $hash_atgpr_info->{'ORF_AMINOACID_LENGTH'} = $first_mutated_kozak->{'ORF_AMINOACID_LENGTH'};
    $hash_atgpr_info->{'STOP_CODON'} = $first_mutated_kozak->{'STOP_CODON'};
    $hash_atgpr_info->{'PROTEIN_SEQUENCE'} = $first_mutated_kozak->{'PROTEIN_SEQUENCE'};
    $hash_atgpr_info->{'KOZAK_MUTATED_SEQUENCE_LENGTH'} = (length($mutated_orf) * 100 / length($original_cds_seq)) . '%';
    # Check if kozak position are not pointing to met.
    if (substr($original_cds_seq, $hash_atgpr_info->{'START'}, 3) ne 'ATG'){
        print ("Error\nfirst kozak pos = " . $hash_atgpr_info->{'START'} . " in $original_cds_seq\nCodon found = " . substr($original_cds_seq, $hash_atgpr_info->{'START'}, 3) . "\n");
    }
    return $hash_atgpr_info;
}

sub get_match_pwm_info{
	my $this = shift;
    my $tva = $_[0];
    my $hash_pwm_info = {};
    $hash_pwm_info->{'met_position'} = '';
    $hash_pwm_info->{'reading_frame'} = '';
    $hash_pwm_info->{'stop_codon_position'} = '';
    $hash_pwm_info->{'seq_length'} = '';
    $hash_pwm_info->{'init_codon'} = '';
    $hash_pwm_info->{'score'} = '';

    my $mutated_cdna = get_variation_cdna_seq($tva);
    if (defined($mutated_cdna)){
        my $cdna = $tva->transcript->seq->seq;
        my $cds = $tva->transcript->translateable_seq;
        my $translation_start_pos = myUtils::SeqUtils::get_translation_start_pos($cdna, $cds);
        my $five_prime_affected = !defined($tva->transcript_variation->cds_start);

        my $matchpwm = myUtils::MatchPWM->instance();
        my $hits = $matchpwm->myUtils::MatchPWM::get_kozak_matches($mutated_cdna);
        my $n = 0;
        if (defined $hits && defined $hits->{INIT_CODON_POS}){
            $n = scalar @{$hits->{INIT_CODON_POS}};
        }
        my $index = 0;
        # Buscamos un match que no este en 5' utr y que enga un ATG como codon de inicio
        while ($index < $n){
            if ($hits->{INIT_CODON_POS}->[$index] >= $translation_start_pos && substr($mutated_cdna, $hits->{INIT_CODON_POS}->[$index], 3) eq "ATG"){
                last;
            }
            $index = $index + 1;
        }
        # If we have found a kozak sequence after original met
        if($index < $n){
           my $found_met_pos = $hits->{INIT_CODON_POS}->[$index];
           #print "pwm met found at pos $found_met_pos (" . substr($cdna, $found_met_pos, 3) . ")\n";
           my $hash_seq_info = myUtils::SeqUtils::get_met_mutation_info($cdna, $cds, $mutated_cdna, $five_prime_affected, $found_met_pos);
           $hash_pwm_info->{'met_position'} = $hash_seq_info->{'met_position'};
           $hash_pwm_info->{'reading_frame'} = $hash_seq_info->{'reading_frame'};
           $hash_pwm_info->{'stop_codon_position'} = $hash_seq_info->{'stop_codon_position'};
           $hash_pwm_info->{'seq_length'} = $hash_seq_info->{'seq_length'};
           $hash_pwm_info->{'score'} = $hits->{SCORE}->[$index];
           $hash_pwm_info->{'init_codon'} = substr($mutated_cdna, $found_met_pos, 3);
        }
    }
    return $hash_pwm_info;
}

sub get_noderer_info {
	my $this = shift;
    my $tva = $_[0];
    my $hash_noderer_info = {};
    $hash_noderer_info->{'met_position'} = '';
    $hash_noderer_info->{'reading_frame'} = '';
    $hash_noderer_info->{'stop_codon_position'} = '';
    $hash_noderer_info->{'seq_length'} = '';
    $hash_noderer_info->{'init_codon'} = '';
    $hash_noderer_info->{'score'} = '';

    my $mutated_cdna = get_variation_cdna_seq($tva);
    if (defined($mutated_cdna)){
        my $cdna = $tva->transcript->seq->seq;
        my $cds = $tva->transcript->translateable_seq;
        my $translation_start_pos = myUtils::SeqUtils::get_translation_start_pos($cdna, $cds);
        my $five_prime_affected = !defined($tva->transcript_variation->cds_start);

        my $noderer = myUtils::Noderer->instance();
        my $min_efficiency = 87; # Kozak sequence lower bound.
        my $hits = $noderer->myUtils::Noderer::get_kozak_matches($mutated_cdna, $min_efficiency);
        my $n = 0;
        if (defined $hits && defined $hits->{INIT_CODON_POS}){
            $n = scalar @{$hits->{INIT_CODON_POS}};
        }
        my $index = 0;
        # We are looking for a kozak match downstream 5' utr.
        while ($index < $n){
            if ($hits->{INIT_CODON_POS}->[$index] >= $translation_start_pos){
                last;
            }
            $index = $index + 1;
        }
        # If we have found a kozak sequence after original met
        if($index < $n){
           my $found_met_pos = $hits->{INIT_CODON_POS}->[$index];
           my $hash_seq_info = myUtils::SeqUtils::get_met_mutation_info($cdna, $cds, $mutated_cdna, $five_prime_affected, $found_met_pos);
           $hash_noderer_info->{'met_position'} = $hash_seq_info->{'met_position'};
           $hash_noderer_info->{'reading_frame'} = $hash_seq_info->{'reading_frame'};
           $hash_noderer_info->{'stop_codon_position'} = $hash_seq_info->{'stop_codon_position'};
           $hash_noderer_info->{'seq_length'} = $hash_seq_info->{'seq_length'};
           $hash_noderer_info->{'score'} = $hits->{SCORE}->[$index];
           $hash_noderer_info->{'init_codon'} = substr($mutated_cdna, $found_met_pos, 3);
        }
    }
    return $hash_noderer_info;
}

# Receives a transcript object and return information
# about the predicted signal peptide using phobius.
# param 0 -> transcript object.
# return a hash reference with 'START' and 'END' keys
# indicating the position in which the signal peptide starts
# and ends in nucleotid coordinates starting by 0.
sub get_signal_peptide_info{
	my $this = shift;
    my $transcript = shift;
    my %signal_peptide_info = ();
    my $protein_seq = $transcript->translate;
    if (defined($protein_seq)){
        my %phobius_input = ($transcript->display_id => $protein_seq->seq);
        my $phobius_output = myUtils::PhobiusService::get_info_signal_peptide_local(\%phobius_input);
        if (defined($phobius_output)){
            my $feature_list = $phobius_output->{$transcript->display_id};
            foreach my $feature (@{$feature_list}){
                if ($feature->{'TYPE'} eq 'SIGNAL'){
                    # Translate from aminoacid to nucleotid coordinates.
                    my $signal_start = ($feature->{'START'} - 1) * 3;
                    my $signal_end = (($feature->{'END'} - 1) * 3) + 2;
                    $signal_peptide_info{'START'} = $signal_start;
                    $signal_peptide_info{'END'} = $signal_end;
                    last;
                }
            }
        }
    }

    # Return hash reference or undef if hash is empty.
    if (%signal_peptide_info){
        return \%signal_peptide_info;
    } else {
        return undef;
    }
}


# Get the percentage of signal peptide conservation
# by using an alternative initiation codon.
# param 0 -> signal peptide start in bp.
# param 1 -> signal peptide end in bp.
# param 2 -> alternative initiation codon position in bp.
sub get_signal_peptide_conservarion{
	my $this = shift;
    my $sp_start = shift;
    my $sp_end = shift;
    my $alt_met_pos = shift;

    if (!defined ($sp_start) || $sp_start eq '' || !defined($sp_end) || $sp_end eq '' || !defined($alt_met_pos) || $alt_met_pos eq ''){
        return '';
    }

    if ($alt_met_pos > $sp_end){
        return '0%';
    }

    my $sp_length = $sp_end - $sp_start + 1;
    my $sp_mutated_length = $sp_end - $alt_met_pos + 1;
    my $conservation = $sp_mutated_length * 100 / $sp_length;
    return $conservation . '%';
}

# Get info about mets in 5' utr region.
# param 0 -> transcript object
# return a formatted String with the following format:
#   met_position1_conserved|lost met_position2_conserved|lost ...
sub get_five_prime_utr_info{
	my $this = shift;
    my $formatted_result = '';
    my $tva = shift;
    my $transcript = $tva->transcript;
    my $mutated_seq = get_variation_cdna_seq($tva);
    my $wild_seq = $transcript->seq->seq;
    if (defined ($transcript->five_prime_utr) && defined ($wild_seq) && defined ($mutated_seq)){
        my $number_of_new_or_deleted_nucleotides = abs(length($wild_seq) - length($mutated_seq));
        my $met_info_list = myUtils::SeqUtils::get_5_utr_mets_info($transcript->five_prime_utr->seq, $number_of_new_or_deleted_nucleotides);
        foreach my $met_info (@{$met_info_list}){
            $formatted_result = $formatted_result . $met_info->{'met_position'} . '_';
            $formatted_result = $formatted_result . $met_info->{'reading_frame'} . '_';
            $formatted_result = $formatted_result . $met_info->{'termination_codon'} . ' ';
        }
        chop($formatted_result)
    }
    return $formatted_result;
}

1;
