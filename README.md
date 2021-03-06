# Lost in Translation
## Structure
This repository contains the scripts developed for the Master Thesis in Bioinformatics by the University of Murcia 'Lost in Translation: Bioinformatic analysis of variations affecting the translation initiation codon in the Human Genome.

The repository contains the following folders:
  * src -> Contains Perl and R scripts. It also has folder called 'myUtils' with used utilities. The principal scripts are:
    * 'ensemblMining.pl', whose function is to generate a TSV file retrieving data from ensembl about variants affecting the initiation codon and to calculate different features for each variant.
    * 'statistics.R', whose function is to analyse statistically the TSV file generated by 'ensemblMining.pl' by means of a comparison between mutations and polymorphisms.
    * 'exomes.R', whose function is to analyse excel files containing complete exomes of patients obtained from ION Torrent patform (not included, available upon request).
  * data/tsv -> Contains the TSV files obtained by 'ensemblMining.pl' in different dates.

## Dependencies
The following libraries are required to execute 'ensemblMining.pl' script:
  * Perl v5.22.1
  * Class::Singleton
  * Fcntl
  * FindBin
  * Getopt::Long
  * HTML::TableExtract
  * IPC::Open2
  * JSON
  * LWP::Simple
  * LWP::Simple::REST
  * List::Util
  * POSIX
  * SOAP::Lite
  * Statistics::R
  * Try::Tiny
  * XML::LibXML
  * base
  * strict
  * warnings
  * Bio::EnsEMBL::Registry
  * DateTime
  * File::Path
  * List::Util
  
Also, R must be installed in the system with bioconductor library and PWMEnrich package, which must be installed through bioconductor (https://www.bioconductor.org/install/).

## Running the script
ensemblMining.pl only has one parameter that indicates the name of the generated TSV file. You should use Perl v5.22.1:
```bash
perl ensemblMining.pl output.csv
```
## TSV format
Each row of the generated TSV files represents a TranscriptVariationAllele (TVA) object from the EnsEMBL api, which means a concrete allele of a concrete variation affecting a concrete transcript. For each of these Transcript Variation Allele objects, the following information is retrieved:

  * __CHROMOSOME__ -> The chromosome in which the tva is located.
  * __GENE_ID__ -> The EnsEMBL identifier of the affected gene.
  * __GENE_NAME__ -> The name of the affected gene.
  * __TRANSCRIPT_ID__ -> The EnsEMBL identifier of the affected transcript.
  * __TRANSCRIPT_REFSEQ_ID__ -> The RefSeq identifier of the affected transcript.
  * __TRANSCRIPT_BIOTYPE__ -> The annotated biotype of the transcript. All transcripts must have 'protein_coding' biotype due to 'ensemblMining' performs a filter to retrieve only this kind of transcript.
  * __METS_IN_5_UTR__ -> A list of strings with the format 'position_(lost|maintained)_(premature-termination|no-premature-termination)' separated by blank spaces. Each string indicates a position in 5' UTR of the wild transcript where an 'AUG' codon is found, if it maintains the original reading frame or not and if a termination codon is found in the same reading frame before the coding region. The position should be negative because the position 0 would be the start of the coding region.
  * __SIGNAL_PEPTIDE_START__ -> The position where the signal peptide starts in the coding region. If the signal peptide have not been recognised, the column will be empty. The scanning for the signal peptide is done by Phobius application (http://phobius.sbc.su.se/).
  * __SIGNAL_PEPTIDE_END__ -> The position where the signal peptide ends. If the signal peptide have not been recognised, the column will be empty.
  * __CDS_ERRORS__ -> The errors found in the CDS sequence. This value must be empty so that 'ensemblMining.pl' performs a filter before adding the entry.
  * __PROTEIN_ID__ -> The EnsEMBL identifier of the affected protein.
  * __VARIATION_NAME__ -> The variation identifier.
  * __VARIATION_TYPE__ -> The variation type (SNP, deletion, insertion...).
  * __SOURCE__ -> The origin of the variation (dbSNP, Cosmic, HGMD...)
  * __TRANSCRIPT_VARIATION_ALLELE_DBID__ -> The identifier of the tva object.
  * __MINOR_ALLELE_FREQUENCY__ -> The MAF associated with the variation.
  * __CODON_CHANGE__ -> A string indicating the change at codon level. For example 'atG/atA' indicates a simple nucleotide polymorphism that changes a G by an A.
  * __CDS_COORDS__ -> The range of position, in base pairs, affected by the variation. The format is '[startPos, endPos]' where the position 1 is the first nucleotide of the coding region of the transcript.
  * __AMINOACID_CHANGE__ -> A string indicating the change at amino acid level. For example, 'M/I' indicates that a methionine has been changed by an isoleucine.
  * __APPROACH1_MET_POSITION__ -> The position of the first 'AUG' codon found in the coding region of the mutated transcript, where position 0 would be occupied by the wild initiation codon. It could be an alternative initiation codon if the wild one has been lost.
  * __APPROACH1_STOP_CODON_POSITION__ -> The position of the first stop codon found in the same reading frame than 'APPROACH1_MET_POSITION'.
  * __APPROACH1_MUTATED_SEQUENCE_LENGTH__ -> The percentage of the open reading frame length provoked by the use of the 'AUG' codon found in the position 'APPROACH1_MET_POSITION' relative to the wild open reading frame.
  * __APPROACH1_READING_FRAME_STATUS__ -> Indicates if the reading frame of the 'AUG' codon located in 'APPROACH1_MET_POSITION' is the same than the wild initiation codon ('Maintained') or not ('Lost').
  * __APPROACH1_SIGNAL_PEPTIDE_CONSERVATION__ -> The percentage of signal peptide conservation by using the alternative initiation codon found in 'APPROACH1_MET_POSITION' if a signal peptide was detected in the transcript.
  * __APPROACH2_MET_POSITION__ -> The position of the first 'AUG' codon found in the coding region of the mutated transcript with an efficiency greater or equal to 87 according to Noderer work (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4299517/), where position 0 would be occupied by the wild initiation codon. It could be an alternative initiation codon if the wild one has been lost.
  * __APPROACH2_STOP_CODON_POSITION__-> The position of the first stop codon found in the same reading frame than 'APPROACH2_MET_POSITION'.
  * __APPROACH2_MUTATED_SEQUENCE_LENGTH__ -> The percentage of the open reading frame length provoked by the use of the 'AUG' codon found in the position 'APPROACH2_MET_POSITION' relative to the wild open reading frame.
  * __APPROACH2_SCORE__ -> The translation efficiency of the 'AUG' codon found in 'APPROACH2_MET_POSITION' according to Noderer work.
  * __APPROACH2_READING_FRAME_STATUS__ -> Indicates if the reading frame of the 'AUG' codon located in 'APPROACH2_MET_POSITION' is the same than the wild initiation codon ('Maintained') or not ('Lost').
  * __APPROACH2_SIGNAL_PEPTIDE_CONSERVATION__ -> The percentage of signal peptide conservation by using the alternative initiation codon found in 'APPROACH2_MET_POSITION' if a signal peptide was detected in the transcript.
  * __APPROACH3_MET_POSITION__ -> The position of the first 'AUG' codon found in the coding region of the mutated transcript by the scanning of Kozak context with a minimum score of 0.8 through PWMEnrich R library, where position 0 would be occupied by the wild initiation codon. It could be an alternative initiation codon if the wild one has been lost.
  * __APPROACH3_INIT_CODON__ -> The initiation codon found at 'APPROACH3_MET_POSITION'. It must be 'atg' because the script filter any other option.
  * __APPROACH3_STOP_CODON_POSITION__ -> The position of the first stop codon found in the same reading frame than 'APPROACH3_MET_POSITION'.
  * __APPROACH3_MUTATED_SEQUENCE_LENGTH__ -> The percentage of the open reading frame length provoked by the use of the 'AUG' codon found in the position 'APPROACH3_MET_POSITION' relative to the wild open reading frame.
  * __APPROACH3_SCORE__ -> The score of the Kozak region found by PWMEnrich.
  * __APPROACH3_READING_FRAME_STATUS__ -> Indicates if the reading frame of the 'aug' codon located at 'APPROACH3_MET_POSITION' is the same than the wild initiation codon ('Maintained') or not ('Lost').
  * __APPROACH3_SIGNAL_PEPTIDE_CONSERVATION__ -> The percentage of signal peptide conservation by using the alternative initiation codon found in 'APPROACH3_MET_POSITION' if a signal peptide was detected in the transcript.
  * __CONSEQUENCE__ -> The described consequences of the variation.
  * __PHENOTYPE__ -> The described phenotype of the variation. This contains a list of strings separated by '-' with the format 'phenotype description(source)'
  * __SO_TERM__ -> Sequence Ontology terms associated with the variation.
  * __SIFT__ -> Consequences of the variation predicted by SIFT.
  * __POLYPHEN__ -> Consequences of the variation predicted by PolyPhen.
  * __PUBLICATIONS__ -> List of strings with the format 'publication title -> url to article' separated with spaces. It contains publications associated with the variation.
