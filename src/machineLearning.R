setwd("~/hemodonacion/data/tsv")


# Leer el csv sin filtros
# csv = read.csv("26_04_2017.csv", sep="\t",stringsAsFactors=FALSE)
csv = read.csv("9_12_2017.csv", sep="\t",stringsAsFactors=FALSE)
csv[,"MUTATED_SEQUENCE_LENGTH_1"]=as.numeric(gsub("%","",csv$MUTATED_SEQUENCE_LENGTH_1))
csv[,"MUTATED_SEQUENCE_LENGTH_2"]=as.numeric(gsub("%","",csv$MUTATED_SEQUENCE_LENGTH_2))
csv[,"MUTATED_SEQUENCE_LENGTH_3"]=as.numeric(gsub("%","",csv$MUTATED_SEQUENCE_LENGTH_3))
csv[,"SIGNAL_PEPTIDE_CONSERVATION_1"]=as.numeric(gsub("%","",csv$SIGNAL_PEPTIDE_CONSERVATION_1))
csv[,"SIGNAL_PEPTIDE_CONSERVATION_2"]=as.numeric(gsub("%","",csv$SIGNAL_PEPTIDE_CONSERVATION_2))
csv[,"SIGNAL_PEPTIDE_CONSERVATION_3"]=as.numeric(gsub("%","",csv$SIGNAL_PEPTIDE_CONSERVATION_3))

for(i in 1:nrow(csv)){
  if(csv$READING_FRAME_STATUS_1[i] == ''){
    csv$READING_FRAME_STATUS_1[i] = NA
  }
  if(csv$READING_FRAME_STATUS_2[i] == ''){
    csv$READING_FRAME_STATUS_2[i] = NA
  }
  if(csv$READING_FRAME_STATUS_3[i] == ''){
    csv$READING_FRAME_STATUS_3[i] = NA
  }
}
csv$READING_FRAME_STATUS_1 = factor(csv$READING_FRAME_STATUS_1)
csv$READING_FRAME_STATUS_2 = factor(csv$READING_FRAME_STATUS_2)
csv$READING_FRAME_STATUS_3 = factor(csv$READING_FRAME_STATUS_3)
csv$CODON_CHANGE = factor(csv$CODON_CHANGE)
csv$AMINOACID_CHANGE = factor(csv$AMINOACID_CHANGE)
# Calcular metioninas en 5' utr
csv$NMETS_5_UTR = sapply(strsplit(csv$METS_IN_5_UTR, " "), length)

# Contar los que mantienen fase de lectura
for (i in 1:nrow(csv)) {
  metUtrCol = csv$METS_IN_5_UTR[i]
  if(is.na(metUtrCol) || metUtrCol == ""){
    csv$CONSERVED_METS_IN_5_UTR[i] = 0
  } else {
    csv$CONSERVED_METS_IN_5_UTR[i] = length(grep("conserved", strsplit(metUtrCol," ")[[1]]))
  }
}
#Contar las que pierden la fase de lectura
for (i in 1:nrow(csv)) {
  metUtrCol = csv$METS_IN_5_UTR[i]
  if(is.na(metUtrCol) || metUtrCol == ""){
    csv$LOST_METS_IN_5_UTR[i] = 0
  } else {
    csv$LOST_METS_IN_5_UTR[i] = length(grep("lost", strsplit(metUtrCol," ")[[1]]))
  }
}
myvars=c('NMETS_5_UTR','METS_IN_5_UTR', 'CONSERVED_METS_IN_5_UTR', 'LOST_METS_IN_5_UTR')
View(csv[myvars])


# Insertar clase
for (i in 1:nrow(csv)) {
  if(is.na(csv$MINOR_ALLELE_FREQUENCY[i]) || !is.numeric(csv$MINOR_ALLELE_FREQUENCY[i])){
    csv$CLASS[i] = NA
  } else if (csv$MINOR_ALLELE_FREQUENCY[i] >= 0.01){
    csv$CLASS[i] = 'BENIGN'
  } else if (csv$MINOR_ALLELE_FREQUENCY[i] < 0.01){
    csv$CLASS[i] = 'DELETERIOUS'
  }
}
csv$CLASS = factor(csv$CLASS)

# Ordenar y eliminar variables sin importancia
#myvars=c('CHROMOSOME','GENE_ID', 'GENE_NAME', 'TRANSCRIPT_ID', 'TRANSCRIPT_REFSEQ_ID', 'TRANSCRIPT_BIOTYPE', 'METS_IN_5_UTR', 'NMETS_5_UTR', 'CONSERVED_METS_IN_5_UTR', 'LOST_METS_IN_5_UTR', 'SIGNAL_PEPTIDE_START', 'SIGNAL_PEPTIDE_END', 'CDS_ERRORS', 'PROTEIN_ID', 'VARIATION_NAME', 'VARIATION_TYPE', 'SOURCE', 'TRANSCRIPT_VARIATION_ALLELE_DBID', 'MINOR_ALLELE_FREQUENCY', 'CODON_CHANGE', 'CDS_COORDS', 'AMINOACID_CHANGE' ,'MET_POSITION_1' ,'STOP_CODON_POSITION_1' ,'MUTATED_SEQUENCE_LENGTH_1' ,'READING_FRAME_STATUS_1', 'SIGNAL_PEPTIDE_CONSERVATION_1', 'MET_POSITION_2', 'STOP_CODON_POSITION_2', 'MUTATED_SEQUENCE_LENGTH_2', 'SCORE_2', 'READING_FRAME_STATUS_2', 'SIGNAL_PEPTIDE_CONSERVATION_2', 'MET_POSITION_3', 'INIT_CODON_3', 'STOP_CODON_POSITION_3', 'MUTATED_SEQUENCE_LENGTH_3', 'SCORE_3', 'READING_FRAME_STATUS_3', 'SIGNAL_PEPTIDE_CONSERVATION_3', 'CONSEQUENCE', 'PHENOTYPE', 'SO_TERM', 'SIFT', 'POLYPHEN', 'PUBLICATIONS')
myvars=c('NMETS_5_UTR', 'CONSERVED_METS_IN_5_UTR', 'LOST_METS_IN_5_UTR', 'SIGNAL_PEPTIDE_START', 'SIGNAL_PEPTIDE_END', 'CODON_CHANGE', 'AMINOACID_CHANGE' ,'MET_POSITION_1' ,'STOP_CODON_POSITION_1' ,'MUTATED_SEQUENCE_LENGTH_1' ,'READING_FRAME_STATUS_1', 'SIGNAL_PEPTIDE_CONSERVATION_1', 'MET_POSITION_2', 'STOP_CODON_POSITION_2', 'MUTATED_SEQUENCE_LENGTH_2', 'READING_FRAME_STATUS_2', 'SIGNAL_PEPTIDE_CONSERVATION_2', 'MET_POSITION_3', 'STOP_CODON_POSITION_3', 'MUTATED_SEQUENCE_LENGTH_3', 'READING_FRAME_STATUS_3', 'SIGNAL_PEPTIDE_CONSERVATION_3', 'CLASS')
csv = csv[myvars]

# Exportar a weka
library("RWeka")
write.arff(csv, "../weka/9_12_2017.arff")
