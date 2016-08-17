setwd("/home/fabad/hemodonacion/src")

# Leer el csv sin filtros
csv = read.csv("final_out_no_filter.csv", sep="\t",stringsAsFactors=FALSE)
csv[,"MUTATED_SEQUENCE_LENGTH"]=as.numeric(gsub("%","",csv$MUTATED_SEQUENCE_LENGTH))
csv[,"KOZAK_MUTATED_SEQUENCE_LENGTH"]=as.numeric(gsub("%","",csv$KOZAK_MUTATED_SEQUENCE_LENGTH))
csv = add_signal_lost_sup_info(csv)

csv$READING_FRAME_STATUS = factor(csv$READING_FRAME_STATUS)
csv$KOZAK_READING_FRAME_STATUS = factor(csv$KOZAK_READING_FRAME_STATUS)
csv$KOZAK_STOP_CODON = factor(csv$KOZAK_STOP_CODON)
csv$STOP_CODON_POSITION = factor(csv$STOP_CODON_POSITION)
csv$GENE_NAME = factor(csv$GENE_NAME)
csv$SIGNAL_FIRST_MET_AFFECTED = factor(csv$SIGNAL_FIRST_MET_AFFECTED)
csv$SIGNAL_FIRST_KOZAK_AFFECTED = factor(csv$SIGNAL_FIRST_KOZAK_AFFECTED)




# Eliminar los casos en los que hay errores en las regiones 5' o 3'
csv = csv[csv$CDS_ERRORS == '',]

# Eliminar los casos en los que el biotipo es "non_stop_decay" o "nonsense_mediated_decay"
csv = csv[csv$TRANSCRIPT_BIOTYPE != 'non_stop_decay' & csv$TRANSCRIPT_BIOTYPE != 'nonsense_mediated_decay',]

# Numero de genes con mutaciones afectando en el codon de inicio
length(unique(csv$GENE_NAME))

# Obtener un conjunto en el que existe MAF definida y otro en el que no.
csvWithMaf = csv[!is.na(csv$MINOR_ALLELE_FREQUENCY),]
csvNoMaf = csv[is.na(csv$MINOR_ALLELE_FREQUENCY),]

# Dividir el subconjunto con MAF definida en MAF alta y baja
highMaf = csvWithMaf[csvWithMaf$MINOR_ALLELE_FREQUENCY >= 0.01,]
lowMaf = csvWithMaf[csvWithMaf$MINOR_ALLELE_FREQUENCY < 0.01,]

# Resumen de cada subconjunto de datos segun la maf
summary(lowMaf$FIRST_MET_POSITION)
summary(lowMaf$SIGNAL_FIRST_MET_AFFECTED)
summary(lowMaf$READING_FRAME_STATUS)
# Marco de lectura conservado si, ademas de tener "Conserved" tiene una longitud mayor al 1% de la seq original
nrow(lowMaf[lowMaf$MUTATED_SEQUENCE_LENGTH > 1 & lowMaf$READING_FRAME_STATUS == 'Conserved',])
nrow(lowMaf[lowMaf$MUTATED_SEQUENCE_LENGTH <= 1 & lowMaf$READING_FRAME_STATUS == 'Conserved',]) + nrow(lowMaf[lowMaf$READING_FRAME_STATUS == 'Lost',])
summary(lowMaf$SIGNAL_FIRST_MET_AFFECTED)

summary(highMaf$FIRST_MET_POSITION)
summary(highMaf$SIGNAL_FIRST_MET_AFFECTED)
summary(highMaf$READING_FRAME_STATUS)
nrow(highMaf[highMaf$MUTATED_SEQUENCE_LENGTH > 1 & highMaf$READING_FRAME_STATUS == 'Conserved',])
nrow(highMaf[highMaf$MUTATED_SEQUENCE_LENGTH <= 1 & highMaf$READING_FRAME_STATUS == 'Conserved',]) + nrow(highMaf[highMaf$READING_FRAME_STATUS == 'Lost',])
summary(highMaf$SIGNAL_FIRST_MET_AFFECTED)

summary(lowMaf$KOZAK_START)
summary(lowMaf$SIGNAL_FIRST_KOZAK_AFFECTED)
summary(lowMaf$KOZAK_READING_FRAME_STATUS)
nrow(lowMaf[lowMaf$KOZAK_MUTATED_SEQUENCE_LENGTH > 1 & lowMaf$KOZAK_READING_FRAME_STATUS == 'Conserved',])
nrow(lowMaf[lowMaf$KOZAK_MUTATED_SEQUENCE_LENGTH <= 1 & !is.na(lowMaf$KOZAK_MUTATED_SEQUENCE_LENGTH) & lowMaf$KOZAK_READING_FRAME_STATUS == 'Conserved',]) + nrow(lowMaf[lowMaf$KOZAK_READING_FRAME_STATUS == 'Lost',])
summary(lowMaf$SIGNAL_FIRST_KOZAK_AFFECTED)

summary(highMaf$KOZAK_START)
summary(highMaf$SIGNAL_FIRST_KOZAK_AFFECTED)
summary(highMaf$KOZAK_READING_FRAME_STATUS)
nrow(highMaf[highMaf$KOZAK_MUTATED_SEQUENCE_LENGTH > 1 & highMaf$KOZAK_READING_FRAME_STATUS == 'Conserved',])
nrow(highMaf[highMaf$KOZAK_MUTATED_SEQUENCE_LENGTH <= 1 & highMaf$KOZAK_READING_FRAME_STATUS == 'Conserved',]) + nrow(highMaf[highMaf$KOZAK_READING_FRAME_STATUS == 'Lost',])
summary(highMaf$SIGNAL_FIRST_KOZAK_AFFECTED)

# Histogramas de la posición de la primera metionina
hist(highMaf$FIRST_MET_POSITION, xlim = c(0,1000))
hist(lowMaf$FIRST_MET_POSITION, xlim = c(0,1000))

# Histogramas de la posición de la metionina de la primera secuencia Kozak
# con puntuacion mayor a 25%
hist(highMaf$KOZAK_START)
hist(lowMaf$KOZAK_START)

# Comprobar homogeneidad de varianzas
var.test(highMaf$FIRST_MET_POSITION, lowMaf$FIRST_MET_POSITION) # Varianzas distintas
var.test(highMaf$KOZAK_START, lowMaf$KOZAK_START) # Varianzas distintas

#Test de normalidad
shapiro.test(highMaf$FIRST_MET_POSITION) # No normal
shapiro.test(lowMaf$FIRST_MET_POSITION) # No normal
shapiro.test(highMaf$KOZAK_START) # No normal
shapiro.test(lowMaf$KOZAK_START) # No normal

# Test de wilcoxon para comparar medias
wilcox.test(highMaf$FIRST_MET_POSITION, lowMaf$FIRST_MET_POSITION, paired = F, conf.level = 0.95) # Distribuciones diferentes
wilcox.test(highMaf$KOZAK_START, lowMaf$KOZAK_START, paired = F, conf.level = 0.95) # Distribuciones diferentes

# Boxplots
boxplot(highMaf$FIRST_MET_POSITION, lowMaf$FIRST_MET_POSITION, ylim=c(0,3000))
boxplot(highMaf$KOZAK_START, lowMaf$KOZAK_START, ylim=c(0,1500))
op <- par(mfrow = c(1, 2))
boxplot(highMaf$FIRST_MET_POSITION, lowMaf$FIRST_MET_POSITION, ylab="Initial codon position",names=c("High MAF", "Low MAF"),
        main="First initial codon position comparative\nbetween high and low MAF")
boxplot(highMaf$KOZAK_START, lowMaf$KOZAK_START, ylim=c(0,4000), ylab="Initial codon position", names = c("High MAF", "Low MAF"),
        main="First initial codon position in strong Kozak sequence\ncomparative between high and low MAF.")
par(op)


###
  #borrar
met_cons=length(lowMaf[lowMaf$MUTATED_SEQUENCE_LENGTH > 1 & lowMaf$READING_FRAME_STATUS == 'Conserved',])
met_nocons=length(lowMaf[lowMaf$MUTATED_SEQUENCE_LENGTH <= 1 & lowMaf$READING_FRAME_STATUS == 'Conserved',]) + length(lowMaf[lowMaf$READING_FRAME_STATUS == 'Lost',])


# select variables v1, v2, v3
myvars <- c("TRANSCRIPT_ID","FIRST_MET_POSITION","STOP_CODON_POSITION","MUTATED_SEQUENCE_LENGTH", "KOZAK_START", "KOZAK_END", "KOZAK_MUTATED_SEQUENCE_LENGTH", "KOZAK_STOP_CODON", "KOZAK_READING_FRAME_STATUS")
myvars=c("MUTATED_SEQ_LENGTH2", "MUTATED_SEQUENCE_LENGTH")
View(csv[myvars])
myvars = c("CHROMOSOME", "GENE_ID", "GENE_NAME", "TRANSCRIPT_ID", "TRANSCRIPT_REFSEQ_ID", "TRANSCRIPT_BIOTYPE", "CDS_ERRORS", "PROTEIN_ID", "VARIATION_NAME", "SOURCE", "TRANSCRIPT_VARIATION_ALLELE_DBID", "MINOR_ALLELE_FREQUENCY", "CODON_CHANGE", "AMINOACID_CHANGE", "FIRST_MET_POSITION", "STOP_CODON_POSITION", "MUTATED_SEQUENCE_LENGTH", "READING_FRAME_STATUS", "KOZAK_START", "KOZAK_END", "KOZAK_STOP_CODON", "KOZAK_MUTATED_SEQUENCE_LENGTH", "KOZAK_ORF_AA_LENGTH", "KOZAK_IDENTITY", "KOZAK_RELIABILITY", "KOZAK_READING_FRAME_STATUS", "KOZAK_PROTEIN_SEQ", "SIGNAL_PEPTIDE_START", "SIGNAL_PEPTIDE_END", "CONSEQUENCE", "PHENOTYPE", "SO_TERM", "SIFT", "POLYPHEN", "PUBLICATIONS")
add_signal_lost_sup_info = function(csv){
  for (i in 1:nrow(csv)){
    signal_start = csv[i,"SIGNAL_PEPTIDE_START"]
    # Sumamos dos ya que en el csv se indica el inicio del ultimo codon
    # De esta forma, signal_end tiene el nucleotido exacto donde
    # termina el peptido señal. 
    signal_end = (csv[i,"SIGNAL_PEPTIDE_END"]) + 2
    first_met = csv[i,"FIRST_MET_POSITION"]
    first_kozak = csv[i,"KOZAK_START"]
    signal_first_met_affected = NA
    signal_first_kozak_affected = NA
    if (!is.na(signal_start) && !is.na(signal_end)){
      if(!is.na(first_met)){
        if (first_met == signal_start){signal_first_met_affected = "Totally conserved"}
        if (first_met > signal_start && first_met <= signal_end){signal_first_met_affected = "Partially conserved"}
        if (first_met == (signal_end + 1)){signal_first_met_affected = "Exactly lost"}
        if (first_met > (signal_end + 1)){signal_first_met_affected = "Lost"}
      }
      if (!is.na(first_kozak)){
        if (first_kozak == signal_start){signal_first_kozak_affected = "Totally conserved"}
        if (first_kozak > signal_start && first_kozak <= signal_end){signal_first_kozak_affected = "Partially conserved"}
        if (first_kozak == (signal_end + 1)){signal_first_kozak_affected = "Exactly lost"}
        if (first_kozak > (signal_end + 1)){signal_first_kozak_affected = "Lost"}
      }
    }
    csv[i, "SIGNAL_FIRST_MET_AFFECTED"] = signal_first_met_affected
    csv[i, "SIGNAL_FIRST_KOZAK_AFFECTED"] = signal_first_kozak_affected
  }
  return(csv)
}

add_kozak_mutated_seq_length = function(csv){
  for (i in 1:nrow(csv)){
  #for (i in 1:3){
    met_start = csv[i, "FIRST_MET_POSITION"]
    met_end = csv[i, "STOP_CODON_POSITION"]
    met_seq_length = as.double(met_end) - as.double(met_start) + 1
    met_percentage = as.numeric(gsub("%","",csv[i, "MUTATED_SEQUENCE_LENGTH"]))
    #met_percentage = as.double(csv[i, "MUTATED_SEQUENCE_LENGTH"])
    
    kozak_start = csv[i, "KOZAK_START"]
    kozak_end = csv[i, "KOZAK_END"]
    kozak_seq_length =as.double(kozak_end) - as.double(kozak_start) + 1
    #cat("\n")
    #cat("met start=",met_start,"\n", "met end=",met_end, "\n", "met length=",met_seq_length,"\n")
    #cat("kozak start=",kozak_start,"\n", "kozak end=",kozak_end, "\n", "kozak length=",kozak_seq_length,"\n")
    
    if(!is.na(kozak_seq_length) && !is.na(met_percentage) && !is.na(met_seq_length)){
      kozak_mutated_seq_length = kozak_seq_length * met_percentage / met_seq_length
      #cat("kozak_mutated_seq_length=", kozak_mutated_seq_length,"\t", "met_percentage=",met_percentage,"\n")
      csv[i, "KOZAK_MUTATED_SEQUENCE_LENGTH"] = paste(kozak_mutated_seq_length, "%", sep="")
    }else{
      csv[i, "KOZAK_MUTATED_SEQUENCE_LENGTH"] = NA
    }
    
  }
  return(csv)
}
write.table(csv[myvars], file = "final_out_no_filter2.csv", na="", sep="\t", row.names = F)
