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


csv$READING_FRAME_STATUS_1 = factor(csv$READING_FRAME_STATUS_1)
csv$READING_FRAME_STATUS_2 = factor(csv$READING_FRAME_STATUS_2)
csv$READING_FRAME_STATUS_3 = factor(csv$READING_FRAME_STATUS_3)
csv$STOP_CODON_POSITION_1 = factor(csv$STOP_CODON_POSITION_1)
csv$STOP_CODON_POSITION_2 = factor(csv$STOP_CODON_POSITION_2)
csv$STOP_CODON_POSITION_3 = factor(csv$STOP_CODON_POSITION_3)
csv$GENE_NAME = factor(csv$GENE_NAME)
#csv$SIGNAL_PEPTIDE_CONSERVATION_1 = factor(csv$SIGNAL_PEPTIDE_CONSERVATION_1)
#csv$SIGNAL_PEPTIDE_CONSERVATION_2 = factor(csv$SIGNAL_PEPTIDE_CONSERVATION_2)
#csv$SIGNAL_PEPTIDE_CONSERVATION_3 = factor(csv$SIGNAL_PEPTIDE_CONSERVATION_3)


# Numero de genes con mutaciones afectando en el codon de inicio
length(unique(csv$GENE_NAME))
View(table(csv$GENE_NAME))
nrow(csv)/length(unique(csv$GENE_NAME))

# Numero de transcritos afectados
length(unique(csv$TRANSCRIPT_ID))

# Numero de variaciones
length(unique(csv$VARIATION_NAME))

# Numero de alelos distintos
length(unique(csv$TRANSCRIPT_VARIATION_ALLELE_DBID))

# Transcritos afectados de TP53, CACNA1C y CDKN2A
length(unique(csv[csv$GENE_NAME=="TP53",]$TRANSCRIPT_ID))
length(unique(csv[csv$GENE_NAME=="CACNA1C",]$TRANSCRIPT_ID))
length(unique(csv[csv$GENE_NAME=="CDKN2A",]$TRANSCRIPT_ID))
length(unique(csv[csv$GENE_NAME=="PAX5",]$TRANSCRIPT_VARIATION_ALLELE_DBID))
length(unique(csv[csv$GENE_NAME=="PAX5",]$TRANSCRIPT_ID))
length(unique(csv[csv$GENE_NAME=="PAX5",]$VARIATION_NAME))
length(unique(csv[csv$GENE_NAME=="DTNA",]$TRANSCRIPT_VARIATION_ALLELE_DBID))
length(unique(csv[csv$GENE_NAME=="DTNA",]$TRANSCRIPT_ID))
length(unique(csv[csv$GENE_NAME=="DTNA",]$VARIATION_NAME))

# Numero de deleciones, inserciones y cambios puntuales
View(table(csv$VARIATION_TYPE))

# El uso de la metionina en kozak fuerte provoca la conservacion del marco
# de lectura en mayor medida que la primera metionina encontrada?
rf_met_conserved = nrow(csv[csv$READING_FRAME_STATUS_1 == "Conserved",]); rf_met_conserved
rf_met_lost = nrow(csv[csv$READING_FRAME_STATUS_1 == "Lost",]); rf_met_lost
rf_kozak_conserved = nrow(csv[csv$READING_FRAME_STATUS_3 == "Conserved",]); rf_kozak_conserved
rf_kozak_lost = nrow(csv[csv$READING_FRAME_STATUS_3 == "Lost",]); rf_kozak_lost
m = as.table(rbind(c(rf_met_conserved, rf_met_lost), c(rf_kozak_conserved, rf_kozak_lost)))
dimnames(m) = list(alt_met=c("First Met", "Kozak Met"),
                   reading_frame=c("Conserved", "Lost"))
m
chisq.test(m) # Muy significativo. El uso de la metionina en una secuencia de kozak fuerte favorece el mantenimiento del marco de lectura.

# Obtener un conjunto en el que existe MAF definida y otro en el que no.
csvWithMaf = csv[!is.na(csv$MINOR_ALLELE_FREQUENCY),]
csvNoMaf = csv[is.na(csv$MINOR_ALLELE_FREQUENCY),]
# Numero de variantes
length(unique(csvWithMaf$VARIATION_NAME))
length(unique(csvNoMaf$VARIATION_NAME))


# Dividir el subconjunto con MAF definida en MAF alta y baja
highMaf = csvWithMaf[csvWithMaf$MINOR_ALLELE_FREQUENCY >= 0.01,]
lowMaf = csvWithMaf[csvWithMaf$MINOR_ALLELE_FREQUENCY < 0.01,]
# Numero de variantes
length(unique(highMaf$VARIATION_NAME))
length(unique(lowMaf$VARIATION_NAME))
# Numero de transcritos afectados
length(unique(highMaf$TRANSCRIPT_ID))
length(unique(lowMaf$TRANSCRIPT_ID))
# Numero de transcritos sin peptido senal detectado
length(unique(highMaf[is.na(highMaf$SIGNAL_PEPTIDE_START),]$TRANSCRIPT_ID))
length(unique(lowMaf[is.na(lowMaf$SIGNAL_PEPTIDE_START),]$TRANSCRIPT_ID))

# Resumen de cada subconjunto de datos segun la maf
summary(lowMaf$MET_POSITION_1)
summary(lowMaf$SIGNAL_PEPTIDE_CONSERVATION_1)
summary(lowMaf$READING_FRAME_STATUS_1)


summary(highMaf$MET_POSITION_1)
summary(highMaf$SIGNAL_PEPTIDE_CONSERVARION_1)
summary(highMaf$READING_FRAME_STATUS_1)

summary(lowMaf$MET_POSITION_2)
summary(lowMaf$SIGNAL_PEPTIDE_CONSERVARION_2)
summary(lowMaf$READING_FRAME_STATUS_2)

summary(highMaf$MET_POSITION_2)
summary(highMaf$SIGNAL_PEPTIDE_CONSERVARION_2)
summary(highMaf$READING_FRAME_STATUS_2)

summary(highMaf$MET_POSITION_3)
summary(highMaf$SIGNAL_PEPTIDE_CONSERVARION_3)
summary(highMaf$READING_FRAME_STATUS_3)

summary(lowMaf$MET_POSITION_3)
summary(lowMaf$SIGNAL_PEPTIDE_CONSERVARION_3)
summary(lowMaf$READING_FRAME_STATUS_3)


# Histogramas de la posición de la primera metionina
hist(highMaf$MET_POSITION_1, xlim = c(0,1000))
hist(lowMaf$MET_POSITION_1, xlim = c(0,1000))

# Histogramas de la posición de la metionina de la primera secuencia Kozak
# con puntuacion mayor a 25%
hist(highMaf$MET_POSITION_2)
hist(lowMaf$MET_POSITION_2)

# Comprobar homogeneidad de varianzas
var.test(highMaf$MET_POSITION_1, lowMaf$MET_POSITION_1) # Varianzas distintas
var.test(highMaf$MET_POSITION_2, lowMaf$MET_POSITION_2) # Varianzas distintas
var.test(highMaf$MET_POSITION_3, lowMaf$MET_POSITION_3) # Varianzas distintas

#Test de normalidad
shapiro.test(highMaf$MET_POSITION_1) # No normal
shapiro.test(lowMaf$MET_POSITION_1) # No normal
shapiro.test(highMaf$MET_POSITION_2) # No normal
shapiro.test(lowMaf$MET_POSITION_2) # No normal

# Test de wilcoxon para comparar medias
wilcox.test(highMaf$MET_POSITION_1, lowMaf$MET_POSITION_1, paired = F, conf.level = 0.95) # Distribuciones diferentes
wilcox.test(highMaf$MET_POSITION_2, lowMaf$MET_POSITION_2, paired = F, conf.level = 0.95) # Distribuciones diferentes
wilcox.test(highMaf$MET_POSITION_3, lowMaf$MET_POSITION_3, paired = F, conf.level = 0.95) # Distribuciones diferentes

# Boxplots
boxplot(highMaf$MET_POSITION_1, lowMaf$MET_POSITION_1, ylim=c(0,3000))
boxplot(highMaf$MET_POSITION_2, lowMaf$MET_POSITION_2, ylim=c(0,1500))

op <- par(mfrow = c(1, 2))
boxplot(highMaf$MET_POSITION_1, lowMaf$MET_POSITION_1,  ylab="First AUG position (bp)",names=c("High MAF", "Low MAF"))
boxplot(highMaf$MET_POSITION_2, lowMaf$MET_POSITION_2, ylim=c(0,600), ylab="Posición del codón inicial (en pares de bases)", names = c("MAF alta", "MAF baja"),
        main="Comparativa de la posición del primer codón de inicio\npredicho por ATGpr con score > 0.25\nentre los grupos de MAF alta y baja.")
boxplot(highMaf$MET_POSITION_3, lowMaf$MET_POSITION_3, ylab="First AUG position in strong Kozak (bp)", names = c("High MAF", "Low MAF"))

par(op)

# Test de Chi Cuadrado para comparar las variables cualitativas
# MET1 READING FRAME STATUS
rf_conserved_met1_low_maf = nrow(lowMaf[lowMaf$READING_FRAME_STATUS_1 == 'Conserved',])
rf_lost_met1_low_maf = nrow(lowMaf[lowMaf$READING_FRAME_STATUS_1 == 'Lost',])
rf_conserved_met1_high_maf = nrow(highMaf[highMaf$READING_FRAME_STATUS_1 == 'Conserved',])
rf_lost_met1_high_maf = nrow(highMaf[highMaf$READING_FRAME_STATUS_1 == 'Lost',])
m = as.table(rbind(c(rf_conserved_met1_low_maf,rf_lost_met1_low_maf), c(rf_conserved_met1_high_maf,rf_lost_met1_high_maf)))
dimnames(m)=list(GRUPO = c("MAF BAJA", "MAF ALTA"),
                 READING_FRAME_STATUS = c("CONSERVED", "LOST"))
m
chisq.test(m) # p-value menor que 0.05 indica que cada grupo de mutaciones presenta diferencias significativas
###

# MET2 READING FRAME STATUS
rf_conserved_met2_low_maf = nrow(lowMaf[lowMaf$READING_FRAME_STATUS_2 == 'Conserved',])
rf_lost_met2_low_maf = nrow(lowMaf[lowMaf$READING_FRAME_STATUS_2 == 'Lost',])
rf_conserved_met2_high_maf = nrow(highMaf[highMaf$READING_FRAME_STATUS_2 == 'Conserved',])
rf_lost_met2_high_maf = nrow(highMaf[highMaf$READING_FRAME_STATUS_2 == 'Lost',])
m = as.table(rbind(c(rf_conserved_met2_low_maf,rf_lost_met2_low_maf), c(rf_conserved_met2_high_maf,rf_lost_met2_high_maf)))
dimnames(m)=list(GRUPO = c("MAF BAJA", "MAF ALTA"),
                 READING_FRAME_STATUS = c("CONSERVED", "LOST"))
m
chisq.test(m) # p-value mayor que 0.05 indica que cada grupo de mutaciones no presenta diferencias significativas

# MET3 READING FRAME STATUS
rf_conserved_met3_high_maf = nrow(highMaf[highMaf$READING_FRAME_STATUS_3 == 'Conserved',])
rf_lost_met3_high_maf = nrow(highMaf[highMaf$READING_FRAME_STATUS_3 == 'Lost',])
rf_conserved_met3_low_maf = nrow(lowMaf[lowMaf$READING_FRAME_STATUS_3 == 'Conserved',])
rf_lost_met3_low_maf = nrow(lowMaf[lowMaf$READING_FRAME_STATUS_3 == 'Lost',])
m = as.table(rbind(c(rf_conserved_met3_low_maf,rf_lost_met3_low_maf), c(rf_conserved_met3_high_maf,rf_lost_met3_high_maf)))
dimnames(m)=list(GRUPO = c("MAF BAJA", "MAF ALTA"),
                 READING_FRAME_STATUS = c("CONSERVED", "LOST"))
m
chisq.test(m) # p-value mayor que 0.05 indica que cada grupo de mutaciones no presenta diferencias significativas

# SIGNAL PEPTIDE AFFECTED MET1
summary(highMaf$SIGNAL_PEPTIDE_CONSERVATION_1)
length(highMaf$SIGNAL_PEPTIDE_CONSERVATION_1)
summary(lowMaf$SIGNAL_PEPTIDE_CONSERVATION_1)
length(lowMaf$SIGNAL_PEPTIDE_CONSERVATION_1)
# SIGNAL PEPTIDE AFFECTED MET2
summary(highMaf$SIGNAL_PEPTIDE_CONSERVATION_2)
length(highMaf$SIGNAL_PEPTIDE_CONSERVATION_2)
summary(lowMaf$SIGNAL_PEPTIDE_CONSERVATION_2)
length(lowMaf$SIGNAL_PEPTIDE_CONSERVATION_2)
# SIGNAL PEPTIDE AFFECTED MET3
summary(highMaf$SIGNAL_PEPTIDE_CONSERVATION_3)
length(highMaf$SIGNAL_PEPTIDE_CONSERVATION_3)
summary(lowMaf$SIGNAL_PEPTIDE_CONSERVATION_3)
length(lowMaf$SIGNAL_PEPTIDE_CONSERVATION_3)

# PEPTIDE SIGNAL BOXPLOTS
op <- par(mfrow = c(1, 2))
boxplot(highMaf$SIGNAL_PEPTIDE_CONSERVATION_1, lowMaf$SIGNAL_PEPTIDE_CONSERVATION_1, ylab="Conservation percentage of the signal peptide",names=c("High MAF", "Low MAF"),main="A", ylim=c(0,120))
boxplot(highMaf$SIGNAL_PEPTIDE_CONSERVATION_2, lowMaf$SIGNAL_PEPTIDE_CONSERVATION_2, ylab="Conservacion del peptido señal (en %)", names = c("MAF alta", "MAF baja"),
        main="Comparativa de la conservacion del peptido señal\ncon la met de ATGpr con score > 0.25\nentre los grupos de MAF alta y baja.")
boxplot(highMaf$SIGNAL_PEPTIDE_CONSERVATION_3, lowMaf$SIGNAL_PEPTIDE_CONSERVATION_3, ylab="Conservation percentage of the signal peptide", names = c("High MAF", "Low MAF"),main="B", ylim=c(0,120))
par(op)

# SIGNAL PEPTIDE COMPARISON
wilcox.test(highMaf$SIGNAL_PEPTIDE_CONSERVATION_1, lowMaf$SIGNAL_PEPTIDE_CONSERVATION_1, paired = F, conf.level = 0.95) # Distribuciones diferentes
wilcox.test(highMaf$SIGNAL_PEPTIDE_CONSERVATION_2, lowMaf$SIGNAL_PEPTIDE_CONSERVATION_2, paired = F, conf.level = 0.95) # Distribuciones diferentes
wilcox.test(highMaf$SIGNAL_PEPTIDE_CONSERVATION_3, lowMaf$SIGNAL_PEPTIDE_CONSERVATION_3, paired = F, conf.level = 0.95) # Distribuciones diferentes

# POSICIONES AFECTADAS
# Nos quedamos con las mutaciones que solo afectan a una posicion.
affectedPosHighMaf = c()
for (i in 1:(nrow(highMaf))){
  if (is.na(highMaf$CDS_COORDS[i]) || highMaf$CDS_COORDS[i] == ""){
    next
  }
  interval = strsplit(highMaf$CDS_COORDS[i], "-")
  begin = as.numeric(interval[[1]][[1]])
  end = as.numeric(interval[[1]][[2]])
  if(begin == end){
    affectedPosHighMaf = rbind(affectedPosHighMaf, cbind(highMaf[i,], AFFECTED_POS = begin))
  }
}

affectedPosLowMaf = c()
for (i in 1:(nrow(lowMaf))){
  if (is.na(lowMaf$CDS_COORDS[i]) || lowMaf$CDS_COORDS[i] == ""){
    next
  }
  interval = strsplit(lowMaf$CDS_COORDS[i], "-")
  begin = as.numeric(interval[[1]][[1]])
  end = as.numeric(interval[[1]][[2]])
  if(begin == end){
    affectedPosLowMaf = rbind(affectedPosLowMaf, cbind(lowMaf[i,], AFFECTED_POS = begin))
  }
}

# Numero de variantes
length(unique(affectedPosHighMaf$VARIATION_NAME))
length(unique(affectedPosLowMaf$VARIATION_NAME))
# Numero de transcritos afectados
length(unique(affectedPosHighMaf$TRANSCRIPT_ID))
length(unique(affectedPosLowMaf$TRANSCRIPT_ID))

tableAffectedPosHM= table(affectedPosHighMaf$AFFECTED_POS)
tableAffectedPosLM = table(affectedPosLowMaf$AFFECTED_POS)
m = as.table(rbind(c(tableAffectedPosHM[1], tableAffectedPosHM[2], tableAffectedPosHM[3]), c(tableAffectedPosLM[1], tableAffectedPosLM[2], tableAffectedPosLM[3])))
dimnames(m)=list(GROUP = c("MAF BAJA", "MAF ALTA"),
                 AFFECTED_CODON = c("A", "T", "G"))
m
chisq.test(m)


# CODONES DE INICIO PREVIOS AL CODON DE INICIO ORIGINAL
op <- par(mfrow = c(2,1))
numberOfMet5UTRLowMAF=sapply(strsplit(lowMaf$METS_IN_5_UTR, " "), length)
plot(density(numberOfMet5UTRLowMAF), col=gray(0.4), main="", xlab = "Number of \"AUG\" codons in 5' UTR")
numberOfMet5UTRHighMAF=sapply(strsplit(highMaf$METS_IN_5_UTR, " "), length)
lines(density(numberOfMet5UTRHighMAF), col =gray(0), lty=5)
legend(34, 0.68, c("Low MAF","High MAF"), lty=c(1,5), col = c(gray(0.4),gray(0)))
summary(highMaf$METS_IN_5_UTR)
summary(lowMaf$METS_IN_5_UTR)

wilcox.test(numberOfMet5UTRLowMAF, numberOfMet5UTRHighMAF, paired = F, conf.level = 0.95)

# Contar los que mantienen fase de lectura
for (i in 1:nrow(highMaf)) {
  metUtrCol = highMaf$METS_IN_5_UTR[i]
  if(is.na(metUtrCol) || metUtrCol == ""){
    highMaf$CONSERVED_METS_IN_5_UTR[i] = 0
  } else {
    highMaf$CONSERVED_METS_IN_5_UTR[i] = length(grep("conserved", strsplit(metUtrCol," ")[[1]]))
  }
}
for (i in 1:nrow(lowMaf)) {
  metUtrCol = lowMaf$METS_IN_5_UTR[i]
  if(is.na(metUtrCol) || metUtrCol == ""){
    lowMaf$CONSERVED_METS_IN_5_UTR[i] = 0
  } else {
    lowMaf$CONSERVED_METS_IN_5_UTR[i] = length(grep("conserved", strsplit(metUtrCol," ")[[1]]))
  }
}

summary(highMaf$CONSERVED_METS_IN_5_UTR)
summary(lowMaf$CONSERVED_METS_IN_5_UTR)
wilcox.test(highMaf$CONSERVED_METS_IN_5_UTR, lowMaf$CONSERVED_METS_IN_5_UTR, paired = F, conf.level = 0.95)
plot(density(lowMaf$CONSERVED_METS_IN_5_UTR), col=gray(0.4), main="B", xlab = "Number of \"AUG\" codons maintaining the reading frame in 5' UTR")
lines(density(highMaf$CONSERVED_METS_IN_5_UTR), col =gray(0), lty=5)
legend(14, 1, c("Low MAF","High MAF"), lty=c(1,5), col = c(gray(0.4),gray(0)))
par(op)

View(highMaf[highMaf$METS_IN_5_UTR,highMaf$CONSERVED_METS_IN_5_UTR])
myvars = c("METS_IN_5_UTR", "CONSERVED_METS_IN_5_UTR")
View(highMaf[myvars])
# select variables v1, v2, v3
myvars <- c("TRANSCRIPT_ID","FIRST_MET_POSITION","STOP_CODON_POSITION","MUTATED_SEQUENCE_LENGTH", "KOZAK_START", "KOZAK_END", "KOZAK_MUTATED_SEQUENCE_LENGTH", "KOZAK_STOP_CODON", "KOZAK_READING_FRAME_STATUS")
myvars=c("TRANSCRIPT_ID", "VARIATION_NAME", "VARIATION_TYPE", "CDS_COORDS", "AFFECTED_POS")
View(affectedPosLowMaf[myvars])
View(affectedPosHighMaf[myvars])
myvars = c("CHROMOSOME", "GENE_ID", "GENE_NAME", "TRANSCRIPT_ID", "TRANSCRIPT_REFSEQ_ID", "TRANSCRIPT_BIOTYPE", "CDS_ERRORS", "PROTEIN_ID", "VARIATION_NAME", "SOURCE", "TRANSCRIPT_VARIATION_ALLELE_DBID", "MINOR_ALLELE_FREQUENCY", "CODON_CHANGE", "AMINOACID_CHANGE", "FIRST_MET_POSITION", "STOP_CODON_POSITION", "MUTATED_SEQUENCE_LENGTH", "READING_FRAME_STATUS", "KOZAK_START", "KOZAK_END", "KOZAK_STOP_CODON", "KOZAK_MUTATED_SEQUENCE_LENGTH", "KOZAK_ORF_AA_LENGTH", "KOZAK_IDENTITY", "KOZAK_RELIABILITY", "KOZAK_READING_FRAME_STATUS", "KOZAK_PROTEIN_SEQ", "SIGNAL_PEPTIDE_START", "SIGNAL_PEPTIDE_END", "CONSEQUENCE", "PHENOTYPE", "SO_TERM", "SIFT", "POLYPHEN", "PUBLICATIONS")

nrow(csv[csv$GENE_NAME=="AADACL3",])
nrow(csv[csv$GENE_NAME=="C6orf7",])
View(csv[csv$GENE_NAME=="GP6",])
nrow(csv[csv$GENE_NAME=="LENG9",])
View(csv[csv$GENE_NAME=="NARS2",])
View(csv[csv$GENE_NAME=="NUDT11",])
View(csv[csv$GENE_NAME=="SAAL1",])
View(csv[csv$GENE_NAME=="TOP3B",])
View(csv[csv$GENE_NAME=="ZNF682",])
View(csv[csv$GENE_NAME=="ZNF827",])


myvars=c("FIRST_MET_POSITION","STOP_CODON_POSIION","MUTATED_SEQUENCE_LENGTH", "KOZAK_START", "KOZAK_END", "KOZAK_MUTATED_SEQUENCE_LENGTH")
View(csv[myvars])

m = as.table(rbind(c(tableAffectedPosHM[1], tableAffectedPosHM[2], tableAffectedPosHM[3]), c(tableAffectedPosLM[1], tableAffectedPosLM[2], tableAffectedPosLM[3])))
dimnames(m)=list(GROUP = c("MAF BAJA", "MAF ALTA"),
                 AFFECTED_CODON = c("A", "T", "G"))
m
chisq.test(m)
