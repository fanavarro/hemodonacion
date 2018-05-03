setwd("~/hemodonacion/data/tsv")


# Leer el csv sin filtros
#csv = read.csv("26_04_2017.csv", sep="\t",stringsAsFactors=FALSE)
csv = read.csv("30_04_2018.tsv", sep="\t",stringsAsFactors=FALSE)
csv[,"APPROACH1_MUTATED_SEQUENCE_LENGTH"]=as.numeric(gsub("%","",csv$APPROACH1_MUTATED_SEQUENCE_LENGTH))
csv[,"APPROACH2_MUTATED_SEQUENCE_LENGTH"]=as.numeric(gsub("%","",csv$APPROACH2_MUTATED_SEQUENCE_LENGTH))
csv[,"APPROACH3_MUTATED_SEQUENCE_LENGTH"]=as.numeric(gsub("%","",csv$APPROACH3_MUTATED_SEQUENCE_LENGTH))
csv[,"APPROACH1_SIGNAL_PEPTIDE_CONSERVATION"]=as.numeric(gsub("%","",csv$APPROACH1_SIGNAL_PEPTIDE_CONSERVATION))
csv[,"APPROACH2_SIGNAL_PEPTIDE_CONSERVATION"]=as.numeric(gsub("%","",csv$APPROACH2_SIGNAL_PEPTIDE_CONSERVATION))
csv[,"APPROACH3_SIGNAL_PEPTIDE_CONSERVATION"]=as.numeric(gsub("%","",csv$APPROACH3_SIGNAL_PEPTIDE_CONSERVATION))

for(i in 1:nrow(csv)){
  if(csv$APPROACH1_READING_FRAME_STATUS[i] == ''){
    csv$APPROACH1_READING_FRAME_STATUS[i] = NA
  }
  if(csv$APPROACH2_READING_FRAME_STATUS[i] == ''){
    csv$APPROACH2_READING_FRAME_STATUS[i] = NA
  }
  if(csv$APPROACH3_READING_FRAME_STATUS[i] == ''){
    csv$APPROACH3_READING_FRAME_STATUS[i] = NA
  }
}
csv$APPROACH1_READING_FRAME_STATUS = factor(csv$APPROACH1_READING_FRAME_STATUS)
csv$APPROACH2_READING_FRAME_STATUS = factor(csv$APPROACH2_READING_FRAME_STATUS)
csv$APPROACH3_READING_FRAME_STATUS = factor(csv$APPROACH3_READING_FRAME_STATUS)
csv$APPROACH1_STOP_CODON_POSITION = factor(csv$APPROACH1_STOP_CODON_POSITION)
csv$APPROACH2_STOP_CODON_POSITION = factor(csv$APPROACH2_STOP_CODON_POSITION)
csv$APPROACH3_STOP_CODON_POSITION = factor(csv$APPROACH3_STOP_CODON_POSITION)
csv$GENE_NAME = factor(csv$GENE_NAME)
#csv$APPROACH1_SIGNAL_PEPTIDE_CONSERVATION = factor(csv$APPROACH1_SIGNAL_PEPTIDE_CONSERVATION)
#csv$APPROACH2_SIGNAL_PEPTIDE_CONSERVATION = factor(csv$APPROACH2_SIGNAL_PEPTIDE_CONSERVATION)
#csv$APPROACH3_SIGNAL_PEPTIDE_CONSERVATION = factor(csv$APPROACH3_SIGNAL_PEPTIDE_CONSERVATION)

# Calcular metioninas en 5' utr
csv$NMETS_5_UTR = sapply(strsplit(csv$METS_IN_5_UTR, " "), length)
# Contar los que mantienen fase de lectura y no existe codon de fin en el mismo reading frame
for (i in 1:nrow(csv)) {
  metUtrCol = csv$METS_IN_5_UTR[i]
  if(is.na(metUtrCol) || metUtrCol == ""){
    csv$CONSERVED_METS_IN_5_UTR[i] = 0
  } else {
    csv$CONSERVED_METS_IN_5_UTR[i] = length(grep("maintained_no-premature-termination", strsplit(metUtrCol," ")[[1]]))
  }
}
#Contar las que no puedan usarse como tis alternativo (reading frame perdido o reading frame mantenido con un codon de fin encontrado en el mismo)
for (i in 1:nrow(csv)) {
  metUtrCol = csv$METS_IN_5_UTR[i]
  if(is.na(metUtrCol) || metUtrCol == ""){
    csv$LOST_METS_IN_5_UTR[i] = 0
  } else {
    lostReadingFrame = length(grep("lost", strsplit(metUtrCol," ")[[1]]))
    maintainedReadingFrameWithTermCodon = length(grep("maintained_premature-termination", strsplit(metUtrCol," ")[[1]]))
    csv$LOST_METS_IN_5_UTR[i] = lostReadingFrame + maintainedReadingFrameWithTermCodon
  }
}
myvars=c('NMETS_5_UTR','METS_IN_5_UTR', 'CONSERVED_METS_IN_5_UTR', 'LOST_METS_IN_5_UTR')
View(csv[myvars])

# Ordenar variables
myvars=c('CHROMOSOME','GENE_ID', 'GENE_NAME', 'TRANSCRIPT_ID', 'TRANSCRIPT_REFSEQ_ID', 'TRANSCRIPT_BIOTYPE', 'METS_IN_5_UTR', 'NMETS_5_UTR', 'CONSERVED_METS_IN_5_UTR', 'LOST_METS_IN_5_UTR', 'SIGNAL_PEPTIDE_START', 'SIGNAL_PEPTIDE_END', 'CDS_ERRORS', 'PROTEIN_ID', 'VARIATION_NAME', 'VARIATION_TYPE', 'SOURCE', 'TRANSCRIPT_VARIATION_ALLELE_DBID', 'MINOR_ALLELE_FREQUENCY', 'CODON_CHANGE', 'CDS_COORDS', 'AMINOACID_CHANGE' ,'APPROACH1_MET_POSITION' ,'APPROACH1_STOP_CODON_POSITION' ,'APPROACH1_MUTATED_SEQUENCE_LENGTH' ,'APPROACH1_READING_FRAME_STATUS', 'APPROACH1_SIGNAL_PEPTIDE_CONSERVATION', 'APPROACH2_MET_POSITION', 'APPROACH2_INIT_CODON', 'APPROACH2_STOP_CODON_POSITION', 'APPROACH2_MUTATED_SEQUENCE_LENGTH', 'APPROACH2_SCORE', 'APPROACH2_READING_FRAME_STATUS', 'APPROACH2_SIGNAL_PEPTIDE_CONSERVATION', 'APPROACH3_MET_POSITION', 'APPROACH3_INIT_CODON', 'APPROACH3_STOP_CODON_POSITION', 'APPROACH3_MUTATED_SEQUENCE_LENGTH', 'APPROACH3_SCORE', 'APPROACH3_READING_FRAME_STATUS', 'APPROACH3_SIGNAL_PEPTIDE_CONSERVATION', 'CONSEQUENCE', 'PHENOTYPE', 'SO_TERM', 'SIFT', 'POLYPHEN', 'PUBLICATIONS')
csv = csv[myvars]
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

myvars=c('MINOR_ALLELE_FREQUENCY','CLASS')
View(csv[myvars])
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
rf_met_conserved = nrow(csv[csv$APPROACH1_READING_FRAME_STATUS == "Conserved",]); rf_met_conserved
rf_met_lost = nrow(csv[csv$APPROACH1_READING_FRAME_STATUS == "Lost",]); rf_met_lost
rf_kozak_conserved = nrow(csv[csv$APPROACH3_READING_FRAME_STATUS == "Conserved",]); rf_kozak_conserved
rf_kozak_lost = nrow(csv[csv$APPROACH3_READING_FRAME_STATUS == "Lost",]); rf_kozak_lost
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
summary(lowMaf$APPROACH1_MET_POSITION)
summary(lowMaf$APPROACH1_SIGNAL_PEPTIDE_CONSERVATION)
summary(lowMaf$APPROACH1_READING_FRAME_STATUS)


summary(highMaf$APPROACH1_MET_POSITION)
summary(highMaf$APPROACH1_SIGNAL_PEPTIDE_CONSERVARION)
summary(highMaf$APPROACH1_READING_FRAME_STATUS)

summary(lowMaf$APPROACH2_MET_POSITION)
summary(lowMaf$APPROACH2_SIGNAL_PEPTIDE_CONSERVARION)
summary(lowMaf$APPROACH2_READING_FRAME_STATUS)

summary(highMaf$APPROACH2_MET_POSITION)
summary(highMaf$APPROACH2_SIGNAL_PEPTIDE_CONSERVARION)
summary(highMaf$APPROACH2_READING_FRAME_STATUS)

summary(highMaf$APPROACH3_MET_POSITION)
summary(highMaf$APPROACH3_SIGNAL_PEPTIDE_CONSERVARION)
summary(highMaf$APPROACH3_READING_FRAME_STATUS)

summary(lowMaf$APPROACH3_MET_POSITION)
summary(lowMaf$APPROACH3_SIGNAL_PEPTIDE_CONSERVARION)
summary(lowMaf$APPROACH3_READING_FRAME_STATUS)

# Comprobar homogeneidad de varianzas
var.test(highMaf$APPROACH1_MET_POSITION, lowMaf$APPROACH1_MET_POSITION) # Varianzas distintas
var.test(highMaf$APPROACH2_MET_POSITION, lowMaf$APPROACH2_MET_POSITION) # Varianzas distintas
var.test(highMaf$APPROACH3_MET_POSITION, lowMaf$APPROACH3_MET_POSITION) # Varianzas distintas

#Test de normalidad
shapiro.test(highMaf$APPROACH1_MET_POSITION) # No normal
shapiro.test(lowMaf$APPROACH1_MET_POSITION) # No normal
shapiro.test(highMaf$APPROACH2_MET_POSITION) # No normal
shapiro.test(lowMaf$APPROACH2_MET_POSITION) # No normal
shapiro.test(highMaf$APPROACH3_MET_POSITION) # No normal
shapiro.test(lowMaf$APPROACH3_MET_POSITION) # No normal

# Test de wilcoxon para comparar medias
wilcox.test(highMaf$APPROACH1_MET_POSITION, lowMaf$APPROACH1_MET_POSITION, paired = F, conf.level = 0.95) # Distribuciones diferentes
wilcox.test(highMaf$APPROACH2_MET_POSITION, lowMaf$APPROACH2_MET_POSITION, paired = F, conf.level = 0.95) # Distribuciones diferentes
wilcox.test(highMaf$APPROACH3_MET_POSITION, lowMaf$APPROACH3_MET_POSITION, paired = F, conf.level = 0.95) # Distribuciones diferentes

# Boxplots met position
op <- par(mfrow = c(1, 3))
highMAFMedian = summary(highMaf$APPROACH1_MET_POSITION)[['Median']]
lowMAFMedian = summary(lowMaf$APPROACH1_MET_POSITION)[['Median']]
pValue = wilcox.test(highMaf$APPROACH1_MET_POSITION, lowMaf$APPROACH1_MET_POSITION, paired = F, conf.level = 0.95)$p.value
pValue = formatC(pValue, format = "e", digits = 2)
boxplot(highMaf$APPROACH1_MET_POSITION, lowMaf$APPROACH1_MET_POSITION, ylab="AUG position (bp)",names=c("High MAF", "Low MAF"), main="A")
legend("topleft", c(paste("High MAF median:",highMAFMedian), paste("Low MAF median:",lowMAFMedian), paste("P =",pValue)), cex=0.85, bty="n")

highMAFMedian = summary(highMaf$APPROACH2_MET_POSITION)[['Median']]
lowMAFMedian = summary(lowMaf$APPROACH2_MET_POSITION)[['Median']]
pValue = wilcox.test(highMaf$APPROACH2_MET_POSITION, lowMaf$APPROACH2_MET_POSITION, paired = F, conf.level = 0.95)$p.value
pValue = formatC(pValue, format = "e", digits = 2)
boxplot(highMaf$APPROACH2_MET_POSITION, lowMaf$APPROACH2_MET_POSITION, ylab="AUG position (bp)", names = c("High MAF", "Low MAF"), main="B", ylim=c(0, 1350))
legend("topleft", c(paste("High MAF median:",highMAFMedian), paste("Low MAF median:",lowMAFMedian), paste("P =",pValue)), cex=0.85, bty="n")

highMAFMedian = summary(highMaf$APPROACH3_MET_POSITION)[['Median']]
lowMAFMedian = summary(lowMaf$APPROACH3_MET_POSITION)[['Median']]
pValue = wilcox.test(highMaf$APPROACH3_MET_POSITION, lowMaf$APPROACH3_MET_POSITION, paired = F, conf.level = 0.95)$p.value
pValue = formatC(pValue, format = "e", digits = 2)
boxplot(highMaf$APPROACH3_MET_POSITION, lowMaf$APPROACH3_MET_POSITION, ylab="AUG position (bp)", names = c("High MAF", "Low MAF"), main="C")
legend("topleft", c(paste("High MAF median:",highMAFMedian), paste("Low MAF median:",lowMAFMedian), paste("P =",pValue)), cex=0.85, bty="n")
par(op)

# Test de Chi Cuadrado para comparar las variables cualitativas
# MET1 READING FRAME STATUS
rf_conserved_met1_low_maf = summary(lowMaf$APPROACH1_READING_FRAME_STATUS)[['Maintained']]
rf_lost_met1_low_maf = summary(lowMaf$APPROACH1_READING_FRAME_STATUS)[['Lost']]
rf_conserved_met1_high_maf = summary(highMaf$APPROACH1_READING_FRAME_STATUS)[['Maintained']]
rf_lost_met1_high_maf = summary(highMaf$APPROACH1_READING_FRAME_STATUS)[['Lost']]
m = as.table(rbind(c(rf_conserved_met1_low_maf,rf_lost_met1_low_maf), c(rf_conserved_met1_high_maf,rf_lost_met1_high_maf)))
dimnames(m)=list(GRUPO = c("MAF BAJA", "MAF ALTA"),
                 READING_FRAME_STATUS = c("MAINTAINED", "LOST"))
m
chisq.test(m) # p-value menor que 0.05 indica que cada grupo de mutaciones presenta diferencias significativas
###

# MET2 READING FRAME STATUS
rf_conserved_met2_low_maf = summary(lowMaf$APPROACH2_READING_FRAME_STATUS)[['Maintained']]
rf_lost_met2_low_maf = summary(lowMaf$APPROACH2_READING_FRAME_STATUS)[['Lost']]
rf_conserved_met2_high_maf = summary(highMaf$APPROACH2_READING_FRAME_STATUS)[['Maintained']]
rf_lost_met2_high_maf = summary(highMaf$APPROACH2_READING_FRAME_STATUS)[['Lost']]
m = as.table(rbind(c(rf_conserved_met2_low_maf,rf_lost_met2_low_maf), c(rf_conserved_met2_high_maf,rf_lost_met2_high_maf)))
dimnames(m)=list(GRUPO = c("MAF BAJA", "MAF ALTA"),
                 READING_FRAME_STATUS = c("MAINTAINED", "LOST"))
m
chisq.test(m) # p-value mayor que 0.05 indica que cada grupo de mutaciones no presenta diferencias significativas

# MET3 READING FRAME STATUS
rf_conserved_met3_high_maf = summary(lowMaf$APPROACH3_READING_FRAME_STATUS)[['Maintained']]
rf_lost_met3_high_maf = summary(lowMaf$APPROACH3_READING_FRAME_STATUS)[['Lost']]
rf_conserved_met3_low_maf = summary(highMaf$APPROACH3_READING_FRAME_STATUS)[['Maintained']]
rf_lost_met3_low_maf = summary(highMaf$APPROACH3_READING_FRAME_STATUS)[['Lost']]
m = as.table(rbind(c(rf_conserved_met3_low_maf,rf_lost_met3_low_maf), c(rf_conserved_met3_high_maf,rf_lost_met3_high_maf)))
dimnames(m)=list(GRUPO = c("MAF BAJA", "MAF ALTA"),
                 READING_FRAME_STATUS = c("MAINTAINED", "LOST"))
m
chisq.test(m) # p-value mayor que 0.05 indica que cada grupo de mutaciones no presenta diferencias significativas

# SIGNAL PEPTIDE AFFECTED MET1
summary(highMaf$APPROACH1_SIGNAL_PEPTIDE_CONSERVATION)
length(highMaf$APPROACH1_SIGNAL_PEPTIDE_CONSERVATION)
summary(lowMaf$APPROACH1_SIGNAL_PEPTIDE_CONSERVATION)
length(lowMaf$APPROACH1_SIGNAL_PEPTIDE_CONSERVATION)
# SIGNAL PEPTIDE AFFECTED MET2
summary(highMaf$APPROACH2_SIGNAL_PEPTIDE_CONSERVATION)
length(highMaf$APPROACH2_SIGNAL_PEPTIDE_CONSERVATION)
summary(lowMaf$APPROACH2_SIGNAL_PEPTIDE_CONSERVATION)
length(lowMaf$APPROACH2_SIGNAL_PEPTIDE_CONSERVATION)
# SIGNAL PEPTIDE AFFECTED MET3
summary(highMaf$APPROACH3_SIGNAL_PEPTIDE_CONSERVATION)
length(highMaf$APPROACH3_SIGNAL_PEPTIDE_CONSERVATION)
summary(lowMaf$APPROACH3_SIGNAL_PEPTIDE_CONSERVATION)
length(lowMaf$APPROACH3_SIGNAL_PEPTIDE_CONSERVATION)

# SIGNAL PEPTIDE COMPARISON
wilcox.test(highMaf$APPROACH1_SIGNAL_PEPTIDE_CONSERVATION, lowMaf$APPROACH1_SIGNAL_PEPTIDE_CONSERVATION, paired = F, conf.level = 0.95) # Distribuciones diferentes
wilcox.test(highMaf$APPROACH2_SIGNAL_PEPTIDE_CONSERVATION, lowMaf$APPROACH2_SIGNAL_PEPTIDE_CONSERVATION, paired = F, conf.level = 0.95) # Distribuciones diferentes
wilcox.test(highMaf$APPROACH3_SIGNAL_PEPTIDE_CONSERVATION, lowMaf$APPROACH3_SIGNAL_PEPTIDE_CONSERVATION, paired = F, conf.level = 0.95) # Distribuciones diferentes

# SIGNAL PEPTIDE BOXPLOTS
op <- par(mfrow = c(1, 3))
pValue = wilcox.test(highMaf$APPROACH1_SIGNAL_PEPTIDE_CONSERVATION, lowMaf$APPROACH1_SIGNAL_PEPTIDE_CONSERVATION, paired = F, conf.level = 0.95)$p.value
pValue = formatC(pValue, format = "e", digits = 2)
boxplot(highMaf$APPROACH1_SIGNAL_PEPTIDE_CONSERVATION, lowMaf$APPROACH1_SIGNAL_PEPTIDE_CONSERVATION, ylab="Conservation percentage of the signal peptide",names=c("High MAF", "Low MAF"), ylim=c(0,105), main="A")
legend("topleft", c(paste("P =",pValue)), cex=0.85, bty="n")

pValue = wilcox.test(highMaf$APPROACH2_SIGNAL_PEPTIDE_CONSERVATION, lowMaf$APPROACH2_SIGNAL_PEPTIDE_CONSERVATION, paired = F, conf.level = 0.95)$p.value
pValue = formatC(pValue, format = "e", digits = 2)
boxplot(highMaf$APPROACH2_SIGNAL_PEPTIDE_CONSERVATION, lowMaf$APPROACH2_SIGNAL_PEPTIDE_CONSERVATION, ylab="Conservation percentage of the signal peptide", names = c("High MAF", "Low MAF"), ylim=c(0,105), main="B")
legend("topleft", c(paste("P =",pValue)), cex=0.85, bty="n")

pValue = wilcox.test(highMaf$APPROACH3_SIGNAL_PEPTIDE_CONSERVATION, lowMaf$APPROACH3_SIGNAL_PEPTIDE_CONSERVATION, paired = F, conf.level = 0.95)$p.value
pValue = formatC(pValue, format = "e", digits = 2)
boxplot(highMaf$APPROACH3_SIGNAL_PEPTIDE_CONSERVATION, lowMaf$APPROACH3_SIGNAL_PEPTIDE_CONSERVATION, ylab="Conservation percentage of the signal peptide", names = c("High MAF", "Low MAF"), ylim=c(0,105), main="C")
legend("topleft", c(paste("P =",pValue)), cex=0.85, bty="n")
par(op)

# POSICIONES AFECTADAS
# Nos quedamos con las mutaciones que solo afectan a una posicion.
affectedPosHighMaf = c()
for (i in 1:(nrow(highMaf))){
  if (is.na(highMaf$CDS_COORDS[i]) || highMaf$CDS_COORDS[i] == ""){
    next
  }
  interval = strsplit(highMaf$CDS_COORDS[i], ",")
  
  begin = interval[[1]][[1]]
  begin = substr(begin, 2, nchar(begin))
  begin = as.numeric(begin)
  
  end = interval[[1]][[2]]
  end = substr(end, 2, nchar(end) - 1)
  end = as.numeric(end)
  if(begin == end){
    affectedPosHighMaf = rbind(affectedPosHighMaf, cbind(highMaf[i,], AFFECTED_POS = begin))
  }
}

affectedPosLowMaf = c()
for (i in 1:(nrow(lowMaf))){
  if (is.na(lowMaf$CDS_COORDS[i]) || lowMaf$CDS_COORDS[i] == ""){
    next
  }
  interval = strsplit(lowMaf$CDS_COORDS[i], ",")
  
  begin = interval[[1]][[1]]
  begin = substr(begin, 2, nchar(begin))
  begin = as.numeric(begin)
  
  end = interval[[1]][[2]]
  end = substr(end, 2, nchar(end) - 1)
  end = as.numeric(end)
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
numberOfMet5UTRLowMAF=lowMaf$NMETS_5_UTR
plot(density(numberOfMet5UTRLowMAF), col=gray(0.4), main="", xlab = "Number of \"AUG\" codons in 5' UTR")
numberOfMet5UTRHighMAF=highMaf$NMETS_5_UTR
lines(density(numberOfMet5UTRHighMAF), col =gray(0), lty=5)
legend(30, 0.68, c("Low MAF","High MAF"), lty=c(1,5), col = c(gray(0.4),gray(0)))
summary(highMaf$METS_IN_5_UTR)
summary(lowMaf$METS_IN_5_UTR)

wilcox.test(numberOfMet5UTRLowMAF, numberOfMet5UTRHighMAF, paired = F, conf.level = 0.95)


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

# Aciertos y errores de sift y polyphen
aciertosSift = 0
fallosSift = 0
aciertosPolyphen = 0
fallosPolyphen = 0
for (i in 1:nrow(csvWithMaf)) {
  class = csvWithMaf$CLASS[i]
  sift = csvWithMaf$SIFT[i]
  polyphen = csvWithMaf$POLYPHEN[i]
  
  if(class == 'BENIGN'){
    if(grepl('tolerated', sift, ignore.case = T)){
      aciertosSift = aciertosSift +1
    } else{
      fallosSift = fallosSift+1
    }
    
    if(grepl('benign', polyphen, ignore.case = T)){
      aciertosPolyphen = aciertosPolyphen + 1
    } else {
      fallosPolyphen = fallosPolyphen + 1
    }
  } else if (class == 'DELETERIOUS'){
    if(grepl('deleterious', sift, ignore.case = T)){
      aciertosSift = aciertosSift +1
    } else{
      fallosSift = fallosSift+1
    }
    
    if(grepl('damaging', polyphen, ignore.case = T)){
      aciertosPolyphen = aciertosPolyphen + 1
    } else {
      fallosPolyphen = fallosPolyphen + 1
    }
  }
}
porcentajeAciertoSift = aciertosSift/(aciertosSift+fallosSift)
porcentajeAciertoPolyphen = aciertosPolyphen/(aciertosPolyphen+fallosPolyphen)
myvars=c( "MINOR_ALLELE_FREQUENCY", "SIFT","POLYPHEN","CLASS")
View(highMaf[myvars])
