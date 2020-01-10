setwd("~/hemodonacion/ML")
library("RWeka")
library("partykit")
library(caret)
set.seed(12345)
# Leer el csv sin filtros
# csv = read.csv("out.csv", sep="\t",stringsAsFactors=FALSE)
csv = read.csv("out-no-filtro-transcripts.csv", sep="\t",stringsAsFactors=FALSE)


# Nos quedamos solo con los SNPs que afectan a una sola posicion
# csv = csv[csv$VARIATION_TYPE=='SNP',]
rowsToDelete = c()
for(i in 1:nrow(csv)){
  if(csv$CDS_COORDS[i] == ''){
    csv$CDS_COORDS[i] = NA
  } else {
    str = csv$CDS_COORDS[i]
    range = as.numeric(trimws(strsplit(substr(str, 2, nchar(str)-1),',')[[1]]))
    if(range[1] != range[2]){
      rowsToDelete = c(rowsToDelete, i)
    }
  }
}
csv = csv[-c(rowsToDelete), ]
row.names(csv) = NULL
nrow(csv[csv$CLASS=='BENIGN',])
nrow(csv[csv$CLASS=='DELETERIOUS',])
csv[,"MUTATED_SEQUENCE_LENGTH"]=as.numeric(gsub("%","",csv$MUTATED_SEQUENCE_LENGTH))

for(i in 1:nrow(csv)){
  if(csv$READING_FRAME_STATUS[i] == ''){
    csv$READING_FRAME_STATUS[i] = NA
  }
}


# Calcular metioninas en 5' utr
csv$NMETS_5_UTR = sapply(strsplit(csv$METS_IN_5_UTR, " "), length)

# Contar los que mantienen fase de lectura
for (i in 1:nrow(csv)) {
  metUtrCol = csv$METS_IN_5_UTR[i]
  if(is.na(metUtrCol) || metUtrCol == ""){
    csv$CONSERVED_METS_IN_5_UTR[i] = 0
  } else {
    csv$CONSERVED_METS_IN_5_UTR[i] = length(grep("maintained", strsplit(metUtrCol," ")[[1]]))
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
#Contar las que mantienen la fase de lectura sin encontrar un codon stop prematuro
for (i in 1:nrow(csv)) {
  metUtrCol = csv$METS_IN_5_UTR[i]
  if(is.na(metUtrCol) || metUtrCol == ""){
    csv$CONSERVED_METS_NO_STOP_IN_5_UTR[i] = 0
  } else {
    csv$CONSERVED_METS_NO_STOP_IN_5_UTR[i] = length(grep("maintained_no-premature-termination", strsplit(metUtrCol," ")[[1]]))
  }
}

# Ordenar y eliminar variables sin importancia

myvars=c('NMETS_5_UTR', 'CONSERVED_METS_IN_5_UTR', 'LOST_METS_IN_5_UTR', 'CONSERVED_METS_NO_STOP_IN_5_UTR', 'CDS_COORDS', 'AMINOACID_CHANGE','CODON_CHANGE','MET_POSITION', 'READING_FRAME_STATUS','NO_STOP_CODON','PREMATURE_STOP_CODON', 'STOP_CODON_POSITION', 'MUTATED_SEQUENCE_LENGTH','CLASS')
csv = csv[myvars]
csv$READING_FRAME_STATUS = factor(csv$READING_FRAME_STATUS)
csv$CODON_CHANGE = factor(csv$CODON_CHANGE)
csv$AMINOACID_CHANGE = factor(csv$AMINOACID_CHANGE)
csv$PREMATURE_STOP_CODON = factor(csv$PREMATURE_STOP_CODON)
csv$NO_STOP_CODON = factor(csv$NO_STOP_CODON)
csv$CDS_COORDS = factor(csv$CDS_COORDS)
csv$CLASS = factor(csv$CLASS)

summary(csv)
# Eliminar NA
csv = na.exclude(csv)
write.table(csv, file = "/home/fabad/mutaciones.tsv", sep="\t", row.names = F)

underSample = downSample(csv[, 1:length(csv)-1], csv[,length(csv)], yname='CLASS')
test = setdiff(csv,underSample)
ctrl <- trainControl(method = "repeatedcv", 
                     number = 10, 
                     repeats = 10, 
                     verboseIter = FALSE)

model <- caret::train(underSample[, 1:length(csv)-1],
                              underSample[,length(csv)],
                              method = "J48",
                              preProcess = c("scale", "center"),
                              trControl = ctrl)


plot(model$finalModel)

# Exportar a weka
write.arff(underSample, "out-no-filtro-transcripts.arff")

View(csv[c('CONSERVED_METS_NO_STOP_IN_5_UTR')])
