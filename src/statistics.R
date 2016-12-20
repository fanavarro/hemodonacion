setwd("~/hemodonacion/data/tsv")
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
        if (first_met >= (signal_end + 1)){signal_first_met_affected = "Lost"}
      }
      if (!is.na(first_kozak)){
        if (first_kozak == signal_start){signal_first_kozak_affected = "Totally conserved"}
        if (first_kozak > signal_start && first_kozak <= signal_end){signal_first_kozak_affected = "Partially conserved"}
        if (first_kozak >= (signal_end + 1)){signal_first_kozak_affected = "Lost"}
      }
    }
    csv[i, "SIGNAL_FIRST_MET_AFFECTED"] = signal_first_met_affected
    csv[i, "SIGNAL_FIRST_KOZAK_AFFECTED"] = signal_first_kozak_affected
  }
  return(csv)
}

add_mutation_type = function(csv){
  for (i in 1:nrow(csv)){
    codon_change = trimws(csv[i, "CODON_CHANGE"])
    if (!is.na(codon_change) && "" != codon_change){
      splitted = strsplit(codon_change, "/")
      original = splitted[[1]][1]
      mutated = splitted[[1]][2]
      if(nchar(original) > nchar(mutated)){
        csv[i, "VARIATION_TYPE"] = "Deletion"
      } else if(nchar(original) < nchar(mutated)){
        csv[i, "VARIATION_TYPE"] = "Insertion"
      } else {
        csv[i, "VARIATION_TYPE"] = "Nucleotid change"
      }
    }
  }
  return(csv)
}

# Leer el csv sin filtros
csv = read.csv("11_12_2016.csv", sep="\t",stringsAsFactors=FALSE)
csv[,"MUTATED_SEQUENCE_LENGTH_1"]=as.numeric(gsub("%","",csv$MUTATED_SEQUENCE_LENGTH_1))
csv[,"MUTATED_SEQUENCE_LENGTH_2"]=as.numeric(gsub("%","",csv$MUTATED_SEQUENCE_LENGTH_2))
csv[,"MUTATED_SEQUENCE_LENGTH_3"]=as.numeric(gsub("%","",csv$MUTATED_SEQUENCE_LENGTH_3))
csv[,"SIGNAL_PEPTIDE_CONSERVATION_1"]=as.numeric(gsub("%","",csv$SIGNAL_PEPTIDE_CONSERVATION_1))
csv[,"SIGNAL_PEPTIDE_CONSERVATION_2"]=as.numeric(gsub("%","",csv$SIGNAL_PEPTIDE_CONSERVATION_2))
csv[,"SIGNAL_PEPTIDE_CONSERVATION_3"]=as.numeric(gsub("%","",csv$SIGNAL_PEPTIDE_CONSERVATION_3))
#csv = add_signal_lost_sup_info(csv)
#csv = add_mutation_type(csv)

# Generar el fichero con las columnas suplementarias
#write.table(csv, file = "final_out_no_filter_sup.csv", na="", sep="\t", row.names = F)


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





# Eliminar los casos en los que hay errores en las regiones 5' o 3'
#csv = csv[csv$CDS_ERRORS == '',]

# Eliminar los casos en los que el biotipo es "non_stop_decay" o "nonsense_mediated_decay"
#csv = csv[csv$TRANSCRIPT_BIOTYPE != 'non_stop_decay' & csv$TRANSCRIPT_BIOTYPE != 'nonsense_mediated_decay',]

# Generar el fichero filtrado
#write.table(csv, file = "final_out_filter.csv", na="", sep="\t", row.names = F)

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

# Transcritos afectados de TP53, CACNA1C y TP53
length(unique(csv[csv$GENE_NAME=="TP53",]$TRANSCRIPT_ID))
length(unique(csv[csv$GENE_NAME=="CACNA1C",]$TRANSCRIPT_ID))
length(unique(csv[csv$GENE_NAME=="CDKN2A",]$TRANSCRIPT_ID))

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
summary(highMaf$SIGNAL_PEPTIDE_CONSERVARION_1)

summary(lowMaf$MET_POSITION_2)
summary(lowMaf$SIGNAL_PEPTIDE_CONSERVARION_2)
summary(lowMaf$READING_FRAME_STATUS_2)
rf_conserved_met2_low_maf = nrow(lowMaf[lowMaf$READING_FRAME_STATUS_2 == 'Conserved',])
rf_lost_met2_low_maf = nrow(lowMaf[lowMaf$READING_FRAME_STATUS_2 == 'Lost',])
summary(lowMaf$SIGNAL_PEPTIDE_CONSERVARION_2)

summary(highMaf$MET_POSITION_2)
summary(highMaf$SIGNAL_PEPTIDE_CONSERVARION_2)
summary(highMaf$READING_FRAME_STATUS_2)
rf_conserved_met2_high_maf = nrow(highMaf[highMaf$READING_FRAME_STATUS_2 == 'Conserved',])
rf_lost_met2_high_maf = nrow(highMaf[highMaf$READING_FRAME_STATUS_2 == 'Lost',])
summary(highMaf$SIGNAL_PEPTIDE_CONSERVARION_2)


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

op <- par(mfrow = c(1, 3))
boxplot(highMaf$MET_POSITION_1, lowMaf$MET_POSITION_1, ylim=c(0,250), ylab="Posición del codón inicial (en pares de bases)",names=c("MAF alta", "MAF baja"),
        main="Comparativa de la posición del primer codón de inicio\nencontrado entre los grupos de MAF alta y baja.")
boxplot(highMaf$MET_POSITION_2, lowMaf$MET_POSITION_2, ylim=c(0,600), ylab="Posición del codón inicial (en pares de bases)", names = c("MAF alta", "MAF baja"),
        main="Comparativa de la posición del primer codón de inicio\npredicho por ATGpr con score > 0.25\nentre los grupos de MAF alta y baja.")
boxplot(highMaf$MET_POSITION_3, lowMaf$MET_POSITION_3, ylim=c(0,600), ylab="Posición del codón inicial (en pares de bases)", names = c("MAF alta", "MAF baja"),
        main="Comparativa de la posición del primer codón de inicio\nencontrado en un contexto de Kozak fuerte\nentre los grupos de MAF alta y baja.")

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

# PEPTIDE SIGNAL AFFECTED MET1
summary(highMaf$SIGNAL_PEPTIDE_CONSERVATION_1)
length(highMaf$SIGNAL_PEPTIDE_CONSERVATION_1)
summary(lowMaf$SIGNAL_PEPTIDE_CONSERVATION_1)
length(lowMaf$SIGNAL_PEPTIDE_CONSERVATION_1)
# PEPTIDE SIGNAL AFFECTED MET2
summary(highMaf$SIGNAL_PEPTIDE_CONSERVATION_2)
length(highMaf$SIGNAL_PEPTIDE_CONSERVATION_2)
summary(lowMaf$SIGNAL_PEPTIDE_CONSERVATION_2)
length(lowMaf$SIGNAL_PEPTIDE_CONSERVATION_2)
# PEPTIDE SIGNAL AFFECTED MET3
summary(highMaf$SIGNAL_PEPTIDE_CONSERVATION_3)
length(highMaf$SIGNAL_PEPTIDE_CONSERVATION_3)
summary(lowMaf$SIGNAL_PEPTIDE_CONSERVATION_3)
length(lowMaf$SIGNAL_PEPTIDE_CONSERVATION_3)

# PEPTIDE SIGNAL BOXPLOTS
op <- par(mfrow = c(1, 3))
boxplot(highMaf$SIGNAL_PEPTIDE_CONSERVATION_1, lowMaf$SIGNAL_PEPTIDE_CONSERVATION_1, ylab="Conservacion del peptido señal (en %)",names=c("MAF alta", "MAF baja"),
        main="Comparativa de la conservacion del peptido senal\nencontrado entre los grupos de MAF alta y baja.")
boxplot(highMaf$SIGNAL_PEPTIDE_CONSERVATION_2, lowMaf$SIGNAL_PEPTIDE_CONSERVATION_2, ylab="Conservacion del peptido señal (en %)", names = c("MAF alta", "MAF baja"),
        main="Comparativa de la conservacion del peptido señal\ncon la met de ATGpr con score > 0.25\nentre los grupos de MAF alta y baja.")
boxplot(highMaf$SIGNAL_PEPTIDE_CONSERVATION_3, lowMaf$SIGNAL_PEPTIDE_CONSERVATION_3, ylab="Conservacion del peptido señal (en %)", names = c("MAF alta", "MAF baja"),
        main="Comparativa de la conservacion del peptido señal\nusando la MET en contexto de Kozak fuerte\nentre los grupos de MAF alta y baja.")
par(op)


# select variables v1, v2, v3
myvars <- c("TRANSCRIPT_ID","FIRST_MET_POSITION","STOP_CODON_POSITION","MUTATED_SEQUENCE_LENGTH", "KOZAK_START", "KOZAK_END", "KOZAK_MUTATED_SEQUENCE_LENGTH", "KOZAK_STOP_CODON", "KOZAK_READING_FRAME_STATUS")
myvars=c("GENE_NAME","VARIATION_NAME","MINOR_ALLELE_FREQUENCY", "FIRST_MET_POSITION","READING_FRAME_STATUS", "KOZAK_START", "KOZAK_READING_FRAME_STATUS", "CODON_CHANGE")
View(highMaf[myvars])
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


csv = add_kozak_mutated_seq_length(csv)
myvars=c("FIRST_MET_POSITION","STOP_CODON_POSIION","MUTATED_SEQUENCE_LENGTH", "KOZAK_START", "KOZAK_END", "KOZAK_MUTATED_SEQUENCE_LENGTH")
View(csv[myvars])
