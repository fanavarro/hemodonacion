setwd("/home/fabad/hemodonacion/src")

# Leer el csv sin filtros
csv = read.csv("kozak.csv", sep="\t",stringsAsFactors=FALSE)
csv$READING_FRAME_STATUS = factor(csv$READING_FRAME_STATUS)
csv$KOZAK_READING_FRAME_STATUS = factor(csv$KOZAK_READING_FRAME_STATUS)
csv$KOZAK_STOP_CODON = factor(csv$KOZAK_STOP_CODON)
csv$STOP_CODON_POSITION = factor(csv$STOP_CODON_POSITION)
csv$GENE_NAME = factor(csv$GENE_NAME)

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
summary(lowMaf)
summary(highMaf)

# Histogramas de la posición de la primera metionina
hist(highMaf$FIRST_MET_POSITION)
hist(lowMaf$FIRST_MET_POSITION)

# Histogramas de la posición de la metionina de la primera secuencia Kozak
hist(highMaf$KOZAK_START)
hist(lowMaf$KOZAK_START)

# Comprobar homogeneidad de varianzas
var.test(highMaf$FIRST_MET_POSITION, lowMaf$FIRST_MET_POSITION) # Varianzas distintas
var.test(highMaf$KOZAK_START, lowMaf$KOZAK_START) # Varianzas iguales

#Test de normalidad
shapiro.test(highMaf$FIRST_MET_POSITION) # No normal
shapiro.test(lowMaf$FIRST_MET_POSITION) # No normal
shapiro.test(highMaf$KOZAK_START) # No normal
shapiro.test(lowMaf$KOZAK_START) # No normal

# Test de wilcoxon para comparar medias
wilcox.test(highMaf$FIRST_MET_POSITION, lowMaf$FIRST_MET_POSITION, paired = F, conf.level = 0.95) # Distribuciones diferentes
wilcox.test(highMaf$KOZAK_START, lowMaf$KOZAK_START, paired = F, conf.level = 0.95) # Distribuciones iguales

# Boxplots
boxplot(highMaf$FIRST_MET_POSITION, lowMaf$FIRST_MET_POSITION)
boxplot(highMaf$KOZAK_START, lowMaf$KOZAK_START)

