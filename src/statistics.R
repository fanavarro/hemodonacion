setwd("/home/fabad/hemodonacion/src")
common=read.table('hist_common.dat')
no_common=read.table('hist_no_common.dat')
plot(common, type='h', xlim=c(0,1000), xlab="Posicion de la primera metionina en la secuencia mutada", ylab="Numero de casos", main="MAF mayor al 1%")
plot(no_common, type='h', xlab="Posicion de la primera metionina en la secuencia mutada", ylab="Numero de casos", main="MAF menor al 1%")


csv = read.csv("by_pos_manual_filter.csv", sep="\t",stringsAsFactors=FALSE)
csvWithMaf = csv[csv$MINOR_ALLELE_FREQUENCY != "-",]

# Eliminar el * en los casos donde se ha obtenido la secuencia
# mutada sin tener informacion del alelo, eliminando directamente
# la primera metionina.
csvWithMaf$FIRST_MET_POSITION = gsub("\\*$", "", csvWithMaf$FIRST_MET_POSITION, perl=TRUE)

highMaf = csvWithMaf[csvWithMaf$MINOR_ALLELE_FREQUENCY >= 0.01,]
lowMaf = csvWithMaf[csvWithMaf$MINOR_ALLELE_FREQUENCY < 0.01,]
noMaf = csv[csv$MINOR_ALLELE_FREQUENCY == "-",]


highMaf$FIRST_MET_POSITION = as.numeric(highMaf$FIRST_MET_POSITION)
lowMaf$FIRST_MET_POSITION = as.numeric(lowMaf$FIRST_MET_POSITION)
mean(highMaf$FIRST_MET_POSITION, na.rm=TRUE)
mean(lowMaf$FIRST_MET_POSITION, na.rm=TRUE)
summary(lowMaf)
summary(highMaf)

hist(highMaf$FIRST_MET_POSITION)
hist(lowMaf$FIRST_MET_POSITION)
highMaf$FIRST_MET_POSITION[1] + highMaf$FIRST_MET_POSITION[2]

# Comprobar homogeneidad de varianzas
var.test(highMaf$FIRST_MET_POSITION, lowMaf$FIRST_MET_POSITION)

#Test de normalidad
shapiro.test(highMaf$FIRST_MET_POSITION)
shapiro.test(lowMaf$FIRST_MET_POSITION)
# Test de wilcoxon para comparar medias
wilcox.test(highMaf$FIRST_MET_POSITION, lowMaf$FIRST_MET_POSITION, paired = F, conf.level = 0.95)

boxplot(highMaf$FIRST_MET_POSITION, lowMaf$FIRST_MET_POSITION)
