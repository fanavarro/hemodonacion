setwd("~/hemodonacion/data/tsv")
#anotaciones
library(GO.db)
library(AnnotationHub)
library(org.Hs.eg.db)

#Enriquecimiento
library(clusterProfiler)
library(topGO)

#graficos
library(ggplot2)
library(plotly)
library(Rgraphviz)

# Funciones
# Dibuja un plot con el resultado de groupGO
plotGGO = function(groupGOResult, countThreshold = 0){
  myggo<-groupGOResult@result[with(groupGOResult@result,	Count > countThreshold),	]
  #myggo = myggo[with(myggo, order("Count", decreasing = TRUE), ),]
  
  p<-ggplot(data=myggo, aes(reorder(Description, Count), Count))+
    geom_bar(stat = "identity", na.rm = TRUE)+
    coord_flip()+theme(axis.title.y=element_blank() )
  p+scale_fill_hue(h=c(10,300), c=100, l=70, direction=1)
}

#Leer csv
csv = read.csv("9_12_2017.csv", sep="\t",stringsAsFactors=FALSE)

# Dividir segun maf
csvWithMaf = csv[!is.na(csv$MINOR_ALLELE_FREQUENCY),]
highMaf = csvWithMaf[csvWithMaf$MINOR_ALLELE_FREQUENCY >= 0.01,]
lowMaf = csvWithMaf[csvWithMaf$MINOR_ALLELE_FREQUENCY < 0.01,]

# Obtener identificadores de genes de cada grupo
highMafEnsemblGenes = unique(as.vector(highMaf$GENE_ID))
lowMafEnsemblGenes = unique(as.vector(lowMaf$GENE_ID))

# Agrupar por terminos de Gene Ontology MF (Molecular Function)
GOLevel = 2
highMafGroupGO = groupGO(gene = highMafEnsemblGenes, OrgDb = org.Hs.eg.db, ont = "MF", level = GOLevel, readable = TRUE, keytype = "ENSEMBL")
lowMafGroupGO = groupGO(gene = lowMafEnsemblGenes, OrgDb = org.Hs.eg.db, ont = "MF", level = GOLevel, readable = TRUE, keytype = "ENSEMBL")

# Visualizacion
plotGGO(highMafGroupGO)
plotGGO(lowMafGroupGO)

# Agrupar por terminos de Gene Ontology BP (Biological Process)
highMafGroupGO = groupGO(gene = highMafEnsemblGenes, OrgDb = org.Hs.eg.db, ont = "BP", level = GOLevel, readable = TRUE, keytype = "ENSEMBL")
lowMafGroupGO = groupGO(gene = lowMafEnsemblGenes, OrgDb = org.Hs.eg.db, ont = "BP", level = GOLevel, readable = TRUE, keytype = "ENSEMBL")

# Visualizacion
# countThreshold = 3
plotGGO(highMafGroupGO)
# countThreshold = 50
plotGGO(lowMafGroupGO)


# Enriquecimiento molecular function
pvalueCutoff = 0.05
qvalueCutoff = 0.05
highMafEnrichment = enrichGO(gene = highMafEnsemblGenes ,OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "MF", pAdjustMethod = "BH", pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff)
highMafEnrichment = setReadable(highMafEnrichment, OrgDb = org.Hs.eg.db)
lowMafEnrichment = enrichGO(gene = lowMafEnsemblGenes ,OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "MF", pAdjustMethod = "BH", pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff)
lowMafEnrichment = setReadable(lowMafEnrichment, OrgDb = org.Hs.eg.db)

# Visualizacion
dotplot(highMafEnrichment)
dotplot(lowMafEnrichment)

# Enriquecimiento biological process
pvalueCutoff = 0.05
qvalueCutoff = 0.05
highMafEnrichment = enrichGO(gene = highMafEnsemblGenes ,OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff)
highMafEnrichment = setReadable(highMafEnrichment, OrgDb = org.Hs.eg.db)
lowMafEnrichment = enrichGO(gene = lowMafEnsemblGenes ,OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff)
lowMafEnrichment = setReadable(lowMafEnrichment, OrgDb = org.Hs.eg.db)

# Visualizacion
dotplot(highMafEnrichment)
dotplot(lowMafEnrichment)
