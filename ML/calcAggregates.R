setwd("~/hemodonacion/ML")
csv = read.csv("SVMEvaluationResults.tsv", sep="\t",stringsAsFactors=FALSE)

agregados = aggregate(csv[,4:ncol(csv)], by=list(csv$Classifier, csv$Gamma, csv$Cost), FUN = 'mean')
colnames(agregados) = colnames(csv)
agregados = agregados[,-4]
write.table(agregados, file="AgregadosSVMEvaluationResults.tsv", sep="\t", row.names = FALSE)
