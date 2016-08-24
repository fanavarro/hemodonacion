library(gdata)
library(hash)
library(descr)

search = function(database, xls){
  #found = data.frame(matrix(ncol=ncol(database)))
  #colnames(found) = colnames(database)
  found = read.table(text="", col.names = colnames(database))
  # Si el xlsx tiene entrada para dbsnp e id de refseq
  if("dbsnp" %in% colnames(xls) && "transcript" %in% colnames(xls)){
    for (i in 1:nrow(xls)){
      dbsnp = toString(xls[i, "dbsnp"])
      refseq_id = toString(xls[i,"transcript"])
      entry = database[database$TRANSCRIPT_REFSEQ_ID == refseq_id & !is.na(database$TRANSCRIPT_REFSEQ_ID) & database$VARIATION_NAME == dbsnp & !is.na(database$VARIATION_NAME), ]
      if (!is.na(entry) && nrow(entry) != 0){
        found = rbind(found, entry)
      }
      
    }
  }
  # Si tiene id de transcrito de ensembl y cambio de nucleotido del tipo aTg/aAg
  else if("Feature" %in% colnames(xls) && "Codons" %in% colnames(xls)){
    for (i in 1:nrow(xls)){
      transcript_id = xls[i, "Feature"]
      codon_change = xls[i,"Codons"]
      entry = database[database$TRANSCRIPT_ID == transcript_id & !is.na(database$TRANSCRIPT_ID) & database$CODON_CHANGE == codon_change & !is.na(database$CODON_CHANGE), ]
      if(!is.na(entry) && nrow(entry) != 0){
        found = rbind(found, entry)
      }
    }
  }
  return(found)
}

View(search(csv, exome_tables[["14-173.xlsx"]]))
View(search(csv, exome_tables[["2064.xlsx"]]))
View(search(csv, exome_tables[["ABGP.xlsx"]]))
View(search(csv, exome_tables[["Exoma 10-164.xlsx"]]))
View(search(csv, exome_tables[["Exoma 11-584.xlsx"]]))
View(search(csv, exome_tables[["Exoma 2166 nuevo.xlsx"]]))
View(search(csv, exome_tables[["Paqui.xlsx"]]))

filter_xls = function(xls){
  minCoverage = 20
  minQuality = 60
  
  # Filter the position affected by the mutation
  if("proteinPos" %in% colnames(xls)){
    xls = xls[xls$proteinPos == 1,]
  } else if ("Protein_position" %in% colnames(xls)){
    xls = xls[xls$Protein_position == 1,]
  } else if ("protein" %in% colnames(xls)){
    index = with(xls, grepl("p.Met1[a-zA-Z]", protein))
    xls = xls[index,]
  }
  
  # Filter by coverage and quality
  if ("coverage" %in% colnames(xls)){
    xls = xls[xls$coverage >= minCoverage,]
  } else if ("QUAL" %in% colnames(xls)){
    xls = xls[xls$QUAL >= minQuality,]
  }
  return(xls)
}

get_gene_names = function(xls){
  if("gene" %in% colnames(xls)){
    return(sort(xls$gene))
  } else if("GeneSymbol" %in% colnames(xls)){
    return(sort(xls$GeneSymbol))
  } else if("Gene.Name" %in% colnames(xls)){
    return(sort(xls$Gene.Name))
  }
}

exomes_dir = "/home/fabad/hemodonacion/data/exomes/"
exome_files = c("14-173.xlsx",
                "2064.xlsx",
                "ABGP.xlsx",
                "Exoma 10-164.xlsx",
                "Exoma 11-584.xlsx",
                "Exoma 2166 nuevo.xlsx",
                "Paqui.xlsx")

exome_tables = hash()

# Read xlsx files and store it in a hash
for (file in exome_files){
  abs_file_dir = paste(exomes_dir, file, sep="")
  exome_tables[file] = read.xls(abs_file_dir)
  write(paste("exoma", file, "leido."), "")
}

# filter data
for (file in exome_files){
  exome_tables[file] = filter_xls(exome_tables[[file]])
}

# get genes
genes = c()
for (file in exome_files){
  genes = c(genes, as.vector(get_gene_names(exome_tables[[file]])))
}

# gene frequency
op = par(las=2) # make label text perpendicular to axis
op = par(mar=c(5,8,4,2)) # increase y-axis margin.
barplot(table(genes), horiz=T)
par(op)
  View(exome_tables[["14-173.xlsx"]])
  View(exome_tables[["2064.xlsx"]])
  View(exome_tables[["ABGP.xlsx"]])
  View(exome_tables[["Exoma 10-164.xlsx"]])
  View(exome_tables[["Exoma 11-584.xlsx"]])
  View(exome_tables[["Exoma 2166 nuevo.xlsx"]])
  View(exome_tables[["Paqui.xlsx"]])
  

  get_gene_names(exome_tables[["14-173.xlsx"]])
  get_gene_names(exome_tables[["2064.xlsx"]])
  get_gene_names(exome_tables[["ABGP.xlsx"]])
  get_gene_names(exome_tables[["Exoma 10-164.xlsx"]])
  get_gene_names(exome_tables[["Exoma 11-584.xlsx"]])
  get_gene_names(exome_tables[["Exoma 2166 nuevo.xlsx"]])
  get_gene_names(exome_tables[["Paqui.xlsx"]])

# Numero de mutaciones encontradas individualmente
nrow(exome_tables[["14-173.xlsx"]])
nrow(exome_tables[["2064.xlsx"]])
nrow(exome_tables[["ABGP.xlsx"]])
nrow(exome_tables[["Exoma 10-164.xlsx"]])
nrow(exome_tables[["Exoma 11-584.xlsx"]])
nrow(exome_tables[["Exoma 2166 nuevo.xlsx"]])
nrow(exome_tables[["Paqui.xlsx"]])