library(gdata)
library(hash)

search = function(database, xls){
  found = read.table(text="", col.names = colnames(database))
  # Si el xlsx tiene entrada para dbsnp e id de refseq
  if("dbsnp" %in% colnames(xls) && "transcript" %in% colnames(xls)){
    for (i in 1:nrow(xls)){
      dbsnpString = toString(xls[i, "dbsnp"])
      dbsnps = strsplit(dbsnpString, ":")
      refseq_id = toString(xls[i,"transcript"])
      for(i in 1:length(dbsnps[[1]])){
        dbsnp = dbsnps[[1]][i]
        entry = database[database$TRANSCRIPT_REFSEQ_ID == refseq_id & 
                           !is.na(database$TRANSCRIPT_REFSEQ_ID) & 
                           database$VARIATION_NAME == dbsnp & 
                           !is.na(database$VARIATION_NAME), ]
        
        if (!is.na(entry) && nrow(entry) != 0){
          found = rbind(found, entry)
          break
        }
      }
    }
  }
  return(found)
}

filter_xls = function(xls){
  xls = filter_coverage_xls(xls)
  xls = filter_protein_pos_xls(xls)
  return(xls)
}
filter_coverage_xls = function(xls){
  minCoverage = 20
  
  # Filter by coverage and quality
  if ("coverage" %in% colnames(xls)){
    xls = xls[!is.na(xls$coverage) & xls$coverage >= minCoverage,]
  }
  return(xls)
}

filter_protein_pos_xls = function(xls){
  # Filter the position affected by the mutation
  if("proteinPos" %in% colnames(xls)){
    xls = xls[!is.na(xls$proteinPos) & xls$proteinPos == 1,]
  } else if ("protein" %in% colnames(xls)){
    index = with(xls, grepl("p.Met1[a-zA-Z]", protein))
    xls = xls[index,]
  }
  return(xls)
}

get_gene_names = function(xls){
  if("gene" %in% colnames(xls)){
    return(sort(paste(xls$gene," (", xls$protein,") ","MAF = ", xls$maf, sep = "")))
  }
}



get_dbsnp_ids = function(xls){
  if("dbsnp" %in% colnames(xls)){
    return(sort(xls[!is.na(xls$dbsnp), "dbsnp"]))
  }
}

exomes_dir = "~/hemodonacion/data/exomes/"
exome_files = c("14-173.xlsx",
                "2064.xlsx",
                "ABGP.xlsx",
                "Exoma 2166 nuevo.xlsx",
                "Paqui.xlsx")

exome_tables = hash()

# Read xlsx files and store it in a hash
for (file in exome_files){
  abs_file_dir = paste(exomes_dir, file, sep="")
  exome_tables[file] = read.xls(abs_file_dir)
  write(paste("exoma", file, "leido."), "")
}



# filter coverage 
for (file in exome_files){
  exome_tables[file] = filter_coverage_xls(exome_tables[[file]])
}

# Mutaciones totales por paciente
nrow(exome_tables[["14-173.xlsx"]])
nrow(exome_tables[["2064.xlsx"]])
nrow(exome_tables[["ABGP.xlsx"]])
nrow(exome_tables[["Exoma 2166 nuevo.xlsx"]])
nrow(exome_tables[["Paqui.xlsx"]])
mean(c(nrow(exome_tables[["14-173.xlsx"]]),nrow(exome_tables[["2064.xlsx"]]),nrow(exome_tables[["ABGP.xlsx"]]),nrow(exome_tables[["Exoma 2166 nuevo.xlsx"]]),nrow(exome_tables[["Paqui.xlsx"]])))
# filter protein pos
for (file in exome_files){
  exome_tables[file] = filter_protein_pos_xls(exome_tables[[file]])
}

# Crear csv con mtuaciones en met1
for (file in exome_files){
  write.table(exome_tables[[file]], file = paste(exomes_dir, file, ".tsv", sep = ""), sep = "\t", row.names = F)
}

# Mutaciones en met1 por paciente
nrow(exome_tables[["14-173.xlsx"]])
nrow(exome_tables[["2064.xlsx"]])
nrow(exome_tables[["ABGP.xlsx"]])
nrow(exome_tables[["Exoma 2166 nuevo.xlsx"]])
nrow(exome_tables[["Paqui.xlsx"]])
median(c(nrow(exome_tables[["14-173.xlsx"]]),nrow(exome_tables[["2064.xlsx"]]),nrow(exome_tables[["ABGP.xlsx"]]),nrow(exome_tables[["Exoma 2166 nuevo.xlsx"]]),nrow(exome_tables[["Paqui.xlsx"]])))

# get genes
genes = c()
for (file in exome_files){
  genes = c(genes, as.vector(get_gene_names(exome_tables[[file]])))
}
# manual normalization
genes[41]="LENG9 (p.Met1Lys) MAF = 0.448" 
genes[53]="LENG9 (p.Met1Lys) MAF = 0.448" 
genes[40] ="KRTAP4-8 (p.Met1Asn) MAF = 0.419" 
genes[52] ="KRTAP4-8 (p.Met1Asn) MAF = 0.419" 
genes[10]="STK11IP (p.Arg1Met) MAF = 0.142"
genes[22]="STK11IP (p.Arg1Met) MAF = 0.142"
genes[42]="OR8K1 (p.Met1Thr) MAF = 0.239"
genes[56]="OR8K1 (p.Met1Thr) MAF = 0.239"
genes[19]="NARS2 (p.Met1Ser)"
genes[43]="POLG2 (p.Met1fs)"
genes[38]="GP6 (p.Met1fs)"
genes[47]="ZNF827 (p.Met1Val)"
genes[57]="SAAL1 (p.Met1fs)"
genes[60]="ZNF682 (p.Met1fs)"
genes[55]="NUDT11 (p.Met1fs)"
# get dbsnp ids
dbsnp_ids = c()
for (file in exome_files){
  dbsnp_ids = c(dbsnp_ids, as.vector(get_dbsnp_ids(exome_tables[[file]])))
}

# Numero de genes afectados entre todos los pacientes.
length(unique(genes))

# gene frequency
op = par(las=2) # make label text perpendicular to axis
op = par(mar=c(5,18,4,2)) # increase y-axis margin.
barplot(table(genes), horiz=T, axes=F, xlab = "Number of carriers")
lablist<-as.vector(c(0:5))
axis(1, at=0:5, labels = FALSE)
text(seq(0, 5, by=1), par("usr")[3] - 0.2, labels = lablist, srt = 0, pos = 1, offset = 0.65, xpd = TRUE)
par(op)
View(exome_tables[["14-173.xlsx"]])
View(exome_tables[["2064.xlsx"]])
View(exome_tables[["ABGP.xlsx"]])
View(exome_tables[["Exoma 2166 nuevo.xlsx"]])
View(exome_tables[["Paqui.xlsx"]])
p1found=(search(csv, exome_tables[["14-173.xlsx"]]))
p2found=(search(csv, exome_tables[["2064.xlsx"]]))
p3found=(search(csv, exome_tables[["ABGP.xlsx"]]))
p4found=(search(csv, exome_tables[["Exoma 2166 nuevo.xlsx"]]))
p5found=(search(csv, exome_tables[["Paqui.xlsx"]]))
mergedVariationsFound = unique(rbind(p1found,p2found,p3found,p4found,p5found))
View(p1found)
View(p2found)
View(p3found)
View(p4found)
View(p5found)
View(data.frame(p1found,p2found))
  

get_gene_names(exome_tables[["14-173.xlsx"]])
get_gene_names(exome_tables[["2064.xlsx"]])
get_gene_names(exome_tables[["ABGP.xlsx"]])
get_gene_names(exome_tables[["Exoma 2166 nuevo.xlsx"]])
get_gene_names(exome_tables[["Paqui.xlsx"]])
myvars=c("coverage", "VARIATION_TYPE")
View(csv[myvars])
# Numero de mutaciones encontradas individualmente
nrow(exome_tables[["14-173.xlsx"]])
nrow(exome_tables[["2064.xlsx"]])
nrow(exome_tables[["ABGP.xlsx"]])
nrow(exome_tables[["Exoma 2166 nuevo.xlsx"]])
nrow(exome_tables[["Paqui.xlsx"]])

# MAF mean
mafs = c(0.318,0.401,0.141,0.233,0.197,0.006,0.104,0.239,0.346,0.06,0.419,0.448,0.309,0.123,0.261,0.37,0.133,0.388,0.081,0.123,0.169)
mean(mafs)
# tests
abs_file_dir = paste(exomes_dir, "PRL.GATK.snp.annovar.hg19_multianno.xls", sep="")
nuevo = read.xls("/home/fabad/hemodonacion/data/exomes/PRL.GATK.snp.annovar.hg19_multianno.xls")
View(table(genes))
