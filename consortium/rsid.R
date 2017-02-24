source("https://bioconductor.org/biocLite.R")
library("BSgenome")
#biocLite("SNPlocs.Hsapiens.dbSNP141.GRCh38")
library("SNPlocs.Hsapiens.dbSNP141.GRCh38")

file.list=list.files(path="../Data/significant/",pattern = "in.txt")
for(i in 12:length(file.list)){
  dat=read.table(paste0("../Data/significant/",file.list[i]))
  rsid = rownames(dat)
  snps = snpid2loc(SNPlocs.Hsapiens.dbSNP141.GRCh38, rsid)
  chroms = gsub('ch([0-9]+)','\\1', names(snps))
  positions = snps
  out = data.frame(rsid, chroms, positions, lfsr = dat[, 1],
                   stringsAsFactors = FALSE)
  write.table(out, file = paste0("../Data/significant/pos", file.list[i]),
              quote = FALSE, row.names = FALSE)
}


##for i =11 and i =16

for(i in 2:length(file.list)){
  dat=read.table(paste0("../Data/significant/",file.list[i]))
  rsid = rownames(dat)
  
  
  snps = snpsById(SNPlocs.Hsapiens.dbSNP141.GRCh38, rownames(dat), 
                  ifnotfound = 'warning')
  chroms = as.character(sapply(snps, function(x) x@seqnames@values))
  chroms = gsub('ch([0-9]+)','\\1', chroms)
  positions = as.character(sapply(snps, function(x) x@ranges@start))
  
  
  
  out = data.frame(rsid, chroms, positions, lfsr = dat[, 1],
                   stringsAsFactors = FALSE)
  write.table(out, file = paste0("../Data/significant/pos-", file.list[i]),
              quote = FALSE, row.names = FALSE)
}
