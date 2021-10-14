# collect all the NCBI Taxonomy ID from kraken output


KrakenDir <- "E://Kraken/output_xx"

krakenFiles <- list.files(KrakenDir,full.names = T)[!dir.exists(list.files(KrakenDir,full.names = T))] 

#krakenFile <- krakenFiles[grepl(sampleName, krakenFiles)]


library(data.table)
all_taxids <- vector("character")
for(krakenFile in krakenFiles){
  data <- fread(krakenFile,select = c(2:3),data.table = F)
  taxids <- unique(sub("^.*\\(taxid\\s(\\d+)\\)","\\1",data$V3))
  all_taxids <- unique(c(all_taxids,taxids))
}

write.table(data.frame(all_taxids,stringsAsFactors = F), "1.out_all_taxids/2021.4.txt", row.names = F, quote = F)
