library(data.table)
library(stringr)
library(dplyr)

fileID_byDate <- "2021.4" 


# create mapping from kraken count file 
krakenCountDir <- "E://Kraken/report_xx"

krakenCountFiles <- list.files(krakenCountDir,full.names = T,pattern = "taxa_count.tsv$")


# tax mapping from kraken count file -----------------------------------------------
taxRankMapping_completeDf <- data.frame()
for(krakenCountFile in krakenCountFiles){
  #krakenCountFile = krakenCountFiles[1]
  writeLines(paste("now processing:",which(krakenCountFiles == krakenCountFile), krakenCountFile,sep = " "))
  
  taxRankMapping_df <- fread(krakenCountFile, select = 1)
  
  taxRankMapping_df$domain <- sapply(taxRankMapping_df$V1,
                                     function(x) {
                                       if(any(grepl("d__",strsplit(x,"|",fixed = T)[[1]]))){
                                         txt = strsplit(x,"|",fixed = T)[[1]][grepl("d__",strsplit(x,"|",fixed = T)[[1]])]
                                         gsub("^d__(.*)","\\1",txt)
                                       }else{
                                         NA
                                       }
                                     })
  taxRankMapping_df$kingdom <- sapply(taxRankMapping_df$V1,
                                      function(x) {
                                        if(any(grepl("k__",strsplit(x,"|",fixed = T)[[1]]))){
                                          txt=strsplit(x,"|",fixed = T)[[1]][grepl("k__",strsplit(x,"|",fixed = T)[[1]])]
                                          sub("^k__(.*)","\\1",txt)
                                        }else{
                                          NA
                                        }
                                      })
  taxRankMapping_df$phylum <- sapply(taxRankMapping_df$V1,
                                     function(x) {
                                       if(any(grepl("p__",strsplit(x,"|",fixed = T)[[1]]))){
                                         txt=strsplit(x,"|",fixed = T)[[1]][grepl("p__",strsplit(x,"|",fixed = T)[[1]])]
                                         sub("^p__(.*)","\\1",txt)
                                       }else{
                                         NA
                                       }
                                     })
  
  taxRankMapping_df$class <- sapply(taxRankMapping_df$V1,
                                    function(x) {
                                      if(any(grepl("c__",strsplit(x,"|",fixed = T)[[1]]))){
                                        txt = strsplit(x,"|",fixed = T)[[1]][grepl("c__",strsplit(x,"|",fixed = T)[[1]])]
                                        sub("^c__(.*)","\\1",txt)
                                      }else{
                                        NA
                                      }
                                    })
  
  taxRankMapping_df$order <- sapply(taxRankMapping_df$V1,
                                    function(x) {
                                      if(any(grepl("o__",strsplit(x,"|",fixed = T)[[1]]))){
                                        txt = strsplit(x,"|",fixed = T)[[1]][grepl("o__",strsplit(x,"|",fixed = T)[[1]])]
                                        sub("^o__(.*)","\\1",txt)
                                      }else{
                                        NA
                                      }
                                    })
  taxRankMapping_df$family <- sapply(taxRankMapping_df$V1,
                                     function(x){
                                       if(any(grepl("f__",strsplit(x,"|",fixed = T)[[1]]))){
                                         txt=strsplit(x,"|",fixed = T)[[1]][grepl("f__",strsplit(x,"|",fixed = T)[[1]])]
                                         sub("^f__(.*)","\\1",txt)
                                       }else{
                                         NA
                                       }
                                     } )
  taxRankMapping_df$genus <- sapply(taxRankMapping_df$V1,
                                    function(x){
                                      if(any(grepl("g__",strsplit(x,"|",fixed = T)[[1]]))){
                                        txt=strsplit(x,"|",fixed = T)[[1]][grepl("g__",strsplit(x,"|",fixed = T)[[1]])]
                                        sub("^g__(.*)","\\1",txt)
                                      }else{
                                        NA
                                      }
                                    })
  taxRankMapping_df$species <- sapply(taxRankMapping_df$V1,
                                      function(x) {
                                        if(any(grepl("s__",strsplit(x,"|",fixed = T)[[1]]))){
                                          txt=strsplit(x,"|",fixed = T)[[1]][grepl("s__",strsplit(x,"|",fixed = T)[[1]])]
                                          sub("^s__(.*)","\\1",txt)
                                        }else{
                                          NA
                                        }
                                      })
  
  taxRankMapping_df <- taxRankMapping_df %>% as.data.frame()%>% unique() 
  
  taxRankMapping_completeDf <- rbind.data.frame(taxRankMapping_completeDf, taxRankMapping_df,stringsAsFactors = F) %>% unique()
}

write.table(taxRankMapping_completeDf, file = paste("4.out_taxRankingMappingDf/",fileID_byDate,".txt",sep = ""), 
            quote = F, row.names = F,sep = "\t")
