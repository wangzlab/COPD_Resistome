#args<-commandArgs(TRUE)


# libraries ---------------------------------------------------------------

library(data.table)
library(stringr)
library(dplyr)
library(ggplot2)

# batch id -------
fileID_byDate <- "2021.4"


# path of files -----------------------------------------------------------

krakenAlignOutDir <- list.files("E://Kraken/output_xx", full.names = T, pattern = fileID_byDate)
ARGsOAPReadDir <- list.files("E://ARGsOAP",full.names = T, pattern = paste("ARGsOAP_out_", fileID_byDate, sep="") )

krakenAlnFiles <- list.files(krakenAlignOutDir, full.names = T, pattern = "kraken_align.out$")
ARGsOAPReadFiles <- list.files(ARGsOAPReadDir,full.names = T,pattern = "blast6out.txt$")

taxidMappingDir <- "E://Kraken/3.out_taxid_mapping"
taxidMappingFile <- list.files(taxidMappingDir,full.names = T, pattern = paste(fileID_byDate,"_out_taxid_mapping",sep = ""))
taxidMapping_df <- fread(taxidMappingFile,data.table = F) #a lot of unidentified kingdom, needs other ways to map to kingdom


taxRankMappingDir <- "E://Kraken/4.out_taxRankingMappingDf"
taxRankMappingFile <- list.files(taxRankMappingDir, full.names = T, pattern = fileID_byDate)
taxRankMapping_df <- fread(taxRankMappingFile,data.table = F)


absentTaxID_File <- list.files(taxidMappingDir,full.names = T, pattern = paste(fileID_byDate,"_absentIDs",sep = ""))
absent_taxids <- fread(absentTaxID_File);absent_taxids<-absent_taxids$V1

ARGsOAPMappingFile <- "E:/ARGsOAP/SARG_mapping_file/SARG_mapping.out"
ARGsOAPMapping_df <- fread(ARGsOAPMappingFile,data.table = F)

saveToPath <- "E://Kraken/5_2.parsed_out_RData_ARGsOAP/"
TaxnomicLevels <- c("domain","kingdom", "phylum","class","order","family","genus","species")


# data processing ---------------------------------------------------------
#f1=as.integer(args[1])
#f2=as.integer(args[2])

ARGsOAPAln_df_batch <-fread(ARGsOAPReadFiles,data.table = F,
      col.names = c("query","target","pident","alnLen","misMatch","gap","queryStart","queryEnd","targetStart","targetEnd","evalue","bitScore"))

ARGsOAPAln_df_batch <- ARGsOAPAln_df_batch %>% 
  filter(alnLen>=25) %>% filter(pident >= 80) %>% filter(evalue <= 1e-7) %>% 
  filter(!duplicated(query))  

f1=1
f2=length(krakenAlnFiles)

krakenAln_df_batch <- NULL
for(krakenAlnFile in krakenAlnFiles[f1:f2]){
  # krakenAlnFile = krakenAlnFiles[17]
  writeLines(paste("current sample/file is: ", krakenAlnFile,sep = ""))
  
  krakenAln_df  <- fread(krakenAlnFile,select=c(2:4))
  
  krakenAln_df$taxid <- sub("^.*\\(taxid\\s(\\d+)\\)","\\1",krakenAln_df$V3) #extract taxid
  krakenAln_df$tax <- sub("^(.*)\\s\\(taxid\\s\\d+\\)","\\1",krakenAln_df$V3) #extract tax from alignment out file
  colnames(krakenAln_df)[which(colnames(krakenAln_df)=="V2")] <- "read_id"
  
  
  sampleName <- sub("(\\S+)\\.extracted.*$","\\1",basename(krakenAlnFile))

  #ARGsOAPAln_df <- ARGsOAPAln_df_batch  %>% select(query, target)
  
  
  #ARGsOAPAln_df$read_id <-gsub("\\.part\\-","part",ARGsOAPAln_df$read_id)
  #krakenAln_df$read_id <- gsub("\\.part\\-","part",krakenAln_df$read_id)
  
  # simplify krakenAln_df to just tax-annotate the reads hit by ARGsOAP
  krakenAln_df <-krakenAln_df %>% dplyr::filter(read_id %in% ARGsOAPAln_df_batch$query)
  
  krakenAln_df_batch <- bind_rows(krakenAln_df_batch, krakenAln_df)
  
}
  

#taxa annotation for krakenAln_df_batch taxids
krakenAln_df_batch <- cbind.data.frame(read_id = krakenAln_df_batch$read_id,
                                 taxid= krakenAln_df_batch$taxid,
                                 tax=krakenAln_df_batch$tax,
                                 domain = vector("character",nrow(krakenAln_df_batch)),
                                 kingdom = vector("character",nrow(krakenAln_df_batch)),
                                 phylum = vector("character",nrow(krakenAln_df_batch)),
                                 class = vector("character",nrow(krakenAln_df_batch)),
                                 order = vector("character",nrow(krakenAln_df_batch)),
                                 family = vector("character",nrow(krakenAln_df_batch)),
                                 genus = vector("character",nrow(krakenAln_df_batch)),
                                 species = vector("character",nrow(krakenAln_df_batch)),
                                 stringsAsFactors =F)

# first map by taxid mapping df, then complete higher tax level info by taxrank mapping df

for(i in c(1:nrow(krakenAln_df_batch))){
  #i=1
  if(i%%1000 == 0) writeLines(paste(i, "out of", nrow(krakenAln_df_batch), "in the for loop", sep = " "))
  tid <- as.integer(krakenAln_df_batch$taxid[i])
  if(tid == 0 | tid %in% absent_taxids){
    next
  }else{
    i_idmap <- which(taxidMapping_df$Original_query_taxid == tid)
    #from taxidMapping_df extract all information 
    krakenAln_df_batch[i, TaxnomicLevels] <-  taxidMapping_df[i_idmap,TaxnomicLevels]
    if(all(krakenAln_df_batch[i, TaxnomicLevels] == "<not present>")) next
    
    highestlvl = colnames(krakenAln_df_batch[i,TaxnomicLevels])[min(which(krakenAln_df_batch[i,TaxnomicLevels] != "<not present>"))] #highest taxo level which is not NA
    
    #from taxRankMapping_df extract information of higher tax levels in case higher levels has "<not present>"
    i_rkmap <- which(taxRankMapping_df[,highestlvl]==taxidMapping_df[i_idmap, highestlvl])[1]
    if(length(i_rkmap)>0){
      for(j_rkmap  in c((which(colnames(taxRankMapping_df)== highestlvl)-1):2)){
        krakenAln_df_batch[i, colnames(taxRankMapping_df)[j_rkmap]] <- taxRankMapping_df[i_rkmap,j_rkmap]
      }
    }
    
  }# this mapping process is much better than the first version (for reads)
  
}#for loop through nrow(krakenAln_df_batch)

  
  
  
  
ARGsOAPAln_df_batch <- as.data.frame(ARGsOAPAln_df_batch,stringsAsFactors=F)
sample_Alndf <- merge(ARGsOAPAln_df_batch %>% select(query, target), 
                      krakenAln_df_batch, by.x = "query", by.y="read_id") 

#sample_Alndf$sample_name <- sub("^(\\S+\\.decontam).*$","\\1",sample_Alndf$query)  #batch 2020.7用到的
sample_Alndf$sample_name <- sub("^(.*)_\\d+$","\\1",sample_Alndf$query)
  

#test = sample_Alndf
sample_Alndf[sample_Alndf == "<not present>"] <- NA
sample_Alndf[sample_Alndf == ""] <- NA



if(!dir.exists(saveToPath)) dir.create(saveToPath)
save( sample_Alndf, ARGsOAPMapping_df, taxidMapping_df, taxRankMapping_df, fileID_byDate, 
      file = paste(saveToPath,fileID_byDate,".RData",sep=""))
  
  
