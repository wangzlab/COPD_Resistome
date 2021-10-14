

# Libraries ---------------------------------------------------------------
library(data.table)
library(dplyr)

# batch id (change every time) -----------------
fileID_byDate <- "2021.4"

# File paths --------------------------------------------------------------
parsedOutDir <- "E://xuexue_data/Kraken/5_2.parsed_out_RData_ARGsOAP/"
RDataFiles <- list.files(parsedOutDir,full.names = T, pattern = fileID_byDate)


ARGsOAP_outDir <- paste("E://xuexue_data/ARGsOAP/ARGsOAP_out_",fileID_byDate,sep = "")
totalReadCountFile <- list.files(ARGsOAP_outDir,full.names = T,pattern = "numOfReads")
totalReadCount_df <- fread(totalReadCountFile,data.table = F)


geneLengthFile <- "E:/ARGsOAP/SARG_mapping_file/SARG.geneLength.txt"
geneLength_df <- fread(geneLengthFile, data.table = F)
geneLength_df$ARGsubtype <- sapply(strsplit(geneLength_df$ARG,"__",fixed = T),"[[",2)


# Data processing ---------------------------------------------------------
taxas <- c("domain", "kingdom", "phylum", "class","order","family","genus","species")


load(RDataFiles)

sample_Alndf$sample_name <-sub("^(.*)\\.part\\-\\d$","\\1",  sub("^(.*)_\\d+$","\\1",sample_Alndf$query)) # if contain part-*, remove

for(taxa in taxas){
  # taxa=taxas[8]
  
  dat <- sample_Alndf[,c("query","target",taxas,"sample_name")] #keep sample_Alndf unchanged
  #dat$sample_name <- sapply(strsplit(dat$query,"_",fixed = T),"[[",1 ) #前面generate sample.RData的环节可能没有弄好，sample_name是错的
  
  #colnames(dat)[which(colnames(dat)=="best-hit")] <-"best_hit"
  colnames(dat)[which(colnames(dat)==taxa)] <- "taxa"
  
  
  ARG_taxa_countDf <- dat %>% dplyr::group_by(sample_name,target,taxa) %>% dplyr::summarise(readCount = n()) %>% as.data.frame()
  ARG_taxa_countDf$totalReadCount_inSample <- sapply(ARG_taxa_countDf$sample_name, 
                                                     function(x) totalReadCount_df$numReads[which(totalReadCount_df$sample_name == x)])
  
  ARG_taxa_countDf$gene_length <- sapply(ARG_taxa_countDf$target, function(x) geneLength_df$len_nt[which(geneLength_df$geneID == x)][1])
  ARG_taxa_countDf$ppm <- ARG_taxa_countDf$readCount/(ARG_taxa_countDf$totalReadCount_inSample/1000000)
  ARG_taxa_countDf$rpkm <- ARG_taxa_countDf$ppm/(ARG_taxa_countDf$gene_length/1000) #calculate based on best-hit because one ARG subtype still have multiple sequences with different gene lengths
  
  ARG_taxa_countDf$ARG <- sapply(ARG_taxa_countDf$target, function(x) geneLength_df$ARGsubtype[which(geneLength_df$geneID == x)][1]) 
  ARG_taxa_countDf <- ARG_taxa_countDf %>% 
    dplyr::group_by(sample_name, ARG,taxa) %>% 
    dplyr::summarise(readCount = sum(readCount), ppm=sum(ppm), rpkm=sum(rpkm))
  
  ARG_taxa_countDf$drug_type <- sapply(ARG_taxa_countDf$ARG, 
                                       function(x) ARGsOAPMapping_df$ARG_type[which(ARGsOAPMapping_df$ARG_subtype == x)[1]]) 
  ARG_taxa_countDf$res_mechanism <- sapply(ARG_taxa_countDf$ARG, 
                                           function(x) ARGsOAPMapping_df$mechanism5[which(ARGsOAPMapping_df$ARG_subtype == x)[1]]) 
  
  
  #colnames(ARG_taxa_countDf)[which(colnames(ARG_taxa_countDf)=="taxa")]<-taxa
  assign(paste("all_ARG_",taxa,"_df",sep = ""), ARG_taxa_countDf, envir = .GlobalEnv)
  
  
} #for loop through Taxas




# organize by gene (no taxon) -------------------------------------
all_ARG_byGene_df <- 
  all_ARG_family_df %>% group_by(sample_name, ARG) %>% 
  summarise(readCount= sum(readCount), ppm = sum(ppm), rpkm=sum(rpkm), drug_type=unique(drug_type), res_mechanism= unique(res_mechanism))


all_ARG_byDrugType_df <-
  all_ARG_family_df %>% group_by(sample_name, drug_type) %>%
  summarise(readCount = sum(readCount), ppm = sum(ppm), rpkm=sum(rpkm))


all_ARG_byResMechanism_df <- 
  all_ARG_family_df %>% group_by(sample_name, res_mechanism) %>%
  summarise(readCount = sum(readCount), ppm = sum(ppm), rpkm=sum(rpkm))


save(all_ARG_byDrugType_df, all_ARG_byGene_df, all_ARG_byResMechanism_df, 
     all_ARG_class_df, all_ARG_domain_df, all_ARG_family_df, all_ARG_genus_df, all_ARG_kingdom_df, all_ARG_order_df, all_ARG_phylum_df, all_ARG_species_df,
     file = paste( "6.out_integratedRData/", fileID_byDate, "_integrated_ARGsOAP.RData", sep = ""))
