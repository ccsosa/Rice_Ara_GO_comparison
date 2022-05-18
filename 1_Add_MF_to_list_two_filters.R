#CALLING LIBRARIES
require(riceidconverter)
require(readxl)
require(tidyr)
require(biomaRt)
require(rentrez)
require(xlsx)
require(plot3D)
require(ggpmisc)
require(gprofiler2)
#require(clusterProfiler)
require(GOCompare)
require(GOfuncR)
require(agricolae)
require(GO.db)
#require(GOSim)
require(dplyr)
require(stringr)
require(topGO)
#require(pivottabler)
require(BSDA)
#require(ggheatmap)
require(pheatmap)
require(treemap)
require(bio3d)
require (factoextra)
require (NbClust)
require (cluster)
require (clustertend)
#require(ViSEAGO)
require(ggpubr)
#require(GOSemSim)

library(rrvgo)
dir <- "D:/TESIS_PHD/CHAPTER1"


dat3 <- as.data.frame(readxl::read_xlsx(paste0(dir,"/","LISTA_GENES_ALUMINIO_2_NO_DUPS_riceconverter_to.xlsx"),
                                        col_names=T))

dat_dups <- as.data.frame(readxl::read_xlsx(paste0(dir,"/","LISTA_GENES_ALUMINIO_2_DUPS.xlsx"),
                                            col_names=T))
dat_dups$...1 <- NULL
colnames(dat_dups)[10] <- "ci_gen"

dat_dups$PLAZA <- NA
dat_dups[,"Candidate gene ID_TRIM...9"] <- dat_dups$`Candidate gene ID_TRIM...10`
dat_dups <- dat_dups[,colnames(dat3)]


x_dups <- unique(dat_dups$ci_gen)
x_dups_list  <- list()

#ids
x <- riceidconverter::RiceIDConvert(trimws(dat_dups$ci_gen),
                                    fromType = "MSU",
                                    toType = "RAP")
for(i in 1:length(x_dups)){
  message(i)
  x_new <- x$RAP[x$MSU==x_dups[[i]]]

  x_new <- bind_rows(replicate(length(x_new),
                               dat_dups[which(dat_dups$ci_gen==x_dups[[i]]),],
                               simplify = FALSE))
  x_new$ci_gen <- x$RAP[x$MSU==x_dups[[i]]]
  x_dups_list[[i]] <- x_new
};rm(i)
x_dups_list <- do.call(rbind,x_dups_list)
x_dups_list$tool_convert <- "DUPS"

###############################################################


dat3 <- rbind(dat3,x_dups_list)

x_ori_ids <- strsplit(dat3$`Candidate gene ID_TRIM...9`,"[.]")
x_ori_ids <- trimws(lapply(x_ori_ids, `[[`, 1))



dat3$MSU <- NA

for(i in 1:nrow(dat3)){
  #message(i)
  if(is.na(str_detect(dat3$ci_gen[[i]],"LOC_"))){
    dat3$MSU[[i]] <- NA
  } else if(str_detect(dat3$ci_gen[[i]],"LOC_")==T){
    dat3$MSU[[i]] <- dat3$ci_gen[[i]]
  } else  {
    if(str_detect(dat3$`Candidate gene ID_TRIM...9`[[i]],"LOC_")==T){
      dat3$MSU[[i]] <- x_ori_ids[[i]]
    } else {
      dat3$MSU[[i]] <- NA
    }
  }
};rm(i)



unique_x <- unique(x$RAP)
unique_x_list <- list()
for(i in 1:length(unique_x)){
  xx_str <- str_subset(x$MSU[which(x$RAP==unique_x[[i]])], "[.]", negate = TRUE)
  if(length(xx_str)>0){
    unique_x_list[[i]] <- data.frame(RAP=unique_x[[i]], MSU=xx_str)
  } else {
    unique_x_list[[i]] <- data.frame(RAP=unique_x[[i]], MSU=NA)
  }

};rm(i)

unique_x_list <- do.call(rbind,unique_x_list)
unique_x_list <- unique_x_list[which(unique_x_list$MSU!="None"),]
unique_x_list <- unique_x_list[which(!is.na(unique_x_list$MSU)),]



for(i in 1:nrow(unique_x_list)){
dat3$MSU[which(dat3$ci_gen==unique_x_list$RAP[[i]])] <- unique_x_list$MSU[[i]]
};rm(i)

 # write.xlsx(dat3, paste0(dir,"/",'LISTA_GENES_ALUMINIO_3a2.xlsx'),
 #           sheetName = "Sheet1",showNA = F,row.names = F)



#mol_func <- as.data.frame(read.csv(paste0(dir,"/MF/","GO_Gropus.csv"),header = T))
#######################################################
#call biomaRt mirror

ensembl_pl = biomaRt::useMart(biomart="plants_mart",host="plants.ensembl.org")
# avail_datasets <- listDatasets(ensembl2)
ensembl_os = biomaRt::useDataset("osativa_eg_gene", mart = ensembl_pl)

#ensembl_oi = biomaRt::useDataset("oindica_eg_gene", mart = ensembl_pl)
ensembl_at = biomaRt::useDataset("athaliana_eg_gene", mart = ensembl_pl)

#######################################################
###first approach MF

x_cand <- strsplit(dat3$ci_gen,"[.]")

############################################################################################
#orthologs

x_orth_sat <- gorth(unique(unlist(x_cand)), source_organism="osativa", target_organism="athaliana",filter_na=T)
#x_orth_ind <- gorth(unlist(x_cand), source_organism="oindica", target_organism="athaliana",filter_na=T)
############################################################################################

#first filter (genes with orthologues in ARATH)
query_orth_sat <- unique(x_orth_sat$input)
#query_orth_ind <- unique(x_orth_ind$input)
dat3_Sat <- dat3[dat3$ci_gen %in% query_orth_sat,]
length(unique(dat3_Sat$ci_gen))
#dat3_Ind <- dat3[dat3$ci_gen %in% query_orth_ind,]
############################################################################################
#get GO term information for O. sativa subsp. japonica
go_ids= biomaRt::getBM(attributes=c('go_id', 'ensembl_gene_id', 'namespace_1003'),
                       filters='ensembl_gene_id',
                       values=unique(unlist(dat3_Sat$ci_gen)),
                       mart=ensembl_os)

go_ids <- go_ids[which(go_ids$namespace_1003=="molecular_function"),]
go_ids <- go_ids[,c(1,2)]

############################################################################################
#get go term information for A. sativa
go_ids_ARATH= biomaRt::getBM(attributes=c('go_id', 'ensembl_gene_id', 'namespace_1003'),
                       filters='ensembl_gene_id',
                       values=unique(unlist(x_orth_sat$ortholog_ensg)),

                                             mart=ensembl_at)

go_ids_ARATH <- go_ids_ARATH[which(go_ids_ARATH$namespace_1003=="molecular_function"),]
go_ids_ARATH <- go_ids_ARATH[,c(1,2)]
#############################################################################################
#second filter (information available in ARATH and OSJ)
x_orth_sat2 <- x_orth_sat[x_orth_sat$ortholog_ensg %in% unique(go_ids_ARATH$ensembl_gene_id),]

dat3_Sat <- dat3_Sat[dat3_Sat$ci_gen %in% unique(x_orth_sat2$input),]
unique(dat3_Sat$ci_gen)
############################################################################################
#save information
write.xlsx(dat3_Sat, paste0(dir,"/",'LISTA_GENES_ALUMINIO_3_Sat.xlsx'),
           sheetName = "Sheet1",showNA = F,row.names = F)

write.xlsx(x_orth_sat2, paste0(dir,"/",'ORTH_SAT.xlsx'),
           sheetName = "Sheet1",showNA = F,row.names = F)
############################################################################################
#get information from ensembl after second filter round
go_ids= biomaRt::getBM(attributes=c('go_id', 'ensembl_gene_id', 'namespace_1003'),
                       filters='ensembl_gene_id',
                       values=unique(dat3_Sat$ci_gen),#unlist(dat3_Sat$ci_gen),
                       mart=ensembl_os)

go_ids <- go_ids[which(go_ids$namespace_1003=="biological_process"),]
go_ids <- go_ids[,c(1,2)]


message("unique OSJ genes: ",length(unique(dat3_Sat$ci_gen)) )
message("unique ARATH genes: ",length(unique(x_orth_sat2$ortholog_ensg)) )
############################################################################################

x_s <-  gprofiler2::gost(query = unique(dat3_Sat$ci_gen),
                         organism = "osativa", ordered_query = FALSE,
                         multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                         measure_underrepresentation = FALSE, evcodes = T,
                         user_threshold = 0.01, correction_method = "g_SCS",
                         #domain_scope = "annotated",
                         custom_bg = NULL,#unique(dat3_Sat$ci_gen),
                         numeric_ns = "", sources = "GO:BP", as_short_link = FALSE)

write.table(as.data.frame(x_s$result[,-c(14)]),paste0(dir,"/MF/TopGO_2_filters/","gprofiler_results.tsv"),na = "",sep="\t",row.names = F,quote = F)
  #x_s$result[,c(1:13)]),paste0(dir,"/MF/TopGO_2_filters/","gprofiler_results.tsv"),na = "",sep="\t",row.names = F,quote = F)

############################################################################################
#preparing data for elim algorithm
# gene_2_GO <- unstack(go_ids[,c(1,2)])
# keep  <- unique(dat3_Sat$ci_gen)  %in% go_ids[,2]#unlist(dat3_Sat$ci_gen) %in% go_ids[,2]
# keep  <- which(keep==TRUE)
# candidate_list <- unique(dat3_Sat$ci_gen)[keep] #unlist(dat3_Sat$ci_gen)[keep]
# geneList <- factor(as.integer(unique(dat3_Sat$ci_gen) %in% candidate_list),levels = c(0,1))
# #geneList=factor(as.integer(unlist(dat3_Sat$ci_gen) %in% candidate_list),levels = c(0,1))
# names(geneList)= unique(dat3_Sat$ci_gen)  #unlist(dat3_Sat$ci_gen)
# summary(geneList)
# #Create topGOdata object
# GOdata=new('topGOdata', ontology='BP', allGenes = geneList, annot = topGO::annFUN.gene2GO,description="OSJ_BP groups",
#            gene2GO = gene_2_GO)
# #https://ucdavis-bioinformatics-training.github.io/2019_August_UCD_mRNAseq_Workshop/differential_expression/enrichment
# #resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")#, cutOff=0.05)
# resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")#, cutOff=0.05)
# #pks <- data.frame(p.adjust(resultKS.elim@score,method = "fdr"))
# sum(resultKS.elim@score <0.05)
#
#
# #https://ucdavis-bioinformatics-training.github.io/2019_August_UCD_mRNAseq_Workshop/differential_expression/enrichment
# #showSigOfNodes(GOdata, score(resultKS.elim), firstSigNodes =27, useInfo = "def")
#
# tab <- GenTable(GOdata, raw.p.value = resultKS.elim, topNodes = length(resultKS.elim@score), numChar = 180)
# tab$p_adj <- p.adjust(tab$raw.p.value,method = "fdr")
# write.table(tab[,c(1,6)][which(tab$raw.p.value<0.05),],paste0(dir,"/MF/TopGO_2_filters/","tab_go_bp.tsv"),na = "",sep="\t",row.names = F,quote = F)
# write.table(tab,paste0(dir,"/MF/TopGO_2_filters/","tab_go_bp_all.tsv"),na = "",sep="\t",row.names = F,quote = F)
#
#
# # tab$FDR <- p.adjust(tab$raw.p.value,method = "fdr")
# # par(cex = 0.3)
# # showSigOfNodes(GOdata, score(resultKS.elim), firstSigNodes = 5, useInfo = "def")
# #p values elim
# score_elim <-  score(resultKS.elim)[score(resultKS.elim) <0.05]
# allGO <- genesInTerm(GOdata)
# x_all <- allGO[names(score_elim)]

######
#gprofiler

x_all <- list()
for(i in 1:nrow(x_s$result)){
x_all[[i]] <- unlist(strsplit(x_s$result$intersection[[i]],","))
}
names(x_all) <- x_s$result$term_id
#https://ucdavis-bioinformatics-training.github.io/2019_August_UCD_mRNAseq_Workshop/differential_expression/enrichment

############################################################################################
###MAP FUNCTIONS IN SATIVA JAPONICA
TERM_GO_list <- list()
for(i in 1:length(x_all)){
  TERM_GO_list[[i]]  <- data.frame(TERM=Term(GOTERM[names(x_all[i])]),GO_ID=names(x_all[i]),genes=x_all[[i]])
};rm(i)

TERM_GO_list <- do.call(rbind,TERM_GO_list)
no_found_df_sat <- data.frame(TERM="UNKNOWN",GO_ID=NA,genes=dat3_Sat$ci_gen[!dat3_Sat$ci_gen %in% TERM_GO_list$genes])

#checking UNKNOWN genes are correct!
intersect((TERM_GO_list$genes),
          (no_found_df_sat$genes))
TERM_GO_list <- rbind(TERM_GO_list,no_found_df_sat)

############################################################################################
# #filter groups
 unique_MF <- unique(TERM_GO_list$TERM) #get uniques

############################################################################################
#loading revigo reults


revigo_data <- data.table::fread(paste0(dir,"/MF/TopGO_2_filters/","revigo.tsv"))#csv"))
revigo_data$final_group <- NA


for(i in 1:nrow(revigo_data)){
  if(revigo_data$Eliminated[[i]]==FALSE){
    revigo_data$final_group[[i]] <- revigo_data$Name[[i]]
  } else {
    revigo_data$final_group[[i]] <- revigo_data$final_group[[i-1]]
  }

};rm(i)

write.table(as.data.frame(revigo_data),paste0(dir,"/MF/TopGO_2_filters/","revigo_data_final.tsv"),na = "",sep="\t",row.names = F,quote = F)



#assign final groups
TERM_GO_list$TERM_FILTERED <- NA

for(i in 1:length(revigo_data$Name)){
  TERM_GO_list$TERM_FILTERED[which(TERM_GO_list$TERM==revigo_data$Name[[i]])] <- revigo_data$final_group[[i]]
};rm(i)


TERM_GO_list$TERM_FILTERED[which( TERM_GO_list$TERM=="UNKNOWN")] <- "UNKNOWN"
length(unique(TERM_GO_list$genes))

unique(TERM_GO_list$TERM_FILTERED)
TERM_GO_list <- TERM_GO_list[which(TERM_GO_list$TERM_FILTERED !="biological_process"),]
length(unique(TERM_GO_list$genes))
#https://www.bioconductor.org/packages/release/bioc/vignettes/rrvgo/inst/doc/rrvgo.html
# library(rrvgo)
# tab$raw.p.value <- as.numeric(tab$raw.p.value)
# tab2 <- tab[which(tab$raw.p.value<0.05),]
# simMatrix <- calculateSimMatrix(tab2$GO.ID,
#                                 orgdb="org.At.tair.db",
#                                 ont="BP",
#                                 method="Rel")
#
# scores <- setNames(-log10(as.numeric(tab2$raw.p.value)), tab2$GO.ID)
# reducedTerms <- reduceSimMatrix(simMatrix,
#                                 scores,
#                                 threshold=0.7,
#                                 orgdb="org.At.tair.db")
#unique(reducedTerms$parentTerm)
#scatterPlot(simMatrix, reducedTerms)


unique(TERM_GO_list$TERM_FILTERED)
unique(revigo_data$Name)
unique(revigo_data$final_group)

############################################################################################
#save data
write.csv(TERM_GO_list,paste0(dir,"/MF/TopGO_2_filters/","BP_Groups_elim2.csv"),na = "",row.names = F)

unique_MF <- unique(TERM_GO_list$TERM_FILTERED)

# <- unique_MF[-5]#[-3]

x_MF_list <- list()
for(j in 1:length(unique_MF)){
  x_MF_list[[j]] <- as.list(unique(TERM_GO_list$genes[which(TERM_GO_list$TERM_FILTERED==unique_MF[[j]])]))
};rm(j)

names(x_MF_list) <- unique_MF

#get genes with info
length(unique(unlist(x_MF_list)))
################################################################################################################
############################################################################################
#get information from ensembl after second filter round
# go_ids= biomaRt::getBM(attributes=c('go_id', 'ensembl_gene_id', 'namespace_1003'),
#                        filters='ensembl_gene_id',
#                        values=unique(dat3_Sat$ci_gen),#unlist(dat3_Sat$ci_gen),
#                        mart=ensembl_os)
#
# go_ids <- go_ids[which(go_ids$namespace_1003=="molecular_function"),]
# go_ids <- go_ids[,c(1,2)]
#
# gene_2_GO=unstack(go_ids[,c(1,2)])
# #all GO terms lists to query
# xx <- as.list(GOTERM)
#
############################################################################################
#preparing data for elim algorithm OSJ

#i <- 1
x_i_go <- lapply(1:length(x_MF_list),function(i){
message(i)
  x_s_i <-  gprofiler2::gost(query = unique(unlist(x_MF_list[[i]])),
                           organism = "osativa", ordered_query = FALSE,
                           multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                           measure_underrepresentation = FALSE, evcodes = T,
                           user_threshold = 0.05, correction_method = "g_SCS",
                           #domain_scope = "annotated",
                           custom_bg = unique(dat3_Sat$ci_gen),#unique(dat3_Sat$ci_gen),
                           numeric_ns = "", sources = "GO:MF", as_short_link = FALSE)


  # candidate <- unique(unlist(x_MF_list[[i]])) #unlist(x_MF_list[[i]])
  # geneList=factor(as.integer(unique(dat3_Sat$ci_gen) %in% candidate),levels = c(0,1))
  # #geneList=factor(as.integer(unlist(dat3_Sat$ci_gen) %in% candidate),levels = c(0,1))
  # names(geneList)= unique(dat3_Sat$ci_gen)#unlist(dat3_Sat$ci_gen)
  #
  #
  # #Create topGOdata object
  # GOdata_i=new('topGOdata', ontology='MF', allGenes = geneList, annot = topGO::annFUN.gene2GO,
  #              gene2GO = gene_2_GO)
  #
  # resultKS.elim_i <- runTest(GOdata_i, algorithm = "elim", statistic = "ks")
  #
  # #p values elim
  #
  # score_elim_i <-  score(resultKS.elim_i)[score(resultKS.elim_i) <0.05]
  # allGO_i <- genesInTerm(GOdata_i)
  # x_all_i <- allGO_i[names(score_elim_i)]

  # Term_i <- as.data.frame(do.call(rbind,lapply(1:length(x_all_i),function(j){
  #   Term_i <- Term(xx_i <- xx[names(xx) %in% names(x_all_i)[j]][[1]])
  #   Term_i <- data.frame(term_name=Term_i,genes= paste(unlist(x_all_i[[j]]),collapse = ","))

      #return(Term_i)
#  })))

  if(!is.null(x_s_i)){
    Term_i <- data.frame(term=x_s_i$result$term_name,genes=x_s_i$result$intersection)

    x_i_go  <- data.frame(feature=names(x_MF_list)[[i]],term_name=Term_i$term,term_id =x_s_i$result$term_id,
                          #names(x_all_i),
                          genes=Term_i$genes)
  } else {
    x_i_go  <- data.frame(feature=names(x_MF_list)[[i]],term_name=NA,term_id =NA,
                        #names(x_all_i),
                        genes=NA)
  }
  return(x_i_go)
})

# length(unique(unlist(x_MF_list[[28]])))
# length(unique(unlist(x_MF_list[1:27])))

x_i_go <- do.call(rbind,x_i_go)
x_i_go2 <- x_i_go[which(x_i_go$term_name!="molecular_function"),]
write.csv(as.data.frame(x_i_go2),paste0(dir,"/MF/TopGO_2_filters/","GOST_OSJ.csv"),na = "",row.names = F)

message("Genes used for enrichment for OSJ: ",length(unique(unlist(x_MF_list))))
############################################################################################
length(unique(unlist(strsplit(x_i_go2$genes,",")))) #number of total genes
#length(unique(unlist(strsplit(x_i_go_arath2$genes,",")))) #number of total genes
#############################################################
#graph one species
#load_all("D:/REPO_GITHUB/GOCompare")
x <- graphGOspecies(df=x_i_go2,#x_s$result,
                               GOterm_field="term_name",
                               option = "GO",
                               numCores=6,
                               saveGraph=T,
                               outdir = paste0(dir,"/MF/TopGO_2_filters"))



x2 <- graphGOspecies(df=x_i_go2,#x_s$result,
                    GOterm_field="term_name",
                    option = "Categories",
                    numCores=6,
                    saveGraph=T,
                    outdir = paste0(dir,"/MF/TopGO_2_filters"))


#Get nodes with values greater than 95%
perc <- x$nodes[which(x$nodes$GO_WEIGHT >= quantile(x$nodes$GO_WEIGHT,probs = 0.9)),]
write.csv(perc,paste0(dir,"/MF/TopGO_2_filters/","PERC_OSJ.csv"),na = "",row.names = F)
# visualize nodes filtered
#View(perc)


###########################arabidopsos mapping




TERM_GO_list_ARA <- list()
for(i in 1:length(x_MF_list)){
v_i <- TERM_GO_list[which(TERM_GO_list$TERM_FILTERED==names(x_MF_list)[i]),]
v_i_i <- TERM_GO_list[which(TERM_GO_list$TERM==names(x_MF_list)[i]),]

TERM_GO_list_ARA[[i]] <- data.frame(TERM=NA,
                     GO_ID=unique(v_i_i$GO_ID),
                      genes=x_orth_sat2$ortholog_ensg[x_orth_sat2$input %in% unique(v_i$genes)],
                     TERM_FILTERED=unique(v_i$TERM_FILTERED))

  };rm(i)

TERM_GO_list_ARA <- do.call(rbind,TERM_GO_list_ARA)

no_found_df_ara <- data.frame(TERM="UNKNOWN",
                              GO_ID=NA,
                              genes=x_orth_sat2$ortholog_ensg[x_orth_sat2$input %in% no_found_df_sat$genes],
                              TERM_FILTERED ="UNKNOWN")
TERM_GO_list_ARA <- rbind(TERM_GO_list_ARA,no_found_df_ara)

write.csv(TERM_GO_list_ARA,paste0(dir,"/MF/TopGO_2_filters/","MF_Groups_elim_ARATH.csv"),na = "",row.names = F)
################################################
#get ARATH gene lists
x_MF_list_ara <- list()
unique_MF_ara <- unique(TERM_GO_list_ARA$TERM_FILTERED)
for(j in 1:length(unique_MF_ara)){
  x_MF_list_ara[[j]] <- as.list(TERM_GO_list_ARA$genes[which(TERM_GO_list_ARA$TERM_FILTERED==unique_MF_ara[[j]])])
};rm(j)

names(x_MF_list_ara) <- unique_MF_ara
############################################################################################
#get information from ensembl after second filter round
# go_ids= biomaRt::getBM(attributes=c('go_id', 'ensembl_gene_id', 'namespace_1003'),
#                        filters='ensembl_gene_id',
#                        values=unique(x_orth_sat2$ortholog_ensg),#unlist(x_orth_sat2$ortholog_ensg),
#                        mart=ensembl_at)
#
# go_ids <- go_ids[which(go_ids$namespace_1003=="molecular_function"),]
# go_ids <- go_ids[,c(1,2)]
#
# gene_2_GO=unstack(go_ids[,c(1,2)])
# #all GO terms lists to query
# xx <- as.list(GOTERM)

###########################################################################
#local enrichment (topGO) ARATH

x_i_go_arath <- lapply(1:length(x_MF_list_ara),function(i){
  # candidate <- unique(x_MF_list_ara[[i]])#unique(unlist(x_MF_list_ara[[i]]))
  # geneList=factor(as.integer(unique(x_orth_sat2$ortholog_ensg) %in% candidate),levels = c(0,1))
  # #geneList=factor(as.integer(unlist(x_orth_sat2$ortholog_ensg) %in% candidate),levels = c(0,1))
  # names(geneList)= unique(x_orth_sat2$ortholog_ensg)#unlist(x_orth_sat2$ortholog_ensg)
  #
  #
  # #Create topGOdata object
  # GOdata_i=new('topGOdata', ontology='MF', allGenes = geneList, annot = topGO::annFUN.gene2GO,
  #              gene2GO = gene_2_GO)
  #
  # resultKS.elim_i <- runTest(GOdata_i, algorithm = "elim", statistic = "ks")
  #
  # #p values elim
  # score_elim_i <-  score(resultKS.elim_i)[score(resultKS.elim_i) <0.05]
  # allGO_i <- genesInTerm(GOdata_i)
  # x_all_i <- allGO_i[names(score_elim_i)]
  #
  # Term_i <- as.data.frame(do.call(rbind,lapply(1:length(x_all_i),function(j){
  #   Term_i <- Term(xx_i <- xx[names(xx) %in% names(x_all_i)[j]][[1]])
  #   Term_i <- data.frame(term_name=Term_i,genes= paste(unlist(x_all_i[[j]]),collapse = ","))
  #   return(Term_i)
  # })))
  message(i)
  x_s_i <-  gprofiler2::gost(query = unique(unlist(x_MF_list_ara[[i]])),
                             organism = "athaliana", ordered_query = FALSE,
                             multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                             measure_underrepresentation = FALSE, evcodes = T,
                             user_threshold = 0.05, correction_method = "g_SCS",
                             #domain_scope = "annotated",
                             custom_bg = unique(x_orth_sat2$ortholog_ensg),#unique(dat3_Sat$ci_gen),
                             numeric_ns = "", sources = "GO:MF", as_short_link = FALSE)

  if(!is.null(x_s_i)){
    Term_i <- data.frame(term=x_s_i$result$term_name,genes=x_s_i$result$intersection)

    x_i_go_arath  <- data.frame(feature=names(x_MF_list)[[i]],term_name=Term_i$term,term_id =x_s_i$result$term_id,
                          #names(x_all_i),
                          genes=Term_i$genes)
  } else {
    x_i_go_arath  <- data.frame(feature=names(x_MF_list)[[i]],term_name=NA,term_id =NA,
                          #names(x_all_i),
                          genes=NA)
  }

  #x_i_go_arath  <- data.frame(feature=names(x_MF_list)[[i]],term_name=Term_i$term_name,term_id = names(x_all_i), genes=Term_i$genes)
  return(x_i_go_arath)
})
x_i_go_arath <- do.call(rbind,x_i_go_arath)
x_i_go_arath2 <- x_i_go_arath[which(x_i_go_arath$term_name!="molecular_function"),]
write.csv(as.data.frame(x_i_go_arath2),paste0(dir,"/MF/TopGO_2_filters/OSJ_ARATH/","GOST_ARATH.csv"),na = "",row.names = F)

message("Genes used for enrichment for ARATH: ",length(unique(unlist(x_MF_list_ara))))

#message("unique MF for enrichment for OSJ: ",length(unique(unlist(x_MF_list_ara))))
message("Genes used for enrichment for ARATH: ",length(unique(unlist(x_MF_list_ara))))

###########################################################################
#get intersection of GO
xx_intersect <- as.data.frame(intersect(x_i_go2$term_name,x_i_go_arath2$term_name))
colnames(xx_intersect) <- "GO_MF"
xx_intersect$ncats_OSJ <- NA
xx_intersect$ngenes_OSJ <- NA
xx_intersect$cats_OSJ <- NA
xx_intersect$genes_OSJ <- NA
#
xx_intersect$ncats_ARATH <- NA
xx_intersect$ngenes_ARATH  <- NA
xx_intersect$cats_ARATH  <- NA
xx_intersect$genes_ARATH  <- NA

for(i in 1:nrow(xx_intersect)){
  xx_intersect$ncats_OSJ[[i]]  <- length(x_i_go2$feature[which(x_i_go2$term_name==xx_intersect$GO_MF[[i]])])
  xx_intersect$cats_OSJ[[i]]  <- paste(x_i_go2$feature[which(x_i_go2$term_name==xx_intersect$GO_MF[[i]])],collapse = "/")
  xx_intersect$ngenes_OSJ[[i]] <- length(unique(unlist(strsplit(x_i_go2$genes[which(x_i_go2$term_name==xx_intersect$GO_MF[[i]])],","))))
  xx_intersect$genes_OSJ[[i]] <- paste(unique(unlist(strsplit(x_i_go2$genes[which(x_i_go2$term_name==xx_intersect$GO_MF[[i]])],","))),collapse = "/")

  ###

  xx_intersect$ncats_ARATH[[i]]  <- length(x_i_go_arath2$feature[which(x_i_go_arath2$term_name==xx_intersect$GO_MF[[i]])])
  xx_intersect$cats_ARATH[[i]]  <- paste(x_i_go_arath2$feature[which(x_i_go_arath2$term_name==xx_intersect$GO_MF[[i]])],collapse = "/")
  xx_intersect$ngenes_ARATH[[i]] <- length(unique(unlist(strsplit(x_i_go_arath2$genes[which(x_i_go_arath2$term_name==xx_intersect$GO_MF[[i]])],","))))
  xx_intersect$genes_ARATH[[i]] <- paste(unique(unlist(strsplit(x_i_go_arath2$genes[which(x_i_go_arath2$term_name==xx_intersect$GO_MF[[i]])],","))),collapse = "/")

  };rm(i)

xx_intersect <- xx_intersect[order(xx_intersect$ncats_OSJ,decreasing = T),]
row.names(xx_intersect) <- xx_intersect$GO_MF
write.csv(as.data.frame(xx_intersect),paste0(dir,"/MF/TopGO_2_filters/","intersect.csv"),na = "",row.names = F)
###########################################################################
#MOST FREQUENT GENES

x_genes <- unlist(strsplit(xx_intersect$genes_OSJ,"/"))
x_genes_genes <- tapply(x_genes,x_genes,length)
x_genes_df <- data.frame(gen=names(x_genes_genes),count=x_genes_genes)
x_genes_df <- x_genes_df[order(x_genes_df$count,decreasing = T),]
write.csv(as.data.frame(x_genes_df),paste0(dir,"/MF/TopGO_2_filters/","OSJ_MF_intersected_gene_counts.csv"),na = "",row.names = F)
###########################################################################
x_genes_arat <- unlist(strsplit(xx_intersect$genes_ARATH,"/"))
x_genes_arat <- tapply(x_genes_arat,x_genes_arat,length)
x_genes_df_arat <- data.frame(gen=names(x_genes_arat),count=x_genes_arat)
x_genes_df_arat <- x_genes_df_arat[order(x_genes_df_arat$count,decreasing =  T),]
write.csv(as.data.frame(x_genes_df_arat),paste0(dir,"/MF/TopGO_2_filters/OSJ_ARATH/","ARATH_MF_intersected__gene_counts.csv"),na = "",row.names = F)
###########################################################################
#require(FactoMineR)
Pcc_int <- corrplot::corrplot(cor(scale(xx_intersect[,-c(1,4,5,8,9)],center = T,scale = T)),method = "shade",type = "lower",diag = F,hclust.method = "ward.D")
plot(xx_intersect$ncats_OSJ,xx_intersect$ncats_ARATH,pch=19, xlab="OSJ functional groups counts per MF", ylab="ARATH functional groups counts per MF")
cor.test(xx_intersect$ncats_OSJ,xx_intersect$ncats_ARATH)
#plot(xx_intersect)
plot(xx_intersect$ngenes_OSJ,xx_intersect$ngenes_ARATH,pch=19,xlab="OSJ genes counts per MF", ylab="ARATH genes counts per MF")
cor.test(xx_intersect$ngenes_OSJ,xx_intersect$ngenes_ARATH)
###########################################################################
#prop test mf and gene counts
ngenes_insect <- xx_intersect[,-c(1,2,4,5,6,8,9)]
ns_int =apply(ngenes_insect,2,sum)

p_values_int=array(NA,dim(ngenes_insect)[1])

for(i in 1:dim(ngenes_insect)[1]){
  xs=as.numeric(ngenes_insect[i,1:2])
  test=prop.test(x=xs,n =ns_int)
  p_values_int[i]=test$p.value
}



table(p_values_int>0.05)/length(p_values_int)
x_df_p_int <- data.frame(BP=row.names(ngenes_insect),P=p_values_int,fdr=p.adjust(p_values_int,method = "fdr"))
x_df_p_int <- x_df_p_int[which(x_df_p_int$fdr<0.05),]
write.csv(as.data.frame(x_df_p_int),paste0(dir,"/MF/TopGO_2_filters/","prop_test_gene_MF.csv"),na = "",row.names = F)

###########################################################################
###########################################################################
#get gene counts per category
###########################################################################
#
unique_MF_to_genes <- names(x_MF_list_ara)
x_counts_cat <- list()
for(i in 1:length(unique_MF_to_genes)){

  # i_osj <- x_MF_list[[i]]
  # i_arath <- x_MF_list_ara[[i]]
  # #j <- 1
  # x_sub_j <- list()
  # for(j in 1:length(x_MF_list_ara[[i]])){
  #   #message(j)
  #   x_sub_j[[j]] <- x_orth_sat2$[which( x_orth_sat2$ortholog_ensg == x_MF_list_ara[[i]][[j]])]
  #
  # };rm(j)
  # #
  x_counts_cat[[i]] <- data.frame(CAT=unique_MF_to_genes[[i]],
             OSJ=length(unique(unlist(x_MF_list[[i]]))),
             ARATH=length(unique(unlist(x_MF_list_ara[[i]])))
             )
};rm(i)
x_counts_cat <- do.call(rbind,x_counts_cat)
write.csv(as.data.frame(x_counts_cat),paste0(dir,"/MF/TopGO_2_filters/","Genes_counts.csv"),na = "",row.names = F)
###########################################################################
ns_cat=apply(x_counts_cat[,-1],2,sum)

p_values=array(NA,dim(x_counts_cat)[1])

for(i in 1:dim(x_counts_cat)[1]){
  xs=as.numeric(x_counts_cat[i,2:3])
  test=prop.test(x=xs,n = ns_cat)
  p_values[i]=test$p.value
}



table(p_values>0.05)/length(p_values)
x_df_p_cat <- data.frame(BP=x_counts_cat$CAT,P=p_values,fdr=p.adjust(p_values,method = "fdr"))
x_df_p_cat <- x_df_p_cat[which(x_df_p_cat$fdr<0.05),]
write.csv(as.data.frame(x_df_p_cat),paste0(dir,"/MF/TopGO_2_filters/","prop_test_gene_cat.csv"),na = "",row.names = F)


###########################################################################
#get ortholog counts per OSJ genes
###############
un_OSJ <- unique(x_orth_sat2$input_ensg)
x_count_genes <- lapply(1:length(un_OSJ), function(i){
#3i <- 1
  #message(i)
  x_cat_j <- lapply(1:length(x_MF_list), function(j){
    sum_x_MF <- sum(unlist(x_MF_list[[j]]) %in% un_OSJ[[i]])
    if(sum_x_MF>0){
      x_cat_j  <- names(x_MF_list)[[j]]
      return(x_cat_j)
    }
      })
  x_cat_j<-x_cat_j[!sapply(x_cat_j,is.null)]
  x_data <- data.frame(gene=un_OSJ[[i]],
                    orth_count=length(x_orth_sat2$ortholog_ensg[which(x_orth_sat2$input==un_OSJ[[i]])]),
                    ncats=length(x_cat_j),
                    cats=paste(x_cat_j,collapse = "/")
  )
  return(x_data)
})
x_count_genes <- do.call(rbind, x_count_genes)
x_count_genes <- x_count_genes[order(x_count_genes$orth_count,decreasing = T),]
write.csv(as.data.frame(x_count_genes),paste0(dir,"/MF/TopGO_2_filters/","orth_counts.csv"),na = "",row.names = F)

#############################################################
#graph one species
x_ara <- graphGOspecies(df=x_i_go_arath2,#x_s_ara$result,
                    GOterm_field="term_name",
                    option = "GO",
                    numCores=6,
                    saveGraph=T,
                    outdir = paste0(dir,"/MF/TopGO_2_filters/OSJ_ARATH"))



x2_area <- graphGOspecies(df=x_i_go_arath2,#x_s_ara$result,
                     GOterm_field="term_name",
                     option = "Categories",
                     numCores=6,
                     saveGraph=T,
                     outdir = paste0(dir,"/MF/TopGO_2_filters/OSJ_ARATH"))


#Get nodes with values greater than 95%
#perc_Ara <- x_ara$nodes[which(x_ara$nodes$GO_WEIGHT > quantile(x_ara$nodes$GO_WEIGHT,probs = 0.9)),]
perc_Ara <- x_ara$nodes[which(x_ara$nodes$GO_WEIGHT >= quantile(x_ara$nodes$GO_WEIGHT,probs = 0.9)),]
write.csv(perc_Ara,paste0(dir,"/MF/TopGO_2_filters/OSJ_ARATH/","PERC_ARATH.csv"),na = "",row.names = F)
# visualize nodes filtered

###########################
#comparing GO terms between OSJ and ARATH

x_COMP <- compareGOspecies(df1 = x_i_go2,#x_s$result,
                           df2 = x_i_go_arath2,#x_s_ara$result,
                           GOterm_field = "term_name",
                           species1="O. sativa (japonica)",
                           species2="A. thaliana")
plot(hclust(x_COMP$distance,method = "ward.D"))
x_COMP$graphics

#nrow(x_COMP$shared_GO_list)
stats_count_un_MF_Status <- data.frame(
  SP1 = nrow(x_COMP$unique_GO_list[which(x_COMP$unique_GO_list$species=="O. sativa (japonica)"),]),
  SP2 = nrow(x_COMP$unique_GO_list[which(x_COMP$unique_GO_list$species=="A. thaliana"),]),
  SHARED = nrow(x_COMP$shared_GO_list[which(x_COMP$shared_GO_list$species=="Shared"),]))
write.csv(stats_count_un_MF_Status,paste0(dir,"/MF/TopGO_2_filters/UNIQUE_SHARED_BALANCE.csv"),na = "",row.names = F)

#length(x_COMP$unique_GO_list$species=="A. thaliana")

library(immunarch)
xx_inm <- immunarch::immunr_hclust(as.matrix(x_COMP$distance), .dist = T,.method = "ward.D")
res <- hcut(x_COMP$distance, k = 16, stand = F,hc_method = "ward.D",isdiss=T ) #k=5

png(
  paste0(dir,"/MF/TopGO_2_filters/JACCARD_CLUSTER.png"),
  width     = 45,
  height    = 29,
  units     = "in",
  res       = 400,
  pointsize = 5
)

fviz_dend(res, rect = TRUE,cex = 2.5,
          #type="phylogenic",
          type="rectangle",
          horiz=T,
          lwd=2.5,
repel = TRUE,color_labels_by_k = F,labels_track_height=5,rect_border = F,main ="")
dev.off()



pheatmap(as.matrix(x_COMP$distance),clustering_method="ward.D",fontsize =10,fontsize_col = 8,angle_col=315,fontsize_row = 8,cluster_rows = F)

x_graph_cat <- graph_two_GOspecies(x=x_COMP,
                               species1="O. sativa (japonica)",
                               species2="A. thaliana",
                               GOterm_field="term_name",
                               numCores=6,
                               saveGraph = T,
                               option= "Categories",
                               outdir = paste0(dir,"/MF/TopGO_2_filters"))

x_graph_go <- graph_two_GOspecies(x=x_COMP,
                                   species1="O. sativa (japonica)",
                                   species2="A. thaliana",
                                   GOterm_field="term_name",
                                   numCores=2,
                                   saveGraph = T,
                                   option= "GO",
                                   outdir = paste0(dir,"/MF/TopGO_2_filters"))

perc_twosp <- x_graph_go$nodes[which(x_graph_go$nodes$GO_WEIGHT >= quantile(x_graph_go$nodes$GO_WEIGHT,probs = 0.9)),]
write.csv(perc_twosp,paste0(dir,"/MF/TopGO_2_filters/","PERC_two_spp.csv"),na = "",row.names = F)
############################################################################################
#get shared term name among percentile 90th
allperc <- c(perc_twosp$GO,perc$GO,perc_Ara$GO)# [duplicated(]
allperc[duplicated(allperc)]
write.csv(allperc[duplicated(allperc)],paste0(dir,"/MF/TopGO_2_filters/","PERC_duplicated_all_spp.csv"),na = "",row.names = F)

# visualize nodes filtered

tapply(x_COMP$unique_GO_list$species,x_COMP$unique_GO_list$species,length)
tapply(x_COMP$shared_GO_list$species,x_COMP$shared_GO_list$species,length)
#View(x_COMP$unique_GO_list)

x_eval <- evaluateGO_species(df1 = x_i_go2,df2 = x_i_go_arath2,#x_s$result, x_s_ara$result,
                             species1="O. sativa (japonica)", species2="A. thaliana",GOterm_field = "term_name",test = "prop")

x_eval_cat <- evaluateCAT_species(x_i_go2,x_i_go_arath2,#x_s$result, x_s_ara$result,
                             species1="O. sativa (japonica)", species2="A. thaliana",GOterm_field = "term_name",test="prop")

x_eval <- x_eval[which(x_eval$FDR <0.05),]
x_eval_cat <- x_eval_cat[which(x_eval_cat$FDR <0.05),]
write.csv(x_eval,paste0(dir,"/MF/TopGO_2_filters/","EVAL_MF.csv"),na = "",row.names = F)
write.csv(x_eval_cat,paste0(dir,"/MF/TopGO_2_filters/","EVAL_GROUPS.csv"),na = "",row.names = F)
#################################################################################

#get graphs for genes and species (OSJ)
group_BP  <- unique(TERM_GO_list$TERM_FILTERED)


x_gene_Graph <- lapply(1:length(group_BP),function(i){
  #message(i)
  x_GO_sub <- TERM_GO_list$genes[which(TERM_GO_list$TERM_FILTERED==group_BP[[i]])]
  if(length(x_GO_sub)==1){
    x <-
      data.frame(x1=x_GO_sub,X2=x_GO_sub,GO=group_BP[[i]])
  } else {
    x <-
    data.frame(t(utils::combn(x_GO_sub, 2)),GO=group_BP[i])
  }
  colnames(x) <- c("source","target","GO")
return(x)
})


x_gene_Graph <- do.call(rbind,x_gene_Graph)
#SAVING GRAPH AS UNDIRECTED IN IGRAPH
x_gene_Graph <- igraph::graph_from_data_frame(d = x_gene_Graph, directed = FALSE)
#V(x1)$name

igraph::write.graph(x_gene_Graph,
                    file= paste0(dir,"/MF/TopGO_2_filters/","genes_OSJ.graphml"),
                    format = "graphml")
#################################################################################
#get graphs for genes and species (ARATH)


x_gene_Graph_ARATH <- lapply(1:length(group_BP),function(i){
  #message(i)
  x_GO_sub <- TERM_GO_list_ARA$genes[which(TERM_GO_list_ARA$TERM_FILTERED==group_BP[[i]])]
  if(length(x_GO_sub)==1){
    x <-
      data.frame(x1=x_GO_sub,X2=x_GO_sub,GO=group_BP[[i]])
  } else {
    x <-
      data.frame(t(utils::combn(x_GO_sub, 2)),GO=group_BP[i])
  }
  colnames(x) <- c("source","target","GO")
  return(x)
})


x_gene_Graph_ARATH <- do.call(rbind,x_gene_Graph_ARATH)
#SAVING GRAPH AS UNDIRECTED IN IGRAPH
x_gene_Graph_ARATH <- igraph::graph_from_data_frame(d = x_gene_Graph_ARATH, directed = FALSE)
#V(x1)$name

igraph::write.graph(x_gene_Graph_ARATH,
                    file= paste0(dir,"/MF/TopGO_2_filters/OSJ_ARATH/","genes_ARATH.graphml"),
                    format = "graphml")

#################################################################################
save.image(paste0(dir,"/","STEP1"))
#load(paste0(dir,"/","STEP1"))


#################################################################################
