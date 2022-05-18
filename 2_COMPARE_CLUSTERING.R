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
devtools::load_all()
#require(GOCompare)
require(GOfuncR)
require(agricolae)
#require(GO.db)
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
require(ggpubr)

dir <- "D:/TESIS_PHD/CHAPTER1"
#############################################################################################
load(paste0(dir,"/","STEP1"))
#############################################################################################
# compare go counts per groups
OSJ_tab <- as.data.frame(tapply(x_i_go2$feature,x_i_go2$feature,length))#x_s$result$feature,x_s$result$feature,length))
OSJ_tab$GROUPS <- row.names(OSJ_tab)
colnames(OSJ_tab) <- c("Freq","Group")
OSJ_tab <- OSJ_tab[,c("Group","Freq")]
row.names(OSJ_tab) <- NULL

####

ARATH_tab <- as.data.frame(tapply(x_i_go_arath2$feature,x_i_go_arath2$feature,length))#x_s_ara$result$feature,x_s_ara$result$feature,length))
ARATH_tab$GROUPS <- row.names(ARATH_tab)
colnames(ARATH_tab) <- c("Freq","Group")
ARATH_tab <- ARATH_tab[,c("Group","Freq")]
row.names(ARATH_tab) <- NULL

UNIQUE_TERMS_GROUPS <- unique(TERM_GO_list$TERM_FILTERED)#[-3]

tab_for_zscores <- list()
for(i in 1:length(UNIQUE_TERMS_GROUPS)){
  message(i)
  sp1 = OSJ_tab$Freq[which(OSJ_tab$Group==UNIQUE_TERMS_GROUPS[[i]])]
  sp2 = ARATH_tab$Freq[which(ARATH_tab$Group==UNIQUE_TERMS_GROUPS[[i]])]
  if(length(sp1)==0){
    sp1 <- 0
  }

  if(length(sp2)==0){
    sp2 <- 0
  }
  tab_for_zscores[[i]] <- data.frame(Group=UNIQUE_TERMS_GROUPS[[i]],
                                     sp1=sp1,
                                     sp2=sp2)

}
tab_for_zscores <- do.call(rbind,tab_for_zscores)
write.csv(tab_for_zscores,paste0(dir,"/MF/TopGO_2_filters/","STATS_GO.csv"),na = "",row.names = F)

wilcox.test(x=tab_for_zscores$sp1,
            y=tab_for_zscores$sp2,paired = T,
            alternative="two.sided",
            estimate=T,
            exact=T,
            correct=T)

hist(tab_for_zscores$sp1)
hist(tab_for_zscores$sp2)

cor.test(tab_for_zscores$sp1,tab_for_zscores$sp2,method = "spearman")

################################################################################
############################################################
# compare go counts per groups uniques per spp


spp <- unique(x_COMP$unique_GO_list$species)

x_GO_list <- list()
for(i in 1:2){
  x_GO_un  <- as.data.frame(table(x_COMP$unique_GO_list$feature[which(x_COMP$unique_GO_list$species==spp[i])],
                                  x_COMP$unique_GO_list$GO[which(x_COMP$unique_GO_list$species==spp[i])]))
  x_GO_un <- x_GO_un[which(x_GO_un$Freq!=0),]
  x_GO_un$spp <- spp[i]
  x_GO_un <- x_GO_un[,c("Var1","Var2","spp")]
  colnames(x_GO_un) <- c("GROUP","MF","SPP")
  x_GO_list[[i]] <- x_GO_un
};rm(i)

x_GO_list <- do.call(rbind,x_GO_list)

unique_groups_un <- unique(x_GO_list$GROUP)
un_un_list <- list()
for(i in 1:2){
  x_j <- list()
  for(j in 1:length(unique_groups_un)){
    x_j[[j]] <- data.frame(#GROUP=unique_groups_un[[j]],
      Freq=
        length(x_GO_list$MF[which(x_GO_list$GROUP==unique_groups_un[[j]] & x_GO_list$SPP==spp[i])])
      #spp=spp[i]
    )
  };rm(j)
  x_j <- do.call(rbind,x_j)
  un_un_list[[i]] <- x_j
  rm(x_j)

};rm(i)

un_un_list <- do.call(cbind,un_un_list)
un_un_list$GROUP <- unique_groups_un
colnames(un_un_list) <- c("sp1","sp2","GROUP")
un_un_list <- un_un_list[,c("GROUP","sp1","sp2")]


write.csv(un_un_list,paste0(dir,"/MF/TopGO_2_filters/","STATS_UNIQUE_GO.csv"),na = "",row.names = F)


wilcox.test(x=un_un_list$sp1,
            y=un_un_list$sp2,paired = T,
            alternative="two.sided",
            estimate=T,
            exact=T,
            correct=T)

hist(un_un_list$sp1)
hist(un_un_list$sp2)

plot(un_un_list$sp1,un_un_list$sp2)
cor.test(un_un_list$sp1,un_un_list$sp2,method = "spearman")
#x_eval$BH <-  stats::p.adjust(x_eval$pvalue, method = "bonferroni")
###################################################################################


############################################################
# compare go counts per groups for shared
unique_groups_shared <- unique(x_COMP$shared_GO_list$feature)
shared_list <- list()

for(j in 1:length(unique_groups_shared)){
  shared_list[[j]] <- data.frame(#GROUP=unique_groups_un[[j]],
    Freq=
      length(x_COMP$shared_GO_list$feature[
        which(x_COMP$shared_GO_list$feature==unique_groups_shared[[j]])]
        #spp=spp[i]
      ))
};rm(j)


shared_list <- do.call(rbind,shared_list)
shared_list$GROUP <- unique_groups_shared
shared_list <- shared_list[,c("GROUP","Freq")]
hist(shared_list$Freq)


write.csv(shared_list,paste0(dir,"/MF/TopGO_2_filters/","STATS_SHARED_GO.csv"),na = "",row.names = F)
#################################################################################
#bind stats
tab_for_zscores2 <- tab_for_zscores
shared_list2 <- shared_list
un_un_list2 <- un_un_list
colnames(tab_for_zscores2)[1] <- "GROUP"
x_stats <- dplyr::full_join(tab_for_zscores2,shared_list,by="GROUP")
x_stats <- dplyr::full_join(x_stats,un_un_list2,by="GROUP")
colnames(x_stats) <- c("BP","OSJ","ARATH","BOTH","OSJ_UNIQUE","ARATH_UNIQUE")
x_stats <- x_stats[order(x_stats$BOTH,decreasing = T),]
write.csv(x_stats,paste0(dir,"/MF/TopGO_2_filters/","STATS_COUNTS.csv"),na = "",row.names = F)
#################################################################################
#i <- 1
un_MF_enriched <-  unique(c(x_i_go2$term_name,x_i_go_arath2$term_name))
xx_mf_en <- lapply(1:length(un_MF_enriched),function(i){
  xx_mf_en <- data.frame(
    MF=un_MF_enriched[[i]],
    sp1=length(x_i_go2$feature[x_i_go2$term_name %in% un_MF_enriched[[i]]]),
    sp2=length(x_i_go_arath2$feature[x_i_go_arath2$term_name %in% un_MF_enriched[[i]]]))
  return(xx_mf_en)
})

xx_mf_en <- do.call(rbind,xx_mf_en)
#############################

ns=apply(xx_mf_en[,-1],2,sum)

p_values=array(NA,dim(xx_mf_en)[1])

for(i in 1:dim(xx_mf_en)[1]){
  xs=as.numeric(xx_mf_en[i,2:3])
  test=prop.test(x=xs,n = ns[1:2])
  p_values[i]=test$p.value
}


table(p_values>0.05)/length(xx_mf_en$MF)
x_df_p <- data.frame(MF=xx_mf_en$MF,P=p_values,fdr=p.adjust(p_values,method = "fdr"))
x_df_p <- x_df_p[which(x_df_p$fdr <0.05),]
x_df_p <- x_df_p[order(x_df_p$P,decreasing = F),]
#write.csv(x_df_p,paste0(dir,"/MF/TopGO_2_filters/","prop_MF.csv"),na = "",row.names = F)
#}
#################################################################################
#x_stats_prop=x_stats[,-1]/apply(x_stats[,-1],2,sum)
#row.names(x_stats_prop) <- x_stats[,1]

x_stats2 <- x_stats[,c(1:3)]
#ns=apply(x_stats2[,-1],2,sum)
ns=apply(x_stats2[,-1],2,sum)
p_values=array(NA,dim(x_stats2)[1])

for(i in 1:dim(x_stats2)[1]){
  xs=as.numeric(x_stats2[i,2:3])
  test=prop.test(x=xs,n =ns)
  p_values[i]=test$p.value

  rm(xs,test)
};rm(i)



table(p_values>0.05)/nrow(x_stats2)
x_df_p2 <- data.frame(BP=x_stats2$BP,P=p_values,fdr=p.adjust(p_values,method = "fdr"))
x_df_p2 <- x_df_p2[which(x_df_p2$fdr<0.05),]
#write.csv(x_df_p2,paste0(dir,"/MF/TopGO_2_filters/","prop_CAT.csv"),na = "",row.names = F)
#################################################################################
#prop unique
x_stats3 <- x_stats[,c(1,5,6)]
#ns=apply(x_stats2[,-1],2,sum)

p_values=array(NA,dim(x_stats3)[1])

for(i in 1:dim(x_stats3)[1]){
  xs=as.numeric(x_stats3[i,2:3])
  test=prop.test(x=xs,n = ns)
  p_values[i]=test$p.value
}



table(p_values>0.05)/21
x_df_p3 <- data.frame(BP=x_stats3$BP,P=p_values,fdr=p.adjust(p_values,method = "fdr"))
x_df_p3 <- x_df_p[which(x_df_p3$fdr<0.05),]
#write.csv(x_df_p3,paste0(dir,"/MF/TopGO_2_filters/","prop_CAT_unique.csv"),na = "",row.names = F)
#################################################################################
# treemaps
#Plot them

png(
  paste0(dir,"/MF/TopGO_2_filters/OSJ_TREE.png"),
  width     = 3.25,
  height    = 3.25,
  units     = "in",
  res       = 1200,
  pointsize = 4
)

treemap(un_un_list[,c(1,2)],
                    palette = "Set1",                        # Select your color palette from the RColorBrewer presets or make your own.
                    title="Unique molecular functions per biological process (O. sativa subsp japonica)",                      # Customize your title
                    fontsize.title=6,
                    index="GROUP",
                    vSize="sp1",
                    type="index"
)


dev.off()

png(
  paste0(dir,"/MF/TopGO_2_filters/ARATH_TREE.png"),
  width     = 3.25,
  height    = 3.25,
  units     = "in",
  res       = 1200,
  pointsize = 4
)
treemap(un_un_list[,c(1,3)],
                       palette = "Set1",                        # Select your color palette from the RColorBrewer presets or make your own.
                       title="Unique molecular functions per biological process (A. thaliana)",                      # Customize your title
                       fontsize.title=6,
                       index="GROUP",
                       vSize="sp2",
                       type="index"
)


dev.off()
png(
  paste0(dir,"/MF/TopGO_2_filters/BOTH_TREE.png"),
  width     = 3.25,
  height    = 3.25,
  units     = "in",
  res       = 1200,
  pointsize = 4
)

treemap(shared_list,
        palette = "Set1",                        # Select your color palette from the RColorBrewer presets or make your own.
        title="Molecular functions per biological process for O.sativa and A. thaliana",                      # Customize your title
        fontsize.title=6,
        index="GROUP",
        vSize="Freq",
        type="index"
)
dev.off()



##################################################################################
#getting node weights for clustering

node_cat <- unique(x2_area$nodes$feature)

node_Cat_list <- list()
for(i in 1:length(node_cat)){
  #  message(i)

  group= node_cat[[i]]
  shared=x_graph_cat$nodes$GO_WEIGHT[which(x_graph_cat$nodes$GO==node_cat[[i]])]
  sp1_f=x2$nodes$WEIGHT[which(x2$nodes$feature==node_cat[[i]])]
  sp2_f=x2_area$nodes$WEIGHT[which(x2_area$nodes$feature==node_cat[[i]])]

  if(length(shared)==0){
    shared <- NA
  }


  if(length(sp1_f)==0){
    sp1_f <- NA
  }

  if(length(sp2_f)==0){
    sp2_f <- NA
  }

  node_Cat_list[[i]] <- data.frame(group= node_cat[[i]],shared=shared,sp1=sp1_f,sp2=sp2_f)
  rm(group,sp1_f,sp2_f)
};rm(i)

node_Cat_list <- do.call(rbind,node_Cat_list)

node_Cat_list2 <-node_Cat_list
node_Cat_list2 <- node_Cat_list2[complete.cases(node_Cat_list2),]
node_Cat_list2[,c(2:4)] <- scale(node_Cat_list2[,-c(1)],center = T,scale=T)
write.csv(node_Cat_list,paste0(dir,"/MF/TopGO_2_filters/","node_weights.csv"),na = "",quote = F,row.names = F)

library("plot3D")
text3D(node_Cat_list2$shared, node_Cat_list2$sp1, node_Cat_list2$sp2, clab = c("Biological process (Scaled node weights)"),
       pch = 19, cex = 0.8,bty = "g",col = ramp.col(c("blue", "yellow", "red")),theta = 45, phi = 20,
       xlab = "Shared",
       ylab = "OSJ", zlab =  "ARATH",
       colvar = rowMeans(node_Cat_list2[,-c(1)]),
       labels = node_Cat_list2$group, adj = 0.1, font = 2,
       ticktype = "detailed", d = 1
)

row.names(node_Cat_list2) <- node_Cat_list2[,c(1)]
plot(hclust(dist(node_Cat_list2[,-1])))

require(agricolae)
melt_to_gg <- reshape2::melt(node_Cat_list2,"group")
x_kruskal_nodes <- kruskal(melt_to_gg$value,melt_to_gg$variable,group = F)
corrplot::corrplot(cor(node_Cat_list2[,-1]),method = "shade",type = "lower",diag = F,order = 'hclust')
t.test(node_Cat_list2$sp1, node_Cat_list2$sp2)

ggplot(data = melt_to_gg, mapping = aes(x = variable, y = value,fill=variable)) +
  geom_boxplot(alpha = 0) +
  geom_jitter(alpha = 0.3, color = "black")

fviz_nbclust(node_Cat_list2[,-1], hcut, method = "silhouette",dist="euclidean",kmax=15) +
  #geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "silhouette method")

row.names(node_Cat_list2) <- node_Cat_list2$group
res <- hcut(node_Cat_list2[,-1], k = 2, stand = TRUE)
png(
  paste0(dir,"/MF/TopGO_2_filters/NODE_CAT_CLUTER.png"),
  width     = 20,
  height    = 10,
  units     = "in",
  res       = 600,
  pointsize = 4
)

fviz_dend(res, rect = TRUE,cex = 0.7,
          #type="phylogenic",
          type="rectangle",
          horiz=T,
          repel = TRUE)
dev.off()
fviz_silhouette(res)
#################################################################################
#getting node weights for MF

node_go <- unique(c(x_i_go2$term_name,x_i_go_arath2$term_name))#x_s$result$term_name,x_s_ara$result$term_name))

node_go_list <- list()
for(i in 1:length(node_go)){
  #  message(i)

  group= node_go[[i]]
  shared=x_graph_go$nodes$GO_WEIGHT[which(x_graph_go$nodes$GO==node_go[[i]])]
  sp1_f=x$nodes$GO_WEIGHT[which(x$nodes$GO==node_go[[i]])]
  sp2_f=x_ara$nodes$GO_WEIGHT[which(x_ara$nodes$GO==node_go[[i]])]

  if(length(shared)==0){
    shared <- NA
  }


  if(length(sp1_f)==0){
    sp1_f <- NA
  }

  if(length(sp2_f)==0){
    sp2_f <- NA
  }

  node_go_list [[i]] <- data.frame(group= node_go[[i]],shared=shared,sp1=sp1_f,sp2=sp2_f)
  rm(group,sp1_f,sp2_f)
};rm(i)

node_go_list <- do.call(rbind,node_go_list)

node_go_list2 <-node_go_list
node_go_list2 <- node_go_list2[complete.cases(node_go_list2),]
node_go_list2[,c(2:4)] <- scale(node_go_list2[,-c(1)],center = T,scale=T)

library("plot3D")
text3D(node_go_list2$shared, node_go_list2$sp1, node_go_list2$sp2, clab = c("Metabolic pathway (Scaled node weights)"),
       pch = 19, cex = 0.8,bty = "g",col = ramp.col(c("blue", "yellow", "red")),theta = 45, phi = 20,
       xlab = "Shared",
       ylab = "OSJ", zlab =  "ARATH",
       colvar = rowMeans(node_go_list2[,-c(1)]),
       labels = node_go_list2$group, adj = 0.1, font = 2,
       ticktype = "detailed", d = 1
)


melt_to_gg2 <- reshape2::melt(node_go_list2,"group")
x_kruskal_nodes2 <- kruskal(melt_to_gg2$value,melt_to_gg2$variable,group = F)

ggplot(data = melt_to_gg2, mapping = aes(x = variable, y = value,fill=variable)) +
  geom_boxplot(alpha = 0) +
  geom_jitter(alpha = 0.3, color = "black")



t.test(node_go_list2$sp1, node_go_list2$sp2)

corrplot::corrplot(cor(node_go_list2[,-1]),method = "shade",type = "lower",diag = F,order = 'hclust')

#plot(1)
fviz_nbclust(node_go_list2[,-1], hcut, method = "silhouette",dist="euclidean",kmax=50) +
  #geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "silhouette")

row.names(node_go_list2) <- node_go_list2$group
res_go <- hcut(node_go_list2[,-1], k = 10, stand = TRUE)#k=3


png(
  paste0(dir,"/MF/TopGO_2_filters/NODE_MF_CLUTER.png"),
  width     = 23000,
  height    = 15000,
  units     = "px",
  res       = 600,
  pointsize = 8
)
fviz_dend(res_go, rect = TRUE,cex = 1.4,horiz =T,
          #type="phylogenic",
          type="rectangle",
          repel = TRUE)
dev.off()


fviz_silhouette(res_go)


################
#3d plot
# text3D(node_go_list2$shared, node_go_list2$sp1, node_go_list2$sp2, clab = c("Metabolic pathways clusters"),
#        pch = 19, cex = 0.9,bty = "g",#col = ramp.col(c("blue", "yellow", "red")),
#        theta = 45, phi = 20,
#        xlab = "Shared",
#        ylab = "OSJ", zlab =  "ARATH",
#        colvar = as.numeric(node_go_list2$cluster),
#        labels = node_go_list2$group, adj = 0.1, font = 2,
#        ticktype = "detailed", d = 10
# )
#################################################################################


png(
  paste0(dir,"/MF/TopGO_2_filters/MF_CLUTER_PCA.png"),
  width     = 21000,
  height    = 12000,
  units     = "px",
  res       = 600,
  pointsize = 4
)

fviz_cluster(res_go, data = node_go_list2[,-1], #palette = c("#FC4E07", "#00AFBB", "#E7B800"),
             ellipse.type = "convex",
             star.plot = F,
             labelsize = 12,
             repel = TRUE,

             ggtheme = theme_minimal() )
dev.off()

####################################
save.image("D:/TESIS_PHD/CHAPTER1/MF/TopGO_2_filters/STEP_3")
#load("D:/TESIS_PHD/CHAPTER1/MF/TopGO_2_filters/STEP_3")
