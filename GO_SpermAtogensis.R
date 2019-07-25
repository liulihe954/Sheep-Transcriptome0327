                   ##################################################################################################
                   ###                                   Set Up  and  Data prep                    #################
                   ################################################################################################
#setwd("/Users/liulihe95/Desktop/GO_data")
library(WGCNA);library(ppcor);library(igraph);require(cowplot);library(ggplot2);library(ggpubr)
library(ggridges)
library(extrafont)
library(dplyr)
library(plotly)
library(geomnet)
library(writexl)

# sessionInfo()
# theme_set(theme_ridges())
## library("corrplot");library("qgraph")
setwd('/Users/liulihe95/Desktop/Network_GO_0327');getwd()
options(stringsAsFactors = FALSE)
networkData_sperm = read.table("Spermatogenesis_GO.txt");
# reroder
a = c(2,4,5,6,9,13,16,18,19); b = c(2:19)[!c(2:19) %in% a]
networkData_sperm = networkData_sperm[,c(1,a,b,20:22)]
table(networkData_sperm$Significant) # 13 significant and 132 non-signif
datExpr_sperm <- as.data.frame(t(networkData_sperm[,c(2:19)]));names(datExpr_sperm) = networkData_sperm$Gene; rownames(datExpr_sperm) = names(networkData_sperm)[c(2:19)]
dim(networkData_sperm)
datExprGO <- datExpr_sperm[c(1:9),];datExprGT <- datExpr_sperm[c(10:18),]
     ###           <datExprGO>    data set ready to use       <datExprGT>     ####

######                    check for missing value                  ##########
gsg_all_sperm = goodSamplesGenes(datExpr_sperm, verbose = 3);
gsg_all_sperm$allOK
####   distance between samples / outliers
sample_Tree_GO = hclust(dist(datExprGO), method = "average")
sample_Tree_GT = hclust(dist(datExprGT), method = "average")
plot(sample_Tree_GO, main = "Sample clustering to detect outliers GO", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
plot(sample_Tree_GT, main = "Sample clustering to detect outliers GT", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
#####  NO NAs or missing values were detected. but distances may be susceptible (?)


                  ##################################################################################################
                 #####             1. Unweighted network analysis -  NC + ranking                            ######
                 ################################################################################################
### 1.Define calculation function - get_adjmat(dataset, rthreshold, pthreshold) --- get unweighted adj matrx --- use PearsonCor
require(ppcor)
# get_NC_cor = function(datExpr1,datExpr2,r_thres,p_thres){
  #calcu of matrx1
  cormatr1 <- cor(datExpr1)
  adjmatr1 = matrix(1,ncol(datExpr1),ncol(datExpr1))
  colnames(adjmatr1) = colnames(datExpr1)
  rownames(adjmatr1) = colnames(datExpr1)
  adjmatr1[abs(cormatr1) < r_thres] = 0
  #calcu of matrx2
  cormatr2 <- cor(datExpr2)
  adjmatr2 = matrix(1,ncol(datExpr2),ncol(datExpr2))
  colnames(adjmatr2) = colnames(datExpr2)
  rownames(adjmatr2) = colnames(datExpr2)
  adjmatr2[abs(cormatr2) < r_thres] = 0
  #use threshold and get pvalue matrix
  pvalmatr1 = matrix(0,ncol(datExpr1),ncol(datExpr1))
  pvalmatr2 = matrix(0,ncol(datExpr2),ncol(datExpr2))
  for(i in 1:(ncol(datExpr1)-1)){
    for (j in c((i+1):ncol(datExpr1))){
      r1 = cor.test(datExpr1[,i],datExpr1[,j])
      r2 = cor.test(datExpr2[,i],datExpr2[,j])
      pvalmatr1[i,j] = pvalmatr1[j,i] = r1$p.value
      pvalmatr2[i,j] = pvalmatr2[j,i] = r2$p.value
      if(r1$p.value >= p_thres){adjmatr1[i,j] = adjmatr1[j,i] = 0}
      if(r2$p.value >= p_thres){adjmatr2[i,j] = adjmatr2[j,i] = 0}
    }
  }
  #get all the basic NC
  NC1 = conformityBasedNetworkConcepts(adjmatr1)
  NC2 = conformityBasedNetworkConcepts(adjmatr2)
  #combine, rank and show
  basic_results = data.frame(density1 = NC1$fundamentalNCs$Density,
                             density2 = NC2$fundamentalNCs$Density,
                             centralization1 = NC1$fundamentalNCs$Centralization,
                             centralization2 = NC2$fundamentalNCs$Centralization,
                             heterogeneity1 = NC1$fundamentalNCs$Heterogeneity,
                             heterogeneity2 = NC2$fundamentalNCs$Heterogeneity)
  change_results = data.frame(#gene = colnames(datExpr1),
                              con_1 = NC1$fundamentalNCs$Connectivity,
                              scl_con_1 = NC1$fundamentalNCs$ScaledConnectivity,
                              con_2 = NC2$fundamentalNCs$Connectivity,
                              scl_con_2 = NC2$fundamentalNCs$ScaledConnectivity,
                              con_change = -(NC1$fundamentalNCs$Connectivity - NC2$fundamentalNCs$Connectivity),
                              scl_con_change = -(NC1$fundamentalNCs$ScaledConnectivity - NC2$fundamentalNCs$ScaledConnectivity),
                              rank_scl_con = rank(-abs(NC1$fundamentalNCs$ScaledConnectivity-NC2$fundamentalNCs$ScaledConnectivity)),
                              cls_coef_1 = NC1$fundamentalNCs$ClusterCoef,
                              cls_coef_2 = NC2$fundamentalNCs$ClusterCoef,
                              clst_coef_change = c(NC1$fundamentalNCs$ClusterCoef - NC2$fundamentalNCs$ClusterCoef),
                              rank_clstcoef = rank(-abs(NC1$fundamentalNCs$ClusterCoef-NC2$fundamentalNCs$ClusterCoef)))
  Results = list(NC1 = NC1, NC2 = NC2, cormatr_ref = cormatr1, cormatr_test = cormatr2, pvalue1 = pvalmatr1,pvalue2 = pvalmatr2 ,adjmatr_ref = adjmatr1, adjmatr_test = adjmatr2, basic = basic_results,change = change_results)
  return(Results)
}
# get_NC_pcor = function(datExpr1,datExpr2,r_thres, p_thres){
  cormatr1 <- cor(datExpr1)
  adjmatr1 = matrix(1,ncol(datExpr1),ncol(datExpr1))
  
  adjmatr1[abs(cormatr1) < r_thres] = 0
  #calcu of matrx2
  cormatr2 <- cor(datExpr2)
  adjmatr2 = matrix(1,ncol(datExpr2),ncol(datExpr2))
  
  adjmatr2[abs(cormatr2) < r_thres] = 0
  #use threshold
  for(i in 2:ncol(datExpr1)){
    r1 = cor.test(datExpr1[,i-1],datExpr1[,i])
    r2 = cor.test(datExpr2[,i-1],datExpr2[,i])
    if(r1$p.value >= p_thres){adjmatr1[i-1,i] = adjmatr1[i,i-1] = 0}
    if(r2$p.value >= p_thres){adjmatr2[i-1,i] = adjmatr2[i,i-1] = 0}
  }
  
  #calcu of matrx1
  cormatr1 <- pcor(datExpr1)$estimate
  adjmatr1 = matrix(1,ncol(datExpr1),ncol(datExpr1))
  colnames(adjmatr1) = colnames(datExpr1)
  rownames(adjmatr1) = colnames(datExpr1)
  adjmatr1[abs(cormatr1) < r_thres] = 0
  # calcu of matrx2
  cormatr2 <- pcor(datExpr2)$estimate
  adjmatr2 = matrix(1,ncol(datExpr2),ncol(datExpr2))
  colnames(adjmatr2) = colnames(datExpr2)
  rownames(adjmatr2) = colnames(datExpr2)
  adjmatr2[abs(cormatr2) < r_thres] = 0
  # use threshold
  for(i in 2:ncol(datExpr1)){
    #  r1 = pcor(datExpr1[,i-1],datExpr1[,i])$p.value
    #  r2 = pcor(datExpr2[,i-1],datExpr2[,i])$p.value
    if(cormatr1[i-1,i] >= p_thres){adjmatr1[i-1,i] = adjmatr1[i,i-1] = 0}
    if(adjmatr2[i-1,i] >= p_thres){adjmatr2[i-1,i] = adjmatr2[i,i-1] = 0}
  }
  #get all the basic NC
  NC1 = conformityBasedNetworkConcepts(adjmatr1)
  NC2 = conformityBasedNetworkConcepts(adjmatr2)
  #combine, rank and show
  basic_results = data.frame(density1 = NC1$fundamentalNCs$Density,
                             density2 = NC2$fundamentalNCs$Density,
                             centralization1 = NC1$fundamentalNCs$Centralization,
                             centralization2 = NC2$fundamentalNCs$Centralization,
                             heterogeneity1 = NC1$fundamentalNCs$Heterogeneity,
                             heterogeneity2 = NC2$fundamentalNCs$Heterogeneity)
  change_results = data.frame(#gene = colnames(datExpr1),
                              con_1 = NC1$fundamentalNCs$Connectivity,
                              scl_con_1 = NC1$fundamentalNCs$ScaledConnectivity,
                              con_2 = NC2$fundamentalNCs$Connectivity,
                              scl_con_2 = NC2$fundamentalNCs$ScaledConnectivity,
                              con_change = -(NC1$fundamentalNCs$Connectivity - NC2$fundamentalNCs$Connectivity),
                              scl_con_change = -(NC1$fundamentalNCs$ScaledConnectivity - NC2$fundamentalNCs$ScaledConnectivity),
                              rank_scl_con = rank(-abs(NC1$fundamentalNCs$ScaledConnectivity-NC2$fundamentalNCs$ScaledConnectivity)),
                              cls_coef_1 = NC1$fundamentalNCs$ClusterCoef,
                              cls_coef_2 = NC2$fundamentalNCs$ClusterCoef,
                              clst_coef_change = c(NC1$fundamentalNCs$ClusterCoef - NC2$fundamentalNCs$ClusterCoef),
                              rank_clstcoef = rank(-abs(NC1$fundamentalNCs$ClusterCoef-NC2$fundamentalNCs$ClusterCoef)))
  Results = list(NC1 = NC1, NC2 = NC2, cormatr_ref = cormatr1, cormatr_test = cormatr2, adjmatr_ref = adjmatr1, adjmatr_test = adjmatr2, basic = basic_results,change = change_results)
  return(Results)
}

######   apply function to the dataset datExprGT is reference and datExprGO is the test testset  #####
######   For the results we have 
######  NC (network concept) corcatr (corcoef matrix) adjmatr ( 0 or 1) and basics(network basic) change (change of the basics)

### results 
Results_sperm_cor = get_NC_cor(datExprGT,datExprGO,0.5,0.05)

#head(Results_sperm_cor$pvalue1)

### 2.connectivity and cluster coef 
######  connectivity (of each node)  vector  ##################
Con_ref= rowSums(Results_sperm_cor$adjmatr_ref) - 1 
Con_test= rowSums(Results_sperm_cor$adjmatr_test) - 1

###### mean connectivity  ##############
meanCon_ref = sum(Con_ref)/ncol(datExpr_sperm);meanCon_ref
meanCon_test = sum(Con_test)/ncol(datExpr_sperm);meanCon_test

######     density      #################
density_ref = sum(vectorizeMatrix(Results_sperm_cor$adjmatr_ref))/(0.5*ncol(datExpr_sperm)*(ncol(datExpr_sperm)-1));density_ref
density_test = sum(vectorizeMatrix(Results_sperm_cor$adjmatr_test))/(0.5*ncol(datExpr_sperm)*(ncol(datExpr_sperm)-1));density_test
######  coef (of each node)  vector  ##################
clstcoef_ref= Results_sperm_cor$NC1$fundamentalNCs$ClusterCoef
clstcoef_test= Results_sperm_cor$NC2$fundamentalNCs$ClusterCoef
###### mean coef  ##################
meanClstcoef_ref = sum(clstcoef_ref)/ncol(datExpr_sperm);meanClstcoef_ref
meanClstcoef_test = sum(clstcoef_test)/ncol(datExpr_sperm);meanClstcoef_test

### 3.top 10 of connectivity/clustercoeffcient ###
########## assemble dataset 1 _ con
topgene_sperm_con = data.frame(
  Con_ref = Results_sperm_cor$change$con_1,
  Con_rank_ref = rank(-Results_sperm_cor$change$con_1,ties.method = "min"),
  #Con_ref_scl = Results_sperm_cor$change$scl_con_1,
  Con_test = Results_sperm_cor$change$con_2,
  #Con_ref_scl = Results_sperm_cor$change$scl_con_1,
  Con_rank_test = rank(-Results_sperm_cor$change$con_2,ties.method = "min"),
  Con_Change = Results_sperm_cor$change$con_change,
  ConChange_rank =rank(-(abs(Results_sperm_cor$change$con_1 - Results_sperm_cor$change$con_2)),ties.method = "min")
)
rownames(topgene_sperm_con) = colnames(datExprGT)
######### assemble dataset 2 _ clscoef
topgene_sperm_clscoef = data.frame(
  Clscoef_ref = Results_sperm_cor$change$cls_coef_1 ,
  Clscoef_rank_ref = rank(-Results_sperm_cor$change$cls_coef_1,ties.method = "min"),
  Clscoef_test = Results_sperm_cor$change$cls_coef_2,
  Clscoef_rank_test = rank(-Results_sperm_cor$change$cls_coef_2,ties.method = "min"),
  Clscoef_Change = Results_sperm_cor$change$clst_coef_change,
  ConChange_rank =  rank(-(abs(Results_sperm_cor$change$cls_coef_1 - Results_sperm_cor$change$cls_coef_2)),ties.method = "min")
)
rownames(topgene_sperm_clscoef) = colnames(datExprGT)

#### define function 
SelectGene_un_cor = function(dataset,topnumber){
  index1 = dataset[,2] %in% c(1:topnumber)
  index2 = dataset[,4] %in% c(1:topnumber)
  index3 = dataset[,6] %in% c(1:topnumber)
  summary1 = dataset[index1,];summary1 = summary1[order(summary1[,2]),]
  summary2 = dataset[index2,];summary2 = summary2[order(summary2[,4]),]
  summary3 = dataset[index3,];summary3 = summary3[order(summary3[,6]),]
  summary = list(
    ref = summary1,
    test = summary2,
    change =summary3
  )
  return(summary)
}
##### example result
top10_con_change = rownames(SelectGene_un_cor(topgene_sperm_con,22)[[3]])
top10_ClusterCoef_change = rownames(SelectGene_un_cor(topgene_sperm_clscoef,22)[[3]])
intersect(top10_con_change,top10_ClusterCoef_change)

Results_sperm_cor$basic
         ##################################################################################################
         ###                                       2. Plotting                                      ######
         ################################################################################################

########        1. generating dataset ####################
#####  con --- test_combine_dataset
ref = data.frame(
  connectivity = as.numeric(Results_sperm_cor$NC1$fundamentalNCs$Connectivity),
  category = rep("ref",ncol(datExprGT)))
test = data.frame(
  connectivity = as.numeric(Results_sperm_cor$NC2$fundamentalNCs$Connectivity),
  category = rep("test",ncol(datExprGT)))
test_combine_dataset <- do.call('rbind', list(ref,test))
str(test_combine_dataset)
# table(test_combine_dataset$category)

### clst coef --- test_combine_dataset_clstcoef
ref_clscoef = data.frame(
  clstcoef = as.numeric(Results_sperm_cor$NC1$fundamentalNCs$ClusterCoef),
  category = rep("ref",ncol(datExprGT)))
test_clscoef = data.frame(
  clstcoef = as.numeric(Results_sperm_cor$NC2$fundamentalNCs$ClusterCoef),
  category = rep("test",ncol(datExprGT)))
test_combine_dataset_clstcoef <- do.call('rbind', list(ref_clscoef,test_clscoef))
str(test_combine_dataset_clstcoef)


# table(test_combine_dataset$category)

#########        2. plotting         ##########################
#############################################################################################
### 1. type one ---- ridge 
#plot1 = ggplot(test_combine_dataset, aes(x = connectivity, y = category)) +
#  geom_density_ridges(aes(fill = category),scale = 3) +
#  scale_fill_manual(values = c("#00AFBB", "#FC4E07"))+
#  theme_gray()+
#  theme(legend.position="None")+
#  labs(title="Distribution of Connectivity", x="Connectivity", y = "Density")+
#  theme(plot.title = element_text(hjust = 0.5))
#print(plot1)
####
#plot2 = 
#ggplot(test_combine_dataset_clstcoef, aes(x = clstcoef, y = category)) +
#  geom_density_ridges(aes(fill = category),scale = 3) +
#  scale_fill_manual(values = c("#00AFBB", "#FC4E07"))+
#  theme_gray()+
#  theme(legend.position="None")+
#  labs(title="Distribution of Cluster Coefficient", x="Connectivity", y = "Density")+
#  theme(plot.title = element_text(hjust = 0.5))

# tiff("Figure1_connectivity_ridge.tiff", width = 14,nrow = 2, height = 12, units = 'in', res = 300)
# plot_grid(plot1, plot2, align = c("h"),labels = c("A","B"), label_size= 20, label_colour = "darkgreen")
# dev.off()
#############################################################################################
#### 2. type two --- normal (ggplot)
#install.packages("extrafont")
#font_import(pattern="[C/c]omic")
#font_import(pattern="[A/a]rial")
#font_import(pattern="[C/c]alibri")
# font_import()
# loadfonts()
# fonts()
ref = data.frame(
  connectivity = as.numeric(Results_sperm_cor$NC1$fundamentalNCs$Connectivity),
  category = rep("Control",ncol(datExprGT)))
test = data.frame(
  connectivity = as.numeric(Results_sperm_cor$NC2$fundamentalNCs$Connectivity),
  category = rep("Cottonseed",ncol(datExprGT)))
test_combine_dataset <- do.call('rbind', list(ref,test))
str(test_combine_dataset)
# table(test_combine_dataset$category)

### clst coef --- test_combine_dataset_clstcoef
ref_clscoef = data.frame(
  clstcoef = as.numeric(Results_sperm_cor$NC1$fundamentalNCs$ClusterCoef),
  category = rep("Control",ncol(datExprGT)))
test_clscoef = data.frame(
  clstcoef = as.numeric(Results_sperm_cor$NC2$fundamentalNCs$ClusterCoef),
  category = rep("Cottonseed",ncol(datExprGT)))
test_combine_dataset_clstcoef <- do.call('rbind', list(ref_clscoef,test_clscoef))
str(test_combine_dataset_clstcoef)

#color = c(rep("Blue",ncol(datExprGT)),rep("Red",ncol(datExprGT)))
                   
# dat <- data.frame(xx = c(runif(100,20,50),runif(100,40,80),runif(100,0,30)),yy = rep(letters[1:3],each = 100))


plot_a = ggplot(test_combine_dataset, aes(x=connectivity, fill=category)) +
  geom_histogram(data=subset(test_combine_dataset,category  == 'Control'),binwidth=1, alpha=.6, position="identity", aes(y = ..count..), color="black") +
  geom_histogram(data=subset(test_combine_dataset,category  == 'Cottonseed'),binwidth=1, alpha=.6, position="identity", aes(y = ..count..), color="black") +
  geom_histogram(binwidth=1,alpha=0.6, position="identity", aes(y = ..count..), color="black") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 20)) +
  #scale_fill_discrete(name="Groups")+
  scale_fill_manual(name="Maternal Diet",values = c("blue", "red")) + ## name="Experimental\nCondition",breaks=c("ctrl", "trt1", "trt2"),labels=c("Control", "Treatment 1", "Treatment 2"))
  theme(legend.title = element_text(size = 8, color = "black", face = "bold",family = "Microsoft Sans Serif",vjust = 0.5, hjust = 0.5))+
  #scale_color_manual(values=c('#999999','#E69F00'),guide=FALSE)+
  # scale_colour_manual(values = c("blue", "red"))+
  # geom_density(alpha=0.6,trim = F) + values = alpha(c("blue", "red"), .3)
  #xlim(0,100)+
  geom_vline(aes(xintercept=meanCon_ref), color="blue", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=meanCon_test), color="red", linetype="dashed", size=1) +
  theme_gray()+
  theme(legend.position="top",legend.title = element_text(hjust = .5,size = 12,face = "bold"),legend.text = element_text(hjust = .5,size=9))+ # "top"  c(.8, .8)
  
  labs(title="Distribution of Connectivity", x="Connectivity", y = "Frequency")+
  theme(axis.text.x = element_text(size = 15, family = "Microsoft Sans Serif",color = "black", vjust = 0.5, hjust = 0.5))+
  theme(axis.text.y = element_text(size = 15,family = "Microsoft Sans Serif",color = "black", vjust = 0.5, hjust = 0.5))+
  theme(axis.title.x = element_text(size = 15,family = "Microsoft Sans Serif",color = "black",vjust = 0.5, hjust = 0.5))+
  theme(axis.title.y = element_text(size = 15, color = "black",family = "Microsoft Sans Serif", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(size = 20, family = "Microsoft Sans Serif",color = "black", face = "bold", vjust = 0.5, hjust = 0.5))

plot_b = ggplot(test_combine_dataset_clstcoef, aes(x=clstcoef, fill=category)) +
  geom_histogram(data=subset(test_combine_dataset_clstcoef,category  == 'Control'),binwidth=.005,alpha=0.6, position="identity", aes(y = ..count..), color="black") +
  geom_histogram(data=subset(test_combine_dataset_clstcoef,category  == 'Cottonseed'),binwidth=.005,alpha=0.6, position="identity", aes(y = ..count..), color="black")+
  geom_histogram(binwidth=.005,alpha=0.6, position="identity", aes(y = ..count..), color="black") +
  scale_fill_manual(name="Maternal Diet",values = c("blue", "red")) + ## name="Experimental\nCondition",breaks=c("ctrl", "trt1", "trt2"),labels=c("Control", "Treatment 1", "Treatment 2"))
  theme(legend.title = element_text(size = 8, color = "black", face = "bold",family = "Microsoft Sans Serif",vjust = 0.5, hjust = 0.5))+
  #  geom_density(alpha=0.6,trim = F) +
  #scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  xlim(0.1,.95)+
  #scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
  geom_vline(aes(xintercept=meanClstcoef_ref), color="blue", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=meanClstcoef_test), color="red", linetype="dashed", size=1) +
  theme_gray()+
  theme(legend.position = "top",legend.title = element_text(hjust = .5,size = 12,face = "bold"),legend.text=element_text(size=9))+ # "top"  
  #theme(legend.position="None")+
  labs(title="Distribution of Cluster Coefficient", x="Cluster Coefficient", y = "Frequency")+
  theme(axis.text.x = element_text(size = 15, family = "Microsoft Sans Serif",color = "black", vjust = 0.5, hjust = 0.5))+
  theme(axis.text.y = element_text(size = 15,family = "Microsoft Sans Serif",color = "black", vjust = 0.5, hjust = 0.5))+
  theme(axis.title.x = element_text(size = 15,family = "Microsoft Sans Serif",color = "black",vjust = 0.5, hjust = 0.5))+
  theme(axis.title.y = element_text(size = 15, color = "black",family = "Microsoft Sans Serif", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(size = 20, family = "Microsoft Sans Serif",color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))

tiff("Sperm_network_plot.tiff", width = 14, height = 12, units = 'in', res = 300)
plot_grid(plot_a, plot_b, align = c("v"),labels = c("A","B"), nrow = 2,label_size= 20, label_colour = "darkgreen")
dev.off()
getwd()



#
plot3 = ggplot(test_combine_dataset, aes(x=connectivity,color = category,fill= category)) + 
  geom_histogram(data=subset(test_combine_dataset,category  == 'Control'),binwidth=1,alpha=0.6, position="identity", aes(y = ..count..), color="black",fill = c("blue")) +
  geom_histogram(data=subset(test_combine_dataset,category  == 'Cottonseed'),binwidth=1,alpha=0.6, position="identity", aes(y = ..count..), color="black",fill = c("red")) +
  #geom_histogram(binwidth=1,alpha=0.6, position="identity", aes(y = ..count..), color="black",fill = c("green")) + #
   #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  #scale_color_hue(l=40, c=35)+
  #scale_color_brewer(palette="Pastel2")+
  #scale_fill_brewer(palette="Pastel2")+
 # geom_density(alpha=0.6,trim = F) +
  #xlim(0,80)+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 20)) +
  geom_vline(aes(xintercept=meanCon_ref), color="black", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=meanCon_test), color="blue", linetype="dashed", size=1) +
  theme_gray()+
  #scale_fill_discrete(name="Groups")+
  theme(legend.position = "top")+ # ,legend.title = element_text(hjust = .5,size = 12,face = "bold"),legend.text=element_text(size=9))+ # "top"  c(.8, .8)
  
  labs(title="Distribution of Connectivity", x="Connectivity", y = "Frequency")+
  theme(axis.text.x = element_text(size = 15, family = "Microsoft Sans Serif",color = "black", vjust = 0.5, hjust = 0.5))+
  theme(axis.text.y = element_text(size = 15,family = "Microsoft Sans Serif",color = "black", vjust = 0.5, hjust = 0.5))+
  theme(axis.title.x = element_text(size = 15,family = "Microsoft Sans Serif",color = "black",vjust = 0.5, hjust = 0.5))+
  theme(axis.title.y = element_text(size = 15, color = "black",family = "Microsoft Sans Serif", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(size = 20, family = "Microsoft Sans Serif",color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))
plot3

plot4 = ggplot(test_combine_dataset_clstcoef, aes(x=clstcoef, fill=category)) +
  #geom_histogram(binwidth=.006,alpha=0.6, position="identity", aes(y = ..count..), color="black") +
  geom_histogram(data=subset(test_combine_dataset_clstcoef,category  == 'Control'),binwidth=.005,alpha=0.6, position="identity", aes(y = ..count..), color="black",fill = c("blue")) +
  geom_histogram(data=subset(test_combine_dataset_clstcoef,category  == 'Cottonseed'),binwidth=.005,alpha=0.6, position="identity", aes(y = ..count..), color="black",fill = c("red"))+
#  geom_density(alpha=0.6,trim = F) +
  #scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  xlim(0.1,.9)+
  #scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
  geom_vline(aes(xintercept=meanClstcoef_ref), color="black", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=meanClstcoef_test), color="blue", linetype="dashed", size=1) +
  theme_gray()+
  theme(legend.position = "top",legend.title = element_text(hjust = .5,size = 12,face = "bold"),legend.text=element_text(size=9))+ # "top"  
  #theme(legend.position="None")+
  labs(title="Distribution of Cluster Coefficient", x="Cluster Coefficient", y = "Frequency")+
  theme(axis.text.x = element_text(size = 15, family = "Microsoft Sans Serif",color = "black", vjust = 0.5, hjust = 0.5))+
  theme(axis.text.y = element_text(size = 15,family = "Microsoft Sans Serif",color = "black", vjust = 0.5, hjust = 0.5))+
  theme(axis.title.x = element_text(size = 15,family = "Microsoft Sans Serif",color = "black",vjust = 0.5, hjust = 0.5))+
  theme(axis.title.y = element_text(size = 15, color = "black",family = "Microsoft Sans Serif", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(size = 20, family = "Microsoft Sans Serif",color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))
plot4
dev.off()
#test1 
ggarrange(plot3, plot4, ncol=1, nrow=2, common.legend = TRUE, legend="top")



tiff("test1.tiff", width = 14, height = 12, units = 'in', res = 300)

#ggarrange(plot3, plot4, ncol=1, nrow=2, common.legend = TRUE, legend="top")

plot_grid(plot3, plot4, align = c("v"),labels = c("A","B"), nrow = 2,label_size= 20, label_colour = "darkgreen")

dev.off()
getwd()
# getwd()

#########        3.differential expresion and differential connected    ##########################
##################################################################################################
########## put togehter all the index for selection  #############
ScreenDataset_change = data.frame(
  LogFC = round(networkData_sperm$logFC,4),
  #scale_FC =  round(abs(networkData_sperm$logFC)/(max(abs(networkData_sperm$logFC))),4),
  Con_Change = round(Results_sperm_cor$change$con_change,4),
  scale_ConChange = round(Results_sperm_cor$change$con_change/(max(abs(Results_sperm_cor$change$con_change))),4),
  Clscoef_Change = round(Results_sperm_cor$change$clst_coef_change,4),
  scale_ClscoefChange = round(Results_sperm_cor$change$clst_coef_change/(max(abs((Results_sperm_cor$change$clst_coef_change)))),4),
  index = rep("No", ncol(datExprGT))
)
rownames(ScreenDataset_change) = colnames(datExprGT)
### get DEs 
DE_index = which(networkData_sperm$Significant == "Yes")
ScreenDataset_change[DE_index,"index"] = "Yes"

write.csv(ScreenDataset_change,file = "ScreenDataset_change2.csv",row.names = T)

# head(ScreenDataset_change,6)

################                      plotly                         ###################
#library(tidyverse)
#library(plotly)
# packageVersion('plotly')

##########
#str(ScreenDataset_change)
ScreenDataset_change$index = factor(ScreenDataset_change$index)
threed_sperm = plot_ly(ScreenDataset_change, x = ~scale_ConChange, y = ~scale_ClscoefChange, z = ~LogFC, 
        type = "scatter3d", mode = "markers",
        marker = list(opacity = 1, size = 3),
        color = ~index, colors = c( '#BF382A','#0C4B8E'), 
        showlegend = T,
        alpha = 0.8) %>%
  add_markers() %>%
  layout(
    scene = list(camera = list(eye = list(x = -1.25, y = 1.25, z = .15)),
                 xaxis = list(title = 'Con_Change',range = c(-1,1)),
                 yaxis = list(title = 'Clscoef_Change',range = c(-1,1)),
                 zaxis = list(title = 'Log(FC)',range = c(1.2*min(ScreenDataset_change$LogFC),1.2*max(ScreenDataset_change$LogFC)))),
    plot_bgcolor=c('rgb(254, 247, 234)'),
    paper_bgcolor=c('rgb(254, 247, 234)'),
    showlegend = FALSE
  )

#chart_link = api_create(threed_sperm, filename = "3dgraph_sperm")
#chart_link
######################

#########    organize dataset and try to plot    ##########
### convert ensembl id to gene name  ###########
ID2Name_all = read.csv("ToFindGeneName.txt")
Gene_Name = ID2Name_all$Gene.name[ID2Name_all$Gene.stable.ID  %in% names(datExpr_sperm)]
allID = data.frame(Ensembl_ID = names(datExpr_sperm),Gene_Names = Gene_Name)
for (i in c(1:length(allID$Gene_Names))){
   if (allID$Gene_Names[i] == "") {
     allID$Gene_Names[i] <- allID$Ensembl_ID[i]
   }
}
write.csv(allID,"ALL id ready to gene name.txt",quote = F)

# get edgelist structure
library(igraph)
adjmatr1 = Results_sperm_cor$adjmatr_ref
adjmatr2 = Results_sperm_cor$adjmatr_test
rownames(adjmatr1) = allID$Ensembl_ID;colnames(adjmatr1) = allID$Ensembl_ID
rownames(adjmatr2 )= allID$Ensembl_ID;colnames(adjmatr2) = allID$Ensembl_ID

Top50_index_order = order(-abs(Results_sperm_cor$change$con_change))[1:32]
Top50_index1 = allID$Ensembl_ID[Top50_index_order]
Top50_index2 = allID$Gene_Names[Top50_index_order]
DE_index_posi = which(networkData_sperm$Significant == "Yes")
DE_index = allID$Gene_Names[DE_index_posi]

# check in the data set
# table(rownames(Results_sperm_cor$change[Top50_index_order,]) == Top50_index1)
adjmatr_top50_1 = adjmatr1[Top50_index_order,Top50_index_order]
adjmatr_top50_2 = adjmatr2[Top50_index_order,Top50_index_order]
diag(adjmatr_top50_1) = 0
diag(adjmatr_top50_2) = 0

Net_top50_1 = graph.adjacency(adjmatr_top50_1,mode = "undirected", weighted = NULL)
Net_top50_2 = graph.adjacency(adjmatr_top50_2,mode = "undirected", weighted = NULL)

write_graph(Net_top50_1, "Net_ref_sperm.text", format = "edgelist")
write_graph(Net_top50_2, "Net_test_sperm.text", format = "edgelist")

#####
allID_sub_annot = cbind(allID[Top50_index_order,], index = factor(seq(1,32,by= 1),ordered = F))

###### re-get the structures and add annotation using merge()    ######
Net_ref = read.csv("Net_ref_sperm.text",sep = "",col.names = c("source","end"))
Net_test = read.csv("Net_test_sperm.text",sep = "",col.names = c("source","end"))
Net_ref = Net_ref + 1 ; Net_test = Net_test +1

# Net_ref$source = as.factor(Net_ref$source);Net_ref$end = as.factor(Net_ref$end)
# Net_test$source = as.factor(Net_test$source);Net_test$end = as.factor(Net_test$end)


## 2 extra dataset just for substitude 
allID_sub_annot_1 = allID_sub_annot = cbind(allID[Top50_index_order,], source = factor(seq(1,32,by= 1),ordered = F))
allID_sub_annot_2 = allID_sub_annot = cbind(allID[Top50_index_order,], end = factor(seq(1,32,by= 1),ordered = F))
allID_sub_annot = cbind(allID[Top50_index_order,], index = factor(seq(1,32,by= 1),ordered = F))

step1 = merge(allID_sub_annot_1, Net_ref, by ="source")
step2 = merge(allID_sub_annot_2, step1, by ="end")
step3 = merge(allID_sub_annot_1, Net_test, by ="source")
step4 = merge(allID_sub_annot_2, step3, by ="end")
Final_edgelist_ref = step2;Final_edgelist_test = step4

# unique(Final_edgelist_ref$Gene_Names.x)
# unique(Final_edgelist_ref$Gene_Names.y)

### finally get final edgelist  ####
write.csv(Final_edgelist_ref, "Final_edgelist_cyto_ref.txt", quote = F)
write.csv(Final_edgelist_test, "Final_edgelist_cyto_test.txt", quote = F)
write.csv(Top50_index1, "Top50_index.txt", quote = F)

DE_inTop50_intersect = DE_index[DE_index_posi %in% Top50_index_order]



#### get subnetwork results to check plot sequence  #####

cyto_fdback_ref = read.csv("Ref_default_node_2.csv")
node_annot_trans = data.frame(Gene_Names.y = cyto_fdback_ref$name, anno_index = cyto_fdback_ref$label_posi)
Final_edgelist_test_index = data.frame(Gene_Names.x = Final_edgelist_test$Gene_Names.x, Gene_Names.y = Final_edgelist_test$Gene_Names.y)
Final_edgelist_test_withindex = merge(Final_edgelist_test_index,node_annot_trans,by ="Gene_Names.y")
write.csv(Final_edgelist_test_withindex, "Final_edgelist_test_withindex.txt", quote = F)

rownames(Results_sperm_cor$change[,c(5,7)][order(abs(Results_sperm_cor$change[7])),])
table(allID$Ensembl_ID[order(abs(Results_sperm_cor$change[7]))] == rownames(Results_sperm_cor$change[,c(5,7)][order(abs(Results_sperm_cor$change[7])),]))
names(Results_sperm_cor$change)
cyto_fdback_ref = cyto_fdback_ref[order(-cyto_fdback_ref$degree.layout),]intersect(DE_index,cyto_fdback_ref$name)
(Net_top50_1)
## check if true
allID$Gene_Names[which(networkData_sperm$Gene=="ENSOARG00000005059")]
# length(unique(c(unique(Final_edgelist_ref$Gene_Names.x),unique(Final_edgelist_ref$Gene_Names.x))))
cyto_fdback_ref = read.csv("Final_ref_default_node.csv")
table(cyto_fdback_ref$name %in% Top50_index2)
  
######## wirte files for all the con/clustercoef/change / pvalue
######## 
# correlation coefficient
corcoef_out1 = Results_sperm_cor$cormatr_ref
corcoef_out2 = Results_sperm_cor$cormatr_test
pval_out1 = Results_sperm_cor$pvalue1
pval_out2 = Results_sperm_cor$pvalue2
adjmatr_out1 = Results_sperm_cor$adjmatr_ref
adjmatr_out2 = Results_sperm_cor$adjmatr_test
Con_rank_out0 = Results_sperm_cor$change[,c(1,3,5)]
Con_change_out = cbind(Gene_names = allID$Gene_Names,Con_rank_out0)
ClsCoef_change_out0 = Results_sperm_cor$change[,c(8:10)]
ClsCoef_change_out = cbind(Gene_names = allID$Gene_Names,ClsCoef_change_out0)

list_df <- list("CorCoef matrix Reference" = corcoef_out1, 
                "CorCoef matrix Test" = corcoef_out2,
                "P value matrix Reference" = pval_out1,
                "p value matrix Test" = pval_out2,
                "Adjacency matrix Reference" = adjmatr_out1,
                "Adjacency matrix Test" = adjmatr_out2,
                "Con basics" = Con_change_out,
                "Cls Coef basics" = ClsCoef_change_out)

write.xlsx(list_df, file = "Network Statistics.xlsx")

######################################################################################################
### (with differential expressed genes)top gene test ###
Topgene_test = function(n){
  Topgege_test_index_order = order(-abs(Results_sperm_cor$change$con_change))[1:n]
  Re = table(DE_index_posi %in% Topgege_test_index_order)
  return(Re)
}
Topgene_test(50)
######################################################################################################


?graph.adjacency()
#E(Net);V(Net)
##plotting
l <- layout_in_circle(Net1)
par(mfrow=c(1,2))
plot(Net1, layout=l);plot(Net2, layout=l)
##
target = 36
store = c((vcount(Net1)-target+2):vcount(Net1),1:(vcount(Net1)+1-target))
l_one <- cbind(1:vcount(Net1), c(1, vcount(Net1):2));l_one = l_one[store,]
par(mfrow=c(1,2))
plot(Net1, layout=l_one);plot(Net2, layout=l_one);

?write_graph()

net1_text = read.table("Net1.text")
head(net1_text,16)

Net2_NeedRank = read.csv("Net2_NeedRank.csv")

str(Net2_NeedRank)

length(table(Net2_raw$V1))

Net1_add_Rank = cbind(Net1_NeedRank[order(Net1_NeedRank$degree.layout,decreasing = T),],rev(seq(69:1)))
Net1_add_Rank 


#Net1_test = Net1_add_Rank[,c("name","rev(seq(69:1))")]
Net1_test = Net1_add_Rank[,c("name","degree.layout")]
rownames(Net1_test) = Net1_test$name

Net2_raw = read.table("Net2.text")
#Net2_new = merge(Net2_raw,Net1_test,by.x ="V1",by.y="name")
Net2_new = merge(Net2_raw,Net1_test,by.x ="V1",by.y="name")
Net2_new

table(Net2_raw$V1)
table(Net2_raw$V2)

table(Net2_new$degree.layout)

write.csv(Net2_new,"Net2.new_add_old_con.csv")

?merge(Net1_raw,Net1_test,by.x ="V1",by.y="name")

?join()
?match()
names(Net1_add_Rank)

test_rank = read.csv("rank_testing_5216.csv",sep = " ")
?read.csv()
getwd()


####     Module preservation aspect  #####
Con_presv = cor(vectorizeMatrix(as.matrix(Results_sperm_cor$adjmatr_ref)),as.matrix(vectorizeMatrix(Results_sperm_cor$adjmatr_test)))
Con_presv
#Clstcoef_presv = cor(Results_sperm_cor$NC1$fundamentalNCs$ClusterCoef,Results_sperm_cor$NC2$fundamentalNCs$ClusterCoef)
#Clstcoef_presv

### cor.cor ###
cor.cor_sperm = cor(vectorizeMatrix(Results_sperm_cor$cormatr_ref),vectorizeMatrix(Results_sperm_cor$cormatr_test))
cor.cor_sperm
?modulePreservation()

MP_sperm$observed$sperm_GT$intra
table(Results_sperm_cor$cormatr_test == Results_sperm_cor$cormatr_ref)

modulePreservation(Results_sperm_cor)
names(Results_sperm_cor)
dim(as.matrix(Results_sperm_cor$adjmatr1))
###############################
#####differences between network properties #######
###############################

#######change of con general ---overall plotting#####
####show the change in hist: distribution #########
par(mfrow=c(3,1))
(Results_sperm_cor$NC1$fundamentalNCs$Connectivity)
Results_sperm_cor$NC1$fundamentalNCs$Density

mean(sum(vectorizeMatrix(Results_sperm_cor$adjmatr_ref)))
sum(vectorizeMatrix(Results_sperm_cor$adjmatr_ref))/145
sum(Results_sperm_cor$NC1$fundamentalNCs$ClusterCoef)/145
sum(Results_sperm_cor$NC2$fundamentalNCs$ClusterCoef)/145

head(Results_sperm_cor$adjmatr_ref)

#hist_C1 = hist(genCon_C1,breaks = 25)
#hist_M1 = hist(genCon_M1,breaks = 25)

(vectorizeMatrix(Results_sperm_cor$adjmatr_ref))
(colnames(Results_sperm_cor$adjmatr_ref))

genCon_GT = rowSums(Results_sperm_cor$adjmatr_ref) - 1 

genCon_GO = rowSums(Results_sperm_cor$adjmatr_test) - 1 
cor(genCon_GT,genCon_GO)


### rest 
hist_ref = hist(Results_sperm_cor$NC1$fundamentalNCs$ClusterCoef,breaks = 20)
hist_test = hist(Results_sperm_cor$NC2$fundamentalNCs$ClusterCoef,breaks = 20)
plot( hist_ref, col=rgb(0,1/2,1,1/4), xlim=c(0,1),ylim = c(1,100), main = paste("clstcoef distbt of control/treatment"));plot( hist_test, col=rgb(1,0,0,1/4), xlim=c(0,1),ylim = c(1,150),add=T);legend("topright", c("Control", "treatment"), col=c(rgb(0,1/2,1,1/4), rgb(1,0,0,1/4)), lwd=10)
#  at=seq(0,1,.1),
hist_ref = hist(Results_sperm_cor$NC1$fundamentalNCs$ClusterCoef,breaks = 20)
hist_test = hist(Results_sperm_cor$NC2$fundamentalNCs$ClusterCoef,breaks = 20)
plot( hist_ref, col=rgb(0,1/2,1,1/4), xlim=c(0,150),ylim = c(1,40), at=seq(1,110,1),main = paste("Con distbt of C3/M3 _S"));plot( hist_test, col=rgb(1,0,0,1/4), xlim=c(0,150) ,ylim = c(1,100),add=T);legend("topright", c("Control", "treatment"), col=c(rgb(0,1/2,1,1/4), rgb(1,0,0,1/4)), lwd=10)


cor(vectorizeMatrix(Results_sperm_cor$adjmatr_ref,diag = F),vectorizeMatrix(Results_sperm_cor$adjmatr_test,diag = F))
cor(Results_sperm_cor$NC1$fundamentalNCs$Connectivity,Results_sperm_cor$NC2$fundamentalNCs$Connectivity)

vectorizeMatrix(Results_sperm_cor$adjmatr_ref)
length(vectorizeMatrix(Results_sperm_cor$adjmatr_ref))


#######     now we have function to get basics (change of the basics) using cor and pcor      ########
######      basic (basic NC) - change of connectivity and clustercoef (by ranking)           #########
######      declare the parameters you want to use in the fucntion parameters - e.g.cor CUTOFF 0.5  ##
######################################################################################################

Results_sperm_pcor = get_NC_pcor(datExprGT,datExprGO,0.5,0.01)
######         Results_sperm_cor$ change has the need for plotting in cyto       ######################


###########      take the "basics" and do formating for the Cytoscape input        ###############
############                add "rank" for plotting                                ###############
###             2.Define formating function - use package "igraph"             ###################

###NEED
###TO 
###BE 
###DONE

tiff("Figure_int_slo_6_14chr_20190313.tiff", width = 14, height = 12, units = 'in', res = 300)
plot_grid(plot1, plot2,plot3,plot4, plot5,plot6, align = c("hv"), nrow = 3,  
          labels = c("A", "B","C","D","E","F"), label_size= 20, label_colour = "darkgreen")
dev.off()


######################################################################################################################

##### overviews of the plots, specific nodes(maybe hub nodes) may be selected for plotting?
##### basic stats : density/centralization/heterogeneity;
##### the change of connectivity (scaled by the mixmum), as well as the corresponding rank (of abs value); also, the same rationale for clustering coefficient (the density of neighboors-conncection of a node)
##### connectivity and clustering coefficient may not necessary to be conformative, potential composite stats maybe selected/proposed?

##########################################################################################
Net = graph.adjacency(dissTOMGO,mode = "undirected",weighted = TRUE)
l <- layout_in_circle(Net)
#par(mfrow=c(1,2))
plot(Net, layout=circle)
plot(Net, layout=l)
#summary(dissTOMC1)
##########################################################################################



################################################################################################################
################################################################################################################
#################################################################################################################
#=====================================================================================
#  Code chunk 5 - weighted correlation network ---set up (automatic)
#=====================================================================================
##########################################################################################
############                     r weighted                    #####################
powers = c(c(1:10), seq(from = 12, to=20, by=1))
# Call the network topology analysis function
sft_GT = pickSoftThreshold(datExprGT, networkType = "unsigned",powerVector = powers, verbose = 5)

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft_GT$fitIndices[,1], -sign(sft_GT$fitIndices[,3])*sft_GT$fitIndices[,2],
     xlab="Soft Threshold (power) GT",ylab="Scale Free Topology Model Fit,signed R^2, GT",type="n",
     main = paste("Scale independence GT"));text(sft_GT$fitIndices[,1], -sign(sft_GT$fitIndices[,3])*sft_GT$fitIndices[,2],
                                                 labels=powers,cex=cex1,col="red");abline(h=0.80,col="red")# this line corresponds to using an R^2 cut-off of h
# Mean connectivity as a function of the soft-thresholding power
plot(sft_GT$fitIndices[,1], sft_GT$fitIndices[,5],xlab="Soft Threshold (power) GT",ylab="Mean Connectivity GT", type="n",main = paste("Mean connectivity GT"));text(sft_GT$fitIndices[,1], sft_GT$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
#######         TAKE power 8, fitting index reach 0.8         ##########################

##############calculate the critical values that will be used for analysis ##########################
softPower = 8
adjacencyGT = adjacency(datExprGT,power=softPower,type="unsigned");
diag(adjacencyGT)=0
dissTOMGT = 1-TOMsimilarity(adjacencyGT, TOMType="unsigned")
adjacencyGO = adjacency(datExprGO,power=softPower,type="unsigned");
diag(adjacencyGO)=0
dissTOMGO = 1-TOMsimilarity(adjacencyGO, TOMType="unsigned")

geneTreeGT = hclust(as.dist(dissTOMGT), method="average")
geneTreeGO = hclust(as.dist(dissTOMGO), method="average")

###############                visualization                           ##############################
#pdf("dendrogram.pdf",height=6,width=16)
par(mfrow=c(1,2))
plot(geneTreeGT,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity (GT)",
     labels=FALSE,hang=0.04);
plot(geneTreeGO,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity (GO)",
     labels=FALSE,hang=0.04);
dev.off()

########### for now we DO NOT need module assignment, all genes in this go is in one module  #######
####        Next we will determine modules based on dataset GT                      ###########
##########  Color --- all gray  (only one module)                              ##########
mColorhGT=NULL
treeGT = cutreeHybrid(dendro = geneTreeGT, pamStage=FALSE,
                      minClusterSize = 5, cutHeight = 0.0001,
                      deepSplit = 0, distM = dissTOMGT)
mColorhGT=cbind(mColorhGT,labels2colors(treeGT$labels))


table(treeGT)

#for (ds in 0:1){
#    treeC3 = cutreeHybrid(dendro = geneTreeC3, pamStage=FALSE,#
#                     minClusterSize = (30-1*ds), cutHeight = 0.99,
#                      deepSplit = ds, distM = dissTOMC3)
#    mColorhC3=cbind(mColorhC3,labels2colors(treeC3$labels));
#}
#   pdf("Module_choices.pdf", height=10,width=25);
plotDendroAndColors(geneTreeGT, mColorhGT, main = "",dendroLabels=F);
dev.off()
modulesGT = mColorhGT[,1] # (Chosen based on plot below)

#length(modulesGT)
###################################################################################################
####                          NC based on weighted measure                                    #####
get_NC_wt = function(datExpr1,datExpr2,softPower,signORunsign){
  #calcu of matrx1
  adjacency1 = adjacency(datExpr1,power=softPower,type=as.character(signORunsign));
  diag(adjacency1)=0
  TOMsimilarity1 = TOMsimilarity(adjacency1, TOMType=as.character(signORunsign))
  dissTOMGT1 = 1 - TOMsimilarity1
  #calcu of matrx2
  adjacency2 = adjacency(datExpr2,power=softPower,type=as.character(signORunsign));
  diag(adjacency2)=0
  TOMsimilarity2 = TOMsimilarity(adjacency2, TOMType=as.character(signORunsign))
  dissTOMGT2 = 1 - TOMsimilarity2
  #get all the basic NC
  NC1 = conformityBasedNetworkConcepts(TOMsimilarity1)
  NC2 = conformityBasedNetworkConcepts(TOMsimilarity2)
  #combine, rank and show
  basic_results = data.frame(density1 = NC1$fundamentalNCs$Density,
                             density2 = NC2$fundamentalNCs$Density,
                             centralization1 = NC1$fundamentalNCs$Centralization,
                             centralization2 = NC2$fundamentalNCs$Centralization,
                             heterogeneity1 = NC1$fundamentalNCs$Heterogeneity,
                             heterogeneity2 = NC2$fundamentalNCs$Heterogeneity)
  change_results = data.frame(gene = colnames(datExpr1),
                              scl_con_1 = NC1$fundamentalNCs$ScaledConnectivity,
                              scl_con_2 = NC2$fundamentalNCs$ScaledConnectivity,
                              scl_con_change =NC1$fundamentalNCs$ScaledConnectivity - NC2$fundamentalNCs$ScaledConnectivity,
                              rank_scl_con = rank(-abs(NC1$fundamentalNCs$ScaledConnectivity-NC2$fundamentalNCs$ScaledConnectivity)),
                              cls_coef_1 = NC1$fundamentalNCs$ClusterCoef,
                              cls_coef_2 = NC2$fundamentalNCs$ClusterCoef,
                              clst_coef_change = c(NC1$fundamentalNCs$ClusterCoef - NC2$fundamentalNCs$ClusterCoef),
                              rank_clstcoef = rank(-abs(NC1$fundamentalNCs$ClusterCoef-NC2$fundamentalNCs$ClusterCoef)))
  Results = list(basic = basic_results,change = change_results, similarity_ref = TOMsimilarity1, similarity_test = TOMsimilarity1, adjmatr1=adjacency1,adjmatr2=adjacency2 )
  return(Results)
}
Results_sperm = get_NC_wt(datExprGT,datExprGO,8,"unsigned")
names(Results_sperm)

# TOMsimilarity_GO = TOMsimilarity(adjacencyGO, TOMType="unsigned")
# TOMsimilarity_GT = TOMsimilarity(adjacencyGT, TOMType="unsigned")
# NC_GO_w=conformityBasedNetworkConcepts(TOMsimilarity_GO)
# NC_GT_w=conformityBasedNetworkConcepts(TOMsimilarity_GT)


####  steps for detail -- although included in function ######
#   adjacencyGT = adjacency(datExprGT,power=softPower,type="unsigned");
#   diag(adjacencyGT)=0
#   dissTOMGT = 1-TOMsimilarity(adjacencyGT, TOMType="unsigned")
#   adjacencyGO = adjacency(datExprGO,power=softPower,type="unsigned");
#   diag(adjacencyGO)=0
#   dissTOMGO = 1-TOMsimilarity(adjacencyGO, TOMType="unsigned")

geneTreeGT = hclust(as.dist(dissTOMGT), method="average")
geneTreeGO = hclust(as.dist(dissTOMGO), method="average")
###############################################################

######    now!!! do the preservatino statistics      ################
##########    To quantify this result    module preservation statistics    ############################

## substitiue gray with blue  
mColorhGT = gsub("grey", "blue", mColorhGT)
rownames(mColorhGT) = networkData_sperm$Gene
#dim(mColorhGT) ; str(mColorhGT)

multiExpr_sperm = list(GT=list(data=datExprGT),GO=list(data=datExprGO))
#multiColor_sperm = list(GT = modulesGT)
multiColor_sperm = list(GT = mColorhGT,GO = mColorhGT)

MP_sperm=modulePreservation(multiExpr_sperm,multiColor_sperm,referenceNetworks=1,verbose=3,networkType="unsigned",
                       nPermutations=30,maxGoldModuleSize=30,maxModuleSize=145)

MP_sperm=modulePreservation(multiExpr_sperm,multiColor_sperm,referenceNetworks=1,verbose=3,networkType="unsigned",
                            nPermutations=0)

????(MP_sperm$observed)

MP_sperm$observed$sperm_GT$netPresent

?modulePreservation()
names(MP_sperm)

stats = mp1$preservation$Z$ref.C1$inColumnsAlsoPresentIn.M1
stats[order(-stats[,2]),c(1:2)]

#write.csv(stats[order(-stats[,2]),c(1:2)],"module size and pres")





#############################################################################################################################################
mergingThresh = 0.25
net = blockwiseModules(datExprGT,corType="pearson",
                       maxBlockSize=5000,networkType="unsigned",power=8,minModuleSize=145,
                       mergeCutHeight=0.0001,numericLabels=TRUE,saveTOMs = F,
                       pamRespectsDendro=FALSE,saveTOMFileBase="TOM_sperm_GT")
moduleLabelsAutomatic=net$colors
moduleLabelsAutomatic = gsub(0,1,moduleLabelsAutomatic)
# Convert labels to colors for plotting
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)


# A data frame with module eigengenes can be obtained as follows
#       MEsAutomatic=net$MEs

#this is the body weight
weight = as.data.frame(datTraits$weight_g)
names(weight)="weight"
# Next use this trait to define a gene significance variable
GS.weight=as.numeric(cor(datExprFemale,weight,use="p"))
# This translates the numeric values into colors
GS.weightColor=numbers2colors(GS.weight,signed=T)
blocknumber=1
datColors=data.frame(moduleColorsAutomatic,GS.weightColor)[net$blockGenes[[blocknumber]],]

# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[blocknumber]],colors=datColors,
                    groupLabels=c("Module colors","GS.weight"),dendroLabels=FALSE,
                    hang=0.03,addGuide=TRUE,guideHang=0.05)
























#############################################################################################################################################



adjmatrx_sperm1 = Results_sperm$adjmatr1
diag(adjmatrx_sperm1) = 1

change = as.data.frame(Results_sperm$change)


labels_sperm = substr(Results_sperm$change$gene,nchar(Results_sperm$change$gene[1])-2,nchar(Results_sperm$change$gene[1]))




#nchar(Results_sperm$change$gene[1])

change$scl_con_1[36]
change$gene[]
order1 = rank(-change$scl_con_1)
order2 = rank(-change$scl_con_2)

circlePlot(adjmatrx_sperm1, labels_sperm, order1, startNewPlot = T, 
           variable.cex.labels = FALSE, center = c(0.5, 0.5), 
           radii = c(0.35, 0.35))

  



circlePlot(Results_sperm$adjmatr1, labels_sperm, order1, startNewPlot = T, 
           variable.cex.labels = FALSE, center = c(0.5, 0.5), 
           radii = c(0.55, 0.55))

dev.off()







#############################################################################################################################################
# We now set up the multi-set expression data
# and corresponding module colors:
setLabels = c("sperm_GT", "sperm_GO")
datExprGT_mtrx = as.matrix(datExprGT)
datExprGO_mtrx = as.matrix(datExprGO)

str()



multiExpr_sperm=list(sperm_GT=list(data=datExprGT_mtrx),
               sperm_GO=list(data=datExprGO_mtrx))
#moduleColorsGT=moduleColorsAutomatic
moduleColorsGT = c(1:145); moduleColorsGO = c(1:145)
multiColor_sperm=list(sperm_GT=moduleColorsGT,sperm_GO=moduleColorsGO)

names(multiExpr_sperm); names(multiColor_sperm)


# The number of permutations drives the computation time
# of the module preservation function. For a publication use 200 permutations.
# But for brevity, let's use a small number
nPermutations1=10
# Set it to a low number (e.g. 3) if only the medianRank statistic
# and other observed statistics are needed.
# Permutations are only needed for calculating Zsummary
# and other permutation test statistics.
# set the random seed of the permutation test analysis
set.seed(1)
system.time({
  MP_sperm = modulePreservation(multiExpr_sperm, multiColor_sperm,
                                referenceNetworks = c(1:2),
                                nPermutations = nPermutations1,
                                randomSeed = 1, 
                                quickCor = 0, 
                                verbose = 3)
})


# Save the results of the module preservation analysis
save(mp, file = "modulePreservation.RData")
# If needed, reload the data:
load(file = "modulePreservation.RData")

# specify the reference and the test networks
ref=1; test = 2


Obs.PreservationStats= MP_sperm$preservation$observed[[ref]][[test]]
Z.PreservationStats=MP_sperm$preservation$Z[[ref]][[test]]
# Look at the observed preservation statistics
Obs.PreservationStats
########################################################################################################################



###############################################################
#    The first PC is referred to as the module eigengene (ME), and is a single value that
#   represents the highest percent of variance for all genes in a module.
#
PCs_C1 = moduleEigengenes(datExprC1, colors=modulesC1)
ME_C1 = PCs_C1$eigengenes
distPCC1 = 1-abs(cor(ME_C1,use="p"))
distPCC1 = ifelse(is.na(distPCC1), 0, distPCC1)
pcTreeC1 = hclust(as.dist(distPCC1),method="a")
MDS_C1 = cmdscale(as.dist(distPCC1),2)
colorsC1 = names(table(modulesC1))
#####
# PCs_M3 = moduleEigengenes(datExprM3, colors=modulesM3)
# ME_M3 = PCs_M3$eigengenes
# distPCM3 = 1-abs(cor(ME_M3,use="p"))
# distPCM3 = ifelse(is.na(distPCM3), 0, distPCM3)
# pcTreeM3 = hclust(as.dist(distPCM3),method="a")
# MDS_M3 = cmdscale(as.dist(distPCM3),2)
# colorsM3 = names(table(modulesM3))


###  serious of plots####  maybe not useful in our ananlysis
#save.image("tutorial.RData")
#pdf("ModuleEigengeneVisualizations.pdf",height=6,width=6)
par(mfrow=c(1,1), mar=c(0, 3, 1, 1) + 0.1, cex=1)
plot(pcTreeC1, xlab="",ylab="",main="",sub="")
plot(MDS_C1, col= colorsC1, main="MDS plot", cex=2, pch=19)
ordergenesC1 = geneTreeC1$order
plotMat(scale(log(t(datExprC1)[ordergenesC1,])) , rlabels= modulesC1[ordergenesC1], clabels=
          colnames(t(datExprC1)), rcols=modulesC1[ordergenesC1])
for (which.module in names(table(modulesC1))){
  ME = ME_C1[, paste("ME",which.module, sep="")]
  barplot(ME, col=which.module, main="", cex.main=2,
          ylab="eigengene expression",xlab="array sample")
};
dev.off()


#####################################################################################################
#####Step 4: Qualitatively and quantitatively measure network preservation at the module level#####
#####
# pdf("Final_modules.pdf",height=8,width=12)
par(mfrow=c(3,1))

plotDendroAndColors(geneTreeC1, modulesC1, "Modules", dendroLabels=F, hang=0.03, addGuide=TRUE,
                    guideHang=0.05, main="Gene dendrogram and module colors (C1)")
plotDendroAndColors(geneTreeM1, modulesC1, "Modules", dendroLabels=FALSE, hang=0.03, addGuide=TRUE,
                    guideHang=0.05, main="Gene dendrogram and module colors (M1)")

dev.off()

####The "grey" module contains uncharacterized gene while the gold module contains random genes.


###  We first will get the kME values, along with their associated p-values for A1 
###  and will then output the resulting
###  table to a file ("kMEtable1.csv").

geneModuleMembershipC1 = signedKME(datExprC1, ME_C1)
colnames(geneModuleMembershipC1)=paste("PC",colorsC1,".cor",sep="");
MMPvalueC1=corPvalueStudent(as.matrix(geneModuleMembershipC1),dim(t(datExprC1))[[2]]);
colnames(MMPvalueC1)=paste("PC",colorsC1,".pval",sep="");
Gene = rownames(t(datExprC1))
kMEtableC1 = cbind(Gene,Gene,modulesC1)
for (i in 1:length(colorsC1))
  kMEtableC1 = cbind(kMEtableC1, geneModuleMembershipC1[,i], MMPvalueC1[,i])
colnames(kMEtableC1)=c("PSID","Gene","Module",sort(c(colnames(geneModuleMembershipC1),
                                                     colnames(MMPvalueC1))))
write.csv(kMEtableC1,"kMEtableC1.csv",row.names=FALSE)

###   Now repeat for A2, using the module assignments from A1 to determine kME values.
###   First calculate MEs for A2, since we haven't done that yet
PCsM1 = moduleEigengenes(datExprM1, colors=modulesC1)
ME_M1 = PCsM1$eigengenes
geneModuleMembershipM1 = signedKME(datExprM1, ME_M1)
colnames(geneModuleMembershipM1)=paste("PC",colorsC1,".cor",sep="");
MMPvalueM1=corPvalueStudent(as.matrix(geneModuleMembershipM1),dim(t(datExprM1))[[2]]);
colnames(MMPvalueM1)=paste("PC",colorsC1,".pval",sep="");
kMEtableM1 = cbind(Gene,Gene,modulesC1)
for (i in 1:length(colorsC1))
  kMEtableM1 = cbind(kMEtableM1, geneModuleMembershipM1[,i], MMPvalueM1[,i])
colnames(kMEtableM1)=colnames(kMEtableM1)
write.csv(kMEtableM1,"kMEtableM1.csv",row.names=FALSE)


### several ways to put Kme into practice #####
####The first thing we can do is plot the kME values

pdf("all_kMEtable2_vs_kMEtable1.pdf",height=8,width=8)

for (c in 1:length(colorsC1)){
  verboseScatterplot(geneModuleMembershipM1[,c],geneModuleMembershipC1[,c],main=colorsC1[c],
                     xlab="kME in M1",ylab="kME in C1")
}

dev.off()

pdf("inModule_kMEtable2_vs_kMEtable1.pdf",height=8,width=8)

###plots for screening ###
plot(geneModuleMembershipM3[,1],geneModuleMembershipC3[,1],main= "diff kME for screening"???
       col = colorsC3[1], xlim = c(-1,1), ylim = c(-1,1),
     xlab="kME in M3",ylab="kME in C3")

for (c in 2:length(colorsC3)){
  points(geneModuleMembershipM3[,c],geneModuleMembershipC3[,c],col = colorsC3[c],
         xlab="kME in M3",ylab="kME in C3")
}


for (c in 1:length(colorsC1)){
  inMod = modulesC1== colorsC1[c]
  verboseScatterplot(geneModuleMembershipM1[inMod,c],geneModuleMembershipC1[inMod,c],main=colorsC1[c],
                     col = colorsC1[c], add = T,
                     xlab="kME in M1",ylab="kME in C1")
}

dev.off()

#save.image("tutorial.RData") #(optional line of code)
#(Similar

Gene_C1 = colnames(datExprC1)
###The second thing we can do is determine which genes are hubs in both networks####
### These genes represent the top 10 genes per module based on kME in both networks.
topGenesKME = NULL
for (c in 1:length(colorsC1)){
  kMErank1 = rank(-geneModuleMembershipC1[,c])
  kMErank2 = rank(-geneModuleMembershipM1[,c])
  maxKMErank = rank(apply(cbind(kMErank1,kMErank2+.00001),1,max))
  topGenesKME = cbind(topGenesKME,Gene[maxKMErank<=10])
}; colnames(topGenesKME) = colorsC1
topGenesKME
### These genes represent the top 10 genes per module based on kME in both networks.
mostkME_3 = NULL
for (c in 1:length(colorsC3)){
  kMErank1 = rank(-geneModuleMembershipC3[,c])
  kMErank2 = rank(-geneModuleMembershipM3[,c])
  dif = kMErank2-kMErank1
  mostkME_3 = cbind(mostkME_3,dif) 
  
  #    maxKMErank = rank(apply(cbind(kMErank1,kMErank2+.00001),1,max))
  #  topGenesKME = cbind(topGenesKME,Gene[maxKMErank<=10])
}; 

change = cbind (Gene_C3, mostkME_3); colnames(change) = c("id",colorsC2)

change = data.frame( ID = change[,1], Gene_name = networkName[c(350:418),2], change[,c(2:3)] )

##example sorting genes in weighted##


order(abs(as.numeric(change[,3])))
##### sort data#####
rank_change = NULL
for (c in c(3:4)){
  changefinal = change[order(abs(as.numeric(change[,c])), decreasing = T),c(1:2,c)]
  write.csv(changefinal,paste("changefinal2_",c,sep = "_"))
}  

#####Third thing is find top genes in each module and compare their changes across condition  ###
###  i.e. top 10 in C3 change to what in M3
topGenesKMEC1_each = NULL
for( b in 1:length(colorsC1)){
  topGenesKMEC1 = NULL
  geneModuleMembershipC1 = geneModuleMembershipC1[order(geneModuleMembershipC1[,b],decreasing = T),]
  geneModuleMembershipM1 = geneModuleMembershipM1[order(geneModuleMembershipM1[,b],decreasing = T),]
  TOP1 = data.frame(geneModuleMembershipC1[c(1:10),b],
                    c(1:10),
                    geneModuleMembershipM1[rownames(geneModuleMembershipC1)[c(1:10)],b],
                    match(rownames(geneModuleMembershipC1)[c(1:10)],rownames(geneModuleMembershipM1)),
                    row.names = rownames(geneModuleMembershipC1)[c(1:10)])
  TOP2 = data.frame(geneModuleMembershipC1[c((length(Gene_C1)-9):length(Gene_C1)),b],
                    c((length(Gene_C1)-9):length(Gene_C1)),
                    geneModuleMembershipM1[rownames(geneModuleMembershipC1)[c((length(Gene_C1)-9):length(Gene_C1))],b],
                    match(rownames(geneModuleMembershipC1)[c((length(Gene_C1)-9):length(Gene_C1))],rownames(geneModuleMembershipM1)),
                    row.names = rownames(geneModuleMembershipC1)[c(1:10)])
  names(TOP1) = names(TOP2)= c(paste(colorsC1[b],"C3"),"inmd_rank1",paste(colorsC1[b],"M3"),"inmd_rank2")
  topGenesKMEC1 = rbind(TOP1,TOP2)
  topGenesKMEC1_each = list(topGenesKMEC1_each,topGenesKMEC1)
}

topGenesKMEC1_each 



##### fourth thing we can do is find genes with relatively high change  #####

### question : rank or value  ???

class(geneModuleMembershipC1)
head(row.names(geneModuleMembershipM1))

#####
ME_change_C1 = data.frame(geneModuleMembershipC1-geneModuleMembershipM1, row.names = row.names(geneModuleMembershipM1))



