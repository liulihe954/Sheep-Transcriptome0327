install.packages("xlsx")

library(xlsx)
install.packages("rJave")


dev.off()

sort(Results_sperm_cor$change$con_change)
Top50_index_order = order(-abs(Results_sperm_cor$change$con_change))[1:32]
subtest_exprGT = datExprGT[,Top50_index_order]
subtest_exprGO = datExprGO[,Top50_index_order]
test_results = get_NC_cor(subtest_exprGT,subtest_exprGO,0.5,0.05)

Top50_index1 = allID$Ensembl_ID[Top50_index_order]
Top50_index2 = allID$Gene_Names[Top50_index_order]


rownames(test_results$change) == rownames(test_results$change)

display_con1 =test_results$change[,(1:4)][order(-test_results$change[,1]),];display_con1
display_con2 =test_results$change[,(1:4)][order(-test_results$change[,3]),];display_con2


Top50_index2[order(-test_results$change[,1])]


cbind(display_con1,Top50_index2[order(-test_results$change[,1])],seq(32,1,by=-1))
cbind(display_con2,Top50_index2[order(-test_results$change[,3])],seq(32,1,by=-1))


order(test_results$change[,1])

dim(datExprGT)

Results_sperm_cor$change$con_2


Results_sperm_cor$change$con_change[order(-abs(Results_sperm_cor$change$con_change))]


install.packages("xlsxjars")
install.packages("rJava")
install.packages("xlsx")

library(grDevices)
library(utils)
library(rJava)
library(xlsxjars)
library(xlsx)
require(rJava)

install.packages('rJava', type='source')

library(rJava)



