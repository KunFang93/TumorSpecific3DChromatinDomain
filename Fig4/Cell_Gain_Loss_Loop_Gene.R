library(dplyr)
library(DESeq2)
library(pheatmap)
library(readxl)
library(writexl)
library(UpSetR)

data = read.table('/Users/KunFang/Documents/lab/Jin_lab/Labmates/lava/BRCA_tissue_Analysis/SLG-DLG/fithic/results_res20k_083121/analysis/GC_loops/cell_gene_normc_mat.txt',
                  header = T, row.names = "Gene")

# identify DLG in MCF7, MCF7TR, T47D, T47DTR, ZR75-1 and ZR75-1TR
data$`MCF7TRsubP` = data$MCF7TR - data$MCF7
data$`T47DTRsubP` = data$T47DTR - data$T47D
data$`ZR75.1TRsubP` = data$ZR75.1TR - data$ZR75.1

data$`MCF7TRdivP` = (data$MCF7TR+1)/(data$MCF7+1)
data$`T47DTRdivP` = (data$T47DTR+1)/(data$T47D+1)
data$`ZR75.1TRdivP` = (data$ZR75.1TR+1)/(data$ZR75.1+1)

data_subdiv_scaled = data.frame(scale(data[,c('MCF7TRsubP','T47DTRsubP','ZR75.1TRsubP',"MCF7TRdivP",'T47DTRdivP','ZR75.1TRdivP')]))
# remove outliers q3 + IQR
quantilecut = apply(data_subdiv_scaled, 2 , quantile , probs = c(0.15,0.2,0.25,0.75,0.8,0.85) , na.rm = TRUE )
IQR =  quantilecut["75%",]-quantilecut["25%",]
quantilecut = rbind("lowout" = quantilecut["25%",] - 1.5*IQR,quantilecut)
quantilecut = rbind("highout" = quantilecut["75%",] + 1.5*IQR,quantilecut)

# DLGs identification
# MCF7TRvsMCF7
MCF7TRP = data_subdiv_scaled[,c('MCF7TRsubP','MCF7TRdivP')]
MCF7TRP = MCF7TRP[MCF7TRP[,'MCF7TRsubP']>=quantilecut["lowout",'MCF7TRsubP']&
                    MCF7TRP[,'MCF7TRdivP']>=quantilecut["lowout",'MCF7TRdivP']&
                    MCF7TRP[,'MCF7TRsubP']<=quantilecut["highout",'MCF7TRsubP']&
                    MCF7TRP[,'MCF7TRdivP']<=quantilecut["highout",'MCF7TRdivP'],]
MCF7TR_Loss = subset(data_subdiv_scaled[,c('MCF7TRsubP','MCF7TRdivP')], MCF7TRsubP<= quantilecut["15%","MCF7TRsubP"] | MCF7TRdivP<= quantilecut["15%","MCF7TRdivP"])
MCF7TR_Gain = subset(data_subdiv_scaled[,c('MCF7TRsubP','MCF7TRdivP')], MCF7TRsubP>= quantilecut["85%","MCF7TRsubP"] | MCF7TRdivP>= quantilecut["85%","MCF7TRdivP"])
par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
with(MCF7TRP, plot(MCF7TRsubP, MCF7TRdivP, pch=20, main="MCF7TRvsMCF7", cex=1.0, xlab="MCF7TR-MCF7", ylab="MCF7TR/MCF7"))
with(subset(MCF7TRP, MCF7TRsubP<= quantilecut["15%","MCF7TRsubP"] | MCF7TRdivP<= quantilecut["15%","MCF7TRdivP"]), points(MCF7TRsubP,MCF7TRdivP, pch=20, col="green", cex=1.0))
with(subset(MCF7TRP, MCF7TRsubP>= quantilecut["85%","MCF7TRsubP"] | MCF7TRdivP>= quantilecut["85%","MCF7TRdivP"]), points(MCF7TRsubP,MCF7TRdivP, pch=20, col="red", cex=1.0))

# T47DTRvsT47D
T47DTRP = data_subdiv_scaled[c('T47DTRsubP','T47DTRdivP')]
T47DTRP = T47DTRP[T47DTRP[,'T47DTRsubP']>=quantilecut["lowout",'T47DTRsubP']&
                  T47DTRP[,'T47DTRdivP']>=quantilecut["lowout",'T47DTRdivP']&
                  T47DTRP[,'T47DTRsubP']<=quantilecut["highout",'T47DTRsubP']&
                  T47DTRP[,'T47DTRdivP']<=quantilecut["highout",'T47DTRdivP'],]
T47DTR_Loss = subset(data_subdiv_scaled[c('T47DTRsubP','T47DTRdivP')], T47DTRsubP<= quantilecut["15%","T47DTRsubP"] | T47DTRdivP<= quantilecut["15%","T47DTRdivP"])
T47DTR_Gain = subset(data_subdiv_scaled[c('T47DTRsubP','T47DTRdivP')], T47DTRsubP>= quantilecut["85%","T47DTRsubP"] | T47DTRdivP>= quantilecut["85%","T47DTRdivP"])
par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
with(T47DTRP, plot(T47DTRsubP, T47DTRdivP, pch=20, main="T47DTRvsT47D", cex=1.0, xlab="T47DTR-T47D", ylab="T47DTR/T47D"))
with(subset(T47DTRP, T47DTRsubP<= quantilecut["15%","T47DTRsubP"] | T47DTRdivP<= quantilecut["15%","T47DTRdivP"]), points(T47DTRsubP,T47DTRdivP, pch=20, col="green", cex=1.0))
with(subset(T47DTRP, T47DTRsubP>= quantilecut["85%","T47DTRsubP"] | T47DTRdivP>= quantilecut["85%","T47DTRdivP"]), points(T47DTRsubP,T47DTRdivP, pch=20, col="red", cex=1.0))

# ZR75.1TRvsZR75.1
ZR75.1TRP = data_subdiv_scaled[c('ZR75.1TRsubP','ZR75.1TRdivP')]
ZR75.1TRP = ZR75.1TRP[ZR75.1TRP[,'ZR75.1TRsubP']>=quantilecut["lowout",'ZR75.1TRsubP']&
                        ZR75.1TRP[,'ZR75.1TRdivP']>=quantilecut["lowout",'ZR75.1TRdivP']&
                        ZR75.1TRP[,'ZR75.1TRsubP']<=quantilecut["highout",'ZR75.1TRsubP']&
                        ZR75.1TRP[,'ZR75.1TRdivP']<=quantilecut["highout",'ZR75.1TRdivP'],]
ZR75.1TR_Loss = subset(data_subdiv_scaled[c('ZR75.1TRsubP','ZR75.1TRdivP')], ZR75.1TRsubP<= quantilecut["15%","ZR75.1TRsubP"] | ZR75.1TRdivP<= quantilecut["15%","ZR75.1TRdivP"])
ZR75.1TR_Gain = subset(data_subdiv_scaled[c('ZR75.1TRsubP','ZR75.1TRdivP')], ZR75.1TRsubP>= quantilecut["85%","ZR75.1TRsubP"] | ZR75.1TRdivP>= quantilecut["85%","ZR75.1TRdivP"])
par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
with(ZR75.1TRP, plot(ZR75.1TRsubP, ZR75.1TRdivP, pch=20, main="ZR75.1TRvsZR75.1", cex=1.0, xlab="ZR75.1TR-ZR75.1", ylab="ZR75.1TR/ZR75.1"))
with(subset(ZR75.1TRP, ZR75.1TRsubP<= quantilecut["15%","ZR75.1TRsubP"] | ZR75.1TRdivP<= quantilecut["15%","ZR75.1TRdivP"]), points(ZR75.1TRsubP,ZR75.1TRdivP, pch=20, col="green", cex=1.0))
with(subset(ZR75.1TRP, ZR75.1TRsubP>= quantilecut["85%","ZR75.1TRsubP"] | ZR75.1TRdivP>= quantilecut["85%","ZR75.1TRdivP"]), points(ZR75.1TRsubP,ZR75.1TRdivP, pch=20, col="red", cex=1.0))

MCF7TR_Gain = cbind("Genes"=rownames(MCF7TR_Gain), data[rownames(MCF7TR_Gain),c("MCF7","MCF7TR",colnames(MCF7TR_Gain))])
MCF7TR_Loss = cbind("Genes"=rownames(MCF7TR_Loss), data[rownames(MCF7TR_Loss),c("MCF7","MCF7TR",colnames(MCF7TR_Loss))])
sheets_mcf7 <- list("MCF7TR_Gain"=MCF7TR_Gain,"MCF7TR_Loss"=MCF7TR_Loss)
write_xlsx(sheets_mcf7,
           '/Users/KunFang/Documents/lab/Jin_lab/Labmates/lava/BRCA_tissue_Analysis/SLG-DLG/fithic/results_res20k_083121/analysis/Individual_specific_loops/MCF7TR_DLs.xlsx')

T47DTR_Gain = cbind("Genes"=rownames(T47DTR_Gain), data[rownames(T47DTR_Gain),c("T47D","T47DTR",colnames(T47DTR_Gain))])
T47DTR_Loss = cbind("Genes"=rownames(T47DTR_Loss), data[rownames(T47DTR_Loss),c("T47D","T47DTR",colnames(T47DTR_Loss))])
sheets_t47d <- list("T47DTR_Gain"=T47DTR_Gain,"T47DTR_Loss"=T47DTR_Loss)
write_xlsx(sheets_t47d,
           '/Users/KunFang/Documents/lab/Jin_lab/Labmates/lava/BRCA_tissue_Analysis/SLG-DLG/fithic/results_res20k_083121/analysis/Individual_specific_loops/T47DTR_DLs.xlsx')


ZR75.1TR_Gain = cbind("Genes"=rownames(ZR75.1TR_Gain), data[rownames(ZR75.1TR_Gain),c("ZR75.1","ZR75.1TR",colnames(ZR75.1TR_Gain))])
ZR75.1TR_Loss = cbind("Genes"=rownames(ZR75.1TR_Loss), data[rownames(ZR75.1TR_Loss),c("ZR75.1","ZR75.1TR",colnames(ZR75.1TR_Loss))])
sheets_zr751 <- list("ZR75.1TR_Gain"=ZR75.1TR_Gain,"ZR75.1TR_Loss"=ZR75.1TR_Loss)
write_xlsx(sheets_zr751,
           '/Users/KunFang/Documents/lab/Jin_lab/Labmates/lava/BRCA_tissue_Analysis/SLG-DLG/fithic/results_res20k_083121/analysis/Individual_specific_loops/ZR75.1TR_DLs.xlsx')

# construct binary matrix
MCF7TR_Gain_bi = data.frame("Gene"=rownames(MCF7TR_Gain),"MCF7TR"=rep(1,dim(MCF7TR_Gain)[1]))
T47DTR_Gain_bi = data.frame("Gene"=rownames(T47DTR_Gain),"T47DTR"=rep(1,dim(T47DTR_Gain)[1]))
ZR75.1TR_Gain_bi = data.frame("Gene"=rownames(ZR75.1TR_Gain),"ZR75.1TR"=rep(1,dim(ZR75.1TR_Gain)[1]))

MCF7TR_Loss_bi = data.frame("Gene"=rownames(MCF7TR_Loss),"MCF7TR"=rep(1,dim(MCF7TR_Loss)[1]))
T47DTR_Loss_bi = data.frame("Gene"=rownames(T47DTR_Loss),"T47DTR"=rep(1,dim(T47DTR_Loss)[1]))
ZR75.1TR_Loss_bi = data.frame("Gene"=rownames(ZR75.1TR_Loss),"ZR75.1TR"=rep(1,dim(ZR75.1TR_Loss)[1]))

binary_Gain = Reduce(function(x, y) merge(x, y, by="Gene", all=TRUE), list(MCF7TR_Gain_bi, T47DTR_Gain_bi, ZR75.1TR_Gain_bi))
binary_Gain[is.na(binary_Gain)] <- 0
binary_Loss = Reduce(function(x, y) merge(x, y, by="Gene", all=TRUE), list(MCF7TR_Loss_bi, T47DTR_Loss_bi, ZR75.1TR_Loss_bi))
binary_Loss[is.na(binary_Loss)] <- 0
rownames(binary_Gain) = binary_Gain$Gene
rownames(binary_Loss) = binary_Loss$Gene

binary_Gain = binary_Gain[,c("MCF7TR","T47DTR","ZR75.1TR")]
binary_Loss = binary_Loss[,c("MCF7TR","T47DTR","ZR75.1TR")]

MCF7TR_specific_Gain = rownames(binary_Gain[rowSums(binary_Gain)==1&binary_Gain$MCF7TR==1,])
T47DTR_specific_Gain = rownames(binary_Gain[rowSums(binary_Gain)==1&binary_Gain$T47DTR==1,])
ZR75.1TR_specific_Gain = rownames(binary_Gain[rowSums(binary_Gain)==1&binary_Gain$ZR75.1TR==1,])

MCF7TR_specific_Loss = rownames(binary_Loss[rowSums(binary_Loss)==1&binary_Loss$MCF7TR==1,])
T47DTR_specific_Loss = rownames(binary_Loss[rowSums(binary_Loss)==1&binary_Loss$T47DTR==1,])
ZR75.1TR_specific_Loss = rownames(binary_Loss[rowSums(binary_Loss)==1&binary_Loss$ZR75.1TR==1,])

MCF7TR_specific_Gain_signal = cbind("Genes"=MCF7TR_specific_Gain,data[MCF7TR_specific_Gain,])
T47DTR_specific_Gain_signal = cbind("Genes"=T47DTR_specific_Gain,data[T47DTR_specific_Gain,])
ZR75.1TR_specific_Gain_signal = cbind("Genes"=ZR75.1TR_specific_Gain,data[ZR75.1TR_specific_Gain,])
MCF7TR_specific_Loss_signal = cbind("Genes"=MCF7TR_specific_Loss,data[MCF7TR_specific_Loss,])
T47DTR_specific_Loss_signal = cbind("Genes"=T47DTR_specific_Loss,data[T47DTR_specific_Loss,])
ZR75.1TR_specific_Loss_signal = cbind("Genes"=ZR75.1TR_specific_Loss,data[ZR75.1TR_specific_Loss,])

sheets_specific <- list("MCF7TR_Gain"=MCF7TR_specific_Gain_signal,"T47DTR_Gain"=T47DTR_specific_Gain_signal,
                        "ZR75-1TR_Gain"=ZR75.1TR_specific_Gain_signal,"MCF7TR_Loss"=MCF7TR_specific_Loss_signal,
                        "T47DTR_Loss"=T47DTR_specific_Loss_signal,"ZR75-1TR_Loss"=ZR75.1TR_specific_Loss_signal)
write_xlsx(sheets_specific,'/Users/KunFang/Documents/lab/Jin_lab/Labmates/lava/BRCA_tissue_Analysis/SLG-DLG/fithic/results_res20k_083121/analysis/Individual_specific_loops/cell_TR_specific_DLGs.xlsx')


categories = c("MCF7TR","T47DTR","ZR75.1TR")
loop_metadata = as.data.frame(cbind("sample"=categories,categories))
rownames(loop_metadata)= loop_metadata$categories
upset(binary_Gain, 
      sets = c("ZR75.1TR","T47DTR","MCF7TR"),
      keep.order = T, matrix.color="black", point.size=2, sets.x.label = 'Total genes number',
      mainbar.y.label = 'number of genes', text.scale = 1.3,
      sets.bar.color=c("black","black","black"),
      set.metadata = list(data = loop_metadata,
                          plots = list(list(type = "matrix_rows", column = "categories",
                                            colors = c(MCF7TR = "wheat", T47DTR = "navy", ZR75.1TR = "purple"),
                                            alpha = 0.5))),
      queries = list(list(query = intersects,params = list("MCF7TR"), color = "red", active = T),
                     list(query = intersects,params = list("T47DTR"), color = "red", active = T),
                     list(query = intersects,params = list("ZR75.1TR"), color = "red", active = T)),
)

upset(binary_Loss, 
      sets = c("ZR75.1TR","T47DTR","MCF7TR"),
      keep.order = T, matrix.color="black", point.size=2, sets.x.label = 'Total genes number',
      mainbar.y.label = 'number of genes', text.scale = 1.3,
      sets.bar.color=c("black","black","black"),
      set.metadata = list(data = loop_metadata,
                          plots = list(list(type = "matrix_rows", column = "categories",
                                            colors = c(MCF7TR = "wheat", T47DTR = "navy", ZR75.1TR = "purple"),
                                            alpha = 0.5))),
      queries = list(list(query = intersects,params = list("MCF7TR"), color = "green", active = T),
                     list(query = intersects,params = list("T47DTR"), color = "green", active = T),
                     list(query = intersects,params = list("ZR75.1TR"), color = "green", active = T)),
)

binary_total = rbind(binary_Gain,binary_Loss)
upset(binary_total, 
      sets = c("ZR75.1TR","T47DTR","MCF7TR"),
      keep.order = T, matrix.color="black", point.size=2, sets.x.label = 'Total genes number',
      mainbar.y.label = 'number of genes', text.scale = 1.3,
      sets.bar.color=c("black","black","black"),
      set.metadata = list(data = loop_metadata,
                          plots = list(list(type = "matrix_rows", column = "categories",
                                            colors = c(MCF7TR = "wheat", T47DTR = "navy", ZR75.1TR = "purple"),
                                            alpha = 0.5))),
      queries = list(list(query = intersects,params = list("MCF7TR"), color = "orange", active = T),
                     list(query = intersects,params = list("T47DTR"), color = "orange", active = T),
                     list(query = intersects,params = list("ZR75.1TR"), color = "orange", active = T)),
)




