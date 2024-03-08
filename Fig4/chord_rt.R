library(tidyverse)
library(viridis)
library(patchwork)
library(hrbrthemes)
library(circlize)
library(chorddiag)
library(readxl)
library(Corbi)
library(UpSetR)
library(stringr)
library(writexl)
library(dendextend)
# Loading data
datdir = '/Users/KunFang/Documents/lab/Jin_lab/Labmates/lava/BRCA_tissue_Analysis/SLG-DLG/fithic/results_res20k_083121/analysis/Individual_specific_loops/'
rt_samples = c('RT1','RT2','RT3','RT4','RT5')
cell_samples = c('MCF7TR','T47DTR','ZR75.1TR')
column_names = c('Genes','PT1','PT2','PT3','PT4','PT5','RT1','RT2','RT3','RT4','RT5')
read_file <- function(sheet_name,colnames){
  df = read_excel(paste(datdir,'RTsIndvsPT.xlsx',sep='/'),sheet=sheet_name,col_names = column_names,skip = 1)
  return(df)
}

RTsIndvsPT_list = list()
for (sample in rt_samples){
  RTsIndvsPT_list[[paste(sample,"Gained",sep="_")]] = read_file(paste(sample,"Gained",sep="_"),column_names)
  RTsIndvsPT_list[[paste(sample,"Lost",sep="_")]] = read_file(paste(sample,"Lost",sep="_"),column_names)
}


cellTR_list = list()
for (sample in cell_samples){
  cellTR_list[[paste(sample,"Gained",sep="_")]] = read_excel(paste(datdir,paste(sample,'DLs.xlsx',sep='_'),sep='/'),sheet=paste(sample,"Gain",sep='_'))
  cellTR_list[[paste(sample,"Lost",sep="_")]] = read_excel(paste(datdir,paste(sample,'DLs.xlsx',sep='_'),sep='/'),sheet=paste(sample,"Loss",sep='_'))
}

binary_Gained_list = list()
for (sample in c(rt_samples,cell_samples)){
  current_key = paste(sample,"Gained",sep="_")
  if (sample %in% rt_samples){
    current_df = RTsIndvsPT_list[[current_key]]
  } else {
    current_df = cellTR_list[[current_key]]
  }
  current_bi = data.frame("Gene"=current_df$Genes,"V1"=rep(1,dim(current_df)[1]))
  colnames(current_bi) = c("Gene",sample)
  binary_Gained_list[[current_key]] = current_bi
}

all_samples = c("RT1","RT2","RT3","RT4","RT5","MCF7TR","T47DTR","ZR75.1TR")
binary_Gain = Reduce(function(x, y) merge(x, y, by="Gene", all=TRUE), 
                     binary_Gained_list)
binary_Gain[is.na(binary_Gain)] <- 0
rownames(binary_Gain)=binary_Gain$Gene
binary_Gain = binary_Gain[,all_samples]

common_gain_list = c()
names_gain_list = c()
adj_gain_mat = data.frame(matrix(ncol=8,nrow=8, dimnames=list(all_samples, all_samples)))
for (sample_x in all_samples){
  for (sample_y in all_samples){
    current_name = paste(sample_x,sample_y,sep='_')
    names_gain_list = c(names_gain_list,current_name)
    if (sample_x == sample_y){
      current_df = binary_Gain[binary_Gain[,sample_x]==1 & rowSums(binary_Gain)==1,]
      common_gain_list = c(common_gain_list,current_df)
      adj_gain_mat[sample_x,sample_y] = dim(current_df)[1]
    } else {
      current_df = binary_Gain[rowSums(binary_Gain[,c(sample_x,sample_y)])==2,]
      common_gain_list = c(common_gain_list,current_df)
      adj_gain_mat[sample_x,sample_y] = dim(current_df)[1]
    }
  }
}

binary_Lost_list = list()
for (sample in c(rt_samples,cell_samples)){
  current_key = paste(sample,"Lost",sep="_")
  if (sample %in% rt_samples){
    current_df = RTsIndvsPT_list[[current_key]]
  } else {
    current_df = cellTR_list[[current_key]]
  }
  current_bi = data.frame("Gene"=current_df$Genes,"V1"=rep(1,dim(current_df)[1]))
  colnames(current_bi) = c("Gene",sample)
  binary_Lost_list[[current_key]] = current_bi
}

binary_Loss = Reduce(function(x, y) merge(x, y, by="Gene", all=TRUE), 
                     binary_Lost_list)
binary_Loss[is.na(binary_Loss)] <- 0
rownames(binary_Loss)=binary_Loss$Gene
binary_Loss = binary_Loss[,all_samples]

common_loss_list = c()
names_loss_list = c()
adj_loss_mat = data.frame(matrix(ncol=length(all_samples),nrow=length(all_samples), dimnames=list(all_samples, all_samples)))
for (sample_x in all_samples){
  for (sample_y in all_samples){
    current_name = paste(sample_x,sample_y,sep='_')
    names_loss_list = c(names_loss_list,current_name)
    if (sample_x == sample_y){
      current_df = binary_Loss[binary_Loss[,sample_x]==1 & rowSums(binary_Loss)==1,]
      common_loss_list = c(common_loss_list,current_df)
      adj_loss_mat[sample_x,sample_y] = dim(current_df)[1]
    } else {
      current_df = binary_Loss[rowSums(binary_Loss[,c(sample_x,sample_y)])==2,]
      common_loss_list = c(common_loss_list,current_df)
      adj_loss_mat[sample_x,sample_y] = dim(current_df)[1]
    }
  }
}

adj2long <- function(adj_mat){
  adj_long_mat <- adj_mat %>%
    rownames_to_column %>%
    gather(key = 'key', value = 'value', -rowname)
  return(adj_long_mat)
}

adj_gain_long_mat <- adj2long(adj_gain_mat)
adj_loss_long_mat <- adj2long(adj_loss_mat)

adj_gain_long_mat <- adj_gain_long_mat[c(1:40,46,55,64),]
adj_loss_long_mat <- adj_loss_long_mat[c(1:40,46,55,64),]

# avoid plot twice
adj_gain_long_mat[adj_gain_long_mat$rowname==adj_gain_long_mat$key,'value']= adj_gain_long_mat[adj_gain_long_mat$rowname==adj_gain_long_mat$key,'value']/2
adj_loss_long_mat[adj_loss_long_mat$rowname==adj_loss_long_mat$key,'value']= adj_loss_long_mat[adj_loss_long_mat$rowname==adj_loss_long_mat$key,'value']/2

mycolors <- viridis(length(all_samples), alpha = 1, begin = 0, end = 1, option = "D")
names(mycolors) <- all_samples
# Base plot
plot_chord <- function(long_mat,samples_list,linkcolors){
  circos.clear()
  circos.par(start.degree = 90, gap.degree = 4,  track.margin = c(-0.1, 0.1),points.overflow.warning = FALSE)
  par(mar = rep(0, 4))
  # color palette
  gridcolors <- mycolors
  chordDiagram(
    x = long_mat, 
    grid.col = gridcolors,
    col = linkcolors,
    transparency = 0.25,
    diffHeight  = -0.04,
    directional = 1,
    direction.type = c("arrows", "diffHeight"),
    annotationTrack = "grid",
    annotationTrackHeight = c(0.05, 0.01),
    link.arr.type = "big.arrow",
    link.sort = TRUE,
    link.largest.ontop = TRUE)
  
  # Add text and axis
  circos.trackPlotRegion(
    track.index = 1,
    bg.border = NA,
    panel.fun = function(x, y) {

      xlim = get.cell.meta.data("xlim")
      sector.index = get.cell.meta.data("sector.index")

      # Add names to the sector.
      # circos.text(
      #   x = mean(xlim),
      #   y = 4,
      #   labels = sector.index,
      #   cex = 1
      # )

      # Add graduation on axis
      circos.axis(
        h = "top",
        major.tick = TRUE,
        # major.at = seq(from = 0, to = xlim[2], by = ifelse(test = xlim[2]>500, yes = 400, no = 200)),
        minor.ticks = 4,
        major.tick.length = 0.5,
        labels.niceFacing = FALSE,
        labels.cex = 0.8,
        labels.pos.adjust = FALSE
      )
    }
  )
}

adj_gain_long_mat$colors <- mycolors[adj_gain_long_mat$rowname]
adj_gain_long_mat[! ((adj_gain_long_mat$key==adj_gain_long_mat$rowname)|adj_gain_long_mat$rowname %in% cell_samples),"colors"]  = "#00000000"
adj_gain_long_mat[adj_gain_long_mat$key %in% cell_samples,"colors"] = "#00000000"
gain_linkcolors = adj_gain_long_mat$colors

adj_loss_long_mat$colors <- mycolors[adj_loss_long_mat$rowname]
adj_loss_long_mat[! ((adj_loss_long_mat$key==adj_gain_long_mat$rowname)|adj_loss_long_mat$rowname %in% cell_samples),"colors"]  = "#00000000"
adj_loss_long_mat[adj_loss_long_mat$key %in% cell_samples,"colors"] = "#00000000"
loss_linkcolors = adj_loss_long_mat$colors

# Figure3A Chord diagram 
plot_chord(adj_gain_long_mat,all_samples,gain_linkcolors)
plot_chord(adj_loss_long_mat,all_samples,loss_linkcolors)

# read Figure2D data
new_col = c("Genes","RT1vsPTs", "RT2vsPTs", "RT3vsPTs", "RT4vsPTs", "RT5vsPTs", "Commons")
RTsvsPTs_common_list = list()
RTsvsPTs_common_list[[">=2RTs_Gained"]] = read_xlsx(paste(datdir,"RTsIndvsPT_commons.xlsx",sep='/'),sheet = ">=2_RTs_Gained",col_names = new_col,skip = 1)
RTsvsPTs_common_list[[">=3RTs_Gained"]] = read_xlsx(paste(datdir,"RTsIndvsPT_commons.xlsx",sep='/'),sheet = ">=3_RTs_Gained",col_names = new_col,skip = 1)
RTsvsPTs_common_list[[">=4RTs_Gained"]] = read_xlsx(paste(datdir,"RTsIndvsPT_commons.xlsx",sep='/'),sheet = ">=4_RTs_Gained",col_names = new_col,skip = 1)

RTsvsPTs_common_list[[">=2RTs_Lost"]] = read_xlsx(paste(datdir,"RTsIndvsPT_commons.xlsx",sep='/'),sheet = ">=2_RTs_Lost",col_names = new_col,skip = 1)
RTsvsPTs_common_list[[">=3RTs_Lost"]] = read_xlsx(paste(datdir,"RTsIndvsPT_commons.xlsx",sep='/'),sheet = ">=3_RTs_Lost",col_names = new_col,skip = 1)
RTsvsPTs_common_list[[">=4RTs_Lost"]] = read_xlsx(paste(datdir,"RTsIndvsPT_commons.xlsx",sep='/'),sheet = ">=4_RTs_Lost",col_names = new_col,skip = 1)
RTsvsPTs_common_list[[">=5RTs_Lost"]] = read_xlsx(paste(datdir,"RTsIndvsPT_commons.xlsx",sep='/'),sheet = ">=5_RTs_Lost",col_names = new_col,skip = 1)

binary_RTTRcommon_Gain_list = list()
RTTR_gain_keys <- c(">=2RTs_Gained",">=3RTs_Gained",">=4RTs_Gained")
new_gain_cols = c()
for (key in RTTR_gain_keys){
  current_df = RTsvsPTs_common_list[[key]]
  current_bi = data.frame("Gene"=current_df$Genes,"V1"=rep(1,dim(current_df)[1]))
  colnames(current_bi) = c("Gene",str_split(key,"_")[[1]][1])
  new_gain_cols = c(new_gain_cols,str_split(key,"_")[[1]][1])
  binary_RTTRcommon_Gain_list[[key]] = current_bi
}
for (sample in cell_samples){
  current_df = data.frame(cbind("Gene"=rownames(binary_Gain[binary_Gain[sample]==1,]),
                                "V1"=binary_Gain[binary_Gain[sample]==1,sample]))
  colnames(current_df) = c("Gene",sample)
  binary_RTTRcommon_Gain_list[[sample]] = current_df
}
binary_RTTRcommon_Gain = Reduce(function(x, y) merge(x, y, by="Gene", all=TRUE), 
                    binary_RTTRcommon_Gain_list)
binary_RTTRcommon_Gain[is.na(binary_RTTRcommon_Gain)] <- 0
rownames(binary_RTTRcommon_Gain)=binary_RTTRcommon_Gain$Gene
binary_RTTRcommon_Gain = binary_RTTRcommon_Gain[,2:dim(binary_RTTRcommon_Gain)[2]]
binary_RTTRcommon_Gain[] <- lapply(binary_RTTRcommon_Gain, function(x) {
 as.numeric(as.character(x))
})

sheets_RTs_TRs_DLGs <- list()
RTs_TRs_gain_types = c()
RTs_TRs_gain_number = c()
for (key in new_gain_cols) {
  for (sample in cell_samples){
    current_df = binary_RTTRcommon_Gain[binary_RTTRcommon_Gain[key] ==1 & binary_RTTRcommon_Gain[sample] ==1,]
    current_key = paste(key,sample,sep="_")
    if (dim(current_df)[1] !=0){
      RTs_TRs_gain_types = c(RTs_TRs_gain_types,current_key)
      RTs_TRs_gain_number = c(RTs_TRs_gain_number,dim(current_df)[1])
      sheets_RTs_TRs_DLGs[[paste(current_key,"Gained",sep="_")]] = cbind("Genes"=rownames(current_df),current_df)
    }
  }
}

# upset
colors_upset = c("darkred","darkorange","gold3","green","blue","purple")
categories = c("ZR75.1TR","T47DTR","MCF7TR",">=4RTs",">=3RTs",">=2RTs")
names(colors_upset) = categories
gain_loop_metadata = as.data.frame(cbind("Type"=categories,categories))
upset(binary_RTTRcommon_Gain,
      sets = c("ZR75.1TR","T47DTR","MCF7TR",">=4RTs",">=3RTs",">=2RTs"),
      keep.order = T, matrix.color="black", point.size=2, sets.x.label = 'Total genes number',
      mainbar.y.label = 'number of genes', text.scale = 1.3,
      sets.bar.color = c("darkred","darkorange","gold3","green","blue","purple"),
      set.metadata = list(data = gain_loop_metadata,
                          plots = list(list(type = "matrix_rows", column = "categories",
                                            colors = c(ZR75.1TR = "darkred", T47DTR = "darkorange", MCF7TR = "gold3",
                                                       `>=4RTs`="green",`>=3RTs`="blue",
                                                       `>=2RTs`="purple"),
                                            alpha = 0.5))),
      queries = list(
        list(query = intersects, params = list(">=2RTs","MCF7TR"), color = "purple", active = T),
        list(query = intersects, params = list(">=2RTs","T47DTR"), color = "purple", active = T),
        list(query = intersects, params = list(">=2RTs","ZR75.1TR"), color = "purple", active = T),
        list(query = intersects, params = list(">=2RTs","MCF7TR","T47DTR"), color = "purple", active = T),
        list(query = intersects, params = list(">=2RTs","MCF7TR","ZR75.1TR"), color = "purple", active = T),
        list(query = intersects, params = list(">=2RTs","T47DTR","ZR75.1TR"), color = "purple", active = T),
        list(query = intersects, params = list(">=2RTs","MCF7TR","T47DTR","ZR75.1TR"), color = "purple", active = T),
        list(query = intersects, params = list(">=2RTs",">=3RTs","MCF7TR"), color = "blue", active = T),
        list(query = intersects, params = list(">=2RTs",">=3RTs","T47DTR"), color = "blue", active = T),
        list(query = intersects, params = list(">=2RTs",">=3RTs","ZR75.1TR"), color = "blue", active = T),
        list(query = intersects, params = list(">=2RTs",">=3RTs","MCF7TR","ZR75.1TR"), color = "blue", active = T),
        list(query = intersects, params = list(">=2RTs",">=3RTs","T47DTR","ZR75.1TR"), color = "blue", active = T),
        list(query = intersects, params = list(">=2RTs",">=3RTs","MCF7TR","T47DTR","ZR75.1TR"), color = "blue", active = T)
                    ),
)

# circular barplot
color_gain = rep(c("gold3","darkorange","darkred"),2)
cir_barplot_df = data.frame(cbind("Type"=RTs_TRs_gain_types,"Counts"=RTs_TRs_gain_number,"colors"=color_gain))
cir_barplot_df$Counts = as.numeric(cir_barplot_df$Counts)
cir_barplot_df = cir_barplot_df[order(cir_barplot_df$Counts),]
circos.clear()
circle_len = length(cir_barplot_df$Type)
circos.par("start.degree" = 90, cell.padding = c(0, 0, 0, 0))
circos.initialize("a", xlim = c(0, 250)) # 'a` just means there is one sector
circos.track(ylim = c(0.5, circle_len+0.5), track.height = 0.8, 
             bg.border = NA, panel.fun = function(x, y) {
               xlim = CELL_META$xlim
               circos.segments(rep(xlim[1], circle_len), 1:circle_len,
                               rep(xlim[2], circle_len), 1:circle_len,
                               col = "#CCCCCC")
               circos.rect(rep(0, circle_len), 1:circle_len - 0.45, cir_barplot_df$Counts, 1:circle_len + 0.45,
                           col = cir_barplot_df$colors, border = "white")
               circos.text(rep(xlim[1], circle_len), 1:circle_len, 
                           paste(cir_barplot_df$Type, ":", cir_barplot_df$Counts), 
                           facing = "inside", adj = c(1.05, 0.5), cex = 1.7) 
               breaks = seq(0, 220, by = 20)
               circos.axis(h = "top", major.at = breaks, labels = breaks, 
                           labels.cex = 1,labels.pos.adjust = FALSE)
             })


# Lost
binary_RTTRcommon_Loss_list = list()
RTTR_loss_keys <- c(">=2RTs_Lost",">=3RTs_Lost",">=4RTs_Lost",">=5RTs_Lost")
new_loss_cols = c()
for (key in RTTR_loss_keys){
  current_df = RTsvsPTs_common_list[[key]]
  current_bi = data.frame("Gene"=current_df$Genes,"V1"=rep(1,dim(current_df)[1]))
  colnames(current_bi) = c("Gene",str_split(key,"_")[[1]][1])
  new_loss_cols = c(new_loss_cols,str_split(key,"_")[[1]][1])
  binary_RTTRcommon_Loss_list[[key]] = current_bi
}
for (sample in cell_samples){
  current_df = data.frame(cbind("Gene"=rownames(binary_Loss[binary_Loss[sample]==1,]),
                                "V1"=binary_Loss[binary_Loss[sample]==1,sample]))
  colnames(current_df) = c("Gene",sample)
  binary_RTTRcommon_Loss_list[[sample]] = current_df
}
binary_RTTRcommon_Loss = Reduce(function(x, y) merge(x, y, by="Gene", all=TRUE), 
                                binary_RTTRcommon_Loss_list)
binary_RTTRcommon_Loss[is.na(binary_RTTRcommon_Loss)] <- 0
rownames(binary_RTTRcommon_Loss)=binary_RTTRcommon_Loss$Gene
binary_RTTRcommon_Loss = binary_RTTRcommon_Loss[,2:dim(binary_RTTRcommon_Loss)[2]]
binary_RTTRcommon_Loss[] <- lapply(binary_RTTRcommon_Loss, function(x) {
  as.numeric(as.character(x))
})

RTs_TRs_loss_types = c()
RTs_TRs_loss_number = c()
for (key in new_loss_cols) {
  for (sample in cell_samples){
    current_df = binary_RTTRcommon_Loss[binary_RTTRcommon_Loss[key] ==1 & binary_RTTRcommon_Loss[sample] ==1,]
    current_key = paste(key,sample,sep="_")
    if (dim(current_df)[1] !=0){
      RTs_TRs_loss_types = c(RTs_TRs_loss_types,current_key)
      RTs_TRs_loss_number = c(RTs_TRs_loss_number,dim(current_df)[1])
      sheets_RTs_TRs_DLGs[[paste(current_key,"Lost",sep="_")]] = cbind("Genes"=rownames(current_df),current_df)
    }
  }
}

# upset
loss_categories = c("ZR75.1TR","T47DTR","MCF7TR",">=5RTs",">=4RTs",">=3RTs",">=2RTs")
loss_loop_metadata = as.data.frame(cbind("Type"=loss_categories,"categories"=loss_categories))
upset(binary_RTTRcommon_Loss,
      sets = c("ZR75.1TR","T47DTR","MCF7TR",">=5RTs",">=4RTs",">=3RTs",">=2RTs"),
      keep.order = T, matrix.color="black", point.size=2, sets.x.label = 'Total genes number',
      mainbar.y.label = 'number of genes', text.scale = 1.3,
      sets.bar.color = c("darkred","darkorange","gold3","deeppink2","green","blue","purple"),
      set.metadata = list(data = loss_loop_metadata,
                          plots = list(list(type = "matrix_rows", column = "categories",
                                            colors = c(ZR75.1TR = "darkred", T47DTR = "darkorange", MCF7TR = "gold3",
                                                       `>=5RTs`="deeppink2",`>=4RTs`="green",`>=3RTs`="blue",`>=2RTs`="purple"),
                                            alpha = 0.5))),
      queries = list(
        list(query = intersects, params = list(">=2RTs","MCF7TR"), color = "purple", active = T),
        list(query = intersects, params = list(">=2RTs","T47DTR"), color = "purple", active = T),
        list(query = intersects, params = list(">=2RTs","ZR75.1TR"), color = "purple", active = T),
        list(query = intersects, params = list(">=2RTs","MCF7TR","T47DTR"), color = "purple", active = T),
        list(query = intersects, params = list(">=2RTs","MCF7TR","ZR75.1TR"), color = "purple", active = T),
        list(query = intersects, params = list(">=2RTs","T47DTR","ZR75.1TR"), color = "purple", active = T),
        list(query = intersects, params = list(">=2RTs","MCF7TR","T47DTR","ZR75.1TR"), color = "purple", active = T),
        list(query = intersects, params = list(">=2RTs",">=3RTs","MCF7TR"), color = "blue", active = T),
        list(query = intersects, params = list(">=2RTs",">=3RTs",">=4RTs","MCF7TR"), color = "green", active = T),
        list(query = intersects, params = list(">=2RTs",">=3RTs",">=4RTs","T47DTR"), color = "green", active = T),
        list(query = intersects, params = list(">=2RTs",">=3RTs","T47DTR"), color = "blue", active = T),
        list(query = intersects, params = list(">=2RTs",">=3RTs","ZR75.1TR"), color = "blue", active = T),
        list(query = intersects, params = list(">=2RTs",">=3RTs","MCF7TR","ZR75.1TR"), color = "blue", active = T),
        list(query = intersects, params = list(">=2RTs",">=3RTs","MCF7TR","T47DTR"), color = "blue", active = T),
        list(query = intersects, params = list(">=2RTs",">=3RTs",">=4RTs","MCF7TR","T47DTR"), color = "green", active = T),
        list(query = intersects, params = list(">=2RTs",">=3RTs","T47DTR","ZR75.1TR"), color = "blue", active = T),
        list(query = intersects, params = list(">=2RTs",">=3RTs","MCF7TR","T47DTR","ZR75.1TR"), color = "blue", active = T),
        list(query = intersects, params = list(">=2RTs",">=3RTs",">=4RTs",">=5RTs","ZR75.1TR"), color = "deeppink2", active = T)
      ),
)

# circular barplot
color_loss = c(rep(c("gold3","darkorange","darkred"),3),"darkred")
cir_barplot_df = data.frame(cbind("Type"=RTs_TRs_loss_types,"Counts"=RTs_TRs_loss_number,"colors"=color_loss))
cir_barplot_df$Counts = as.numeric(cir_barplot_df$Counts)
cir_barplot_df = cir_barplot_df[order(cir_barplot_df$Counts),]
circos.clear()
circle_len = length(cir_barplot_df$Type)
circos.par("start.degree" = 90, cell.padding = c(0, 0, 0, 0))
circos.initialize("a", xlim = c(0, 150)) # 'a` just means there is one sector
circos.track(ylim = c(0.5, circle_len+0.5), track.height = 0.8, 
             bg.border = NA, panel.fun = function(x, y) {
               xlim = CELL_META$xlim
               circos.segments(rep(xlim[1], circle_len), 1:circle_len,
                               rep(xlim[2], circle_len), 1:circle_len,
                               col = "#CCCCCC")
               circos.rect(rep(0, circle_len), 1:circle_len - 0.45, cir_barplot_df$Counts, 1:circle_len + 0.45,
                           col = cir_barplot_df$colors, border = "white")
               circos.text(rep(xlim[1], circle_len), 1:circle_len, 
                           paste(cir_barplot_df$Type, ":", cir_barplot_df$Counts), 
                           facing = "inside", adj = c(1.05, 0.5), cex = 1.7) 
               breaks = seq(0, 130, by = 20)
               circos.axis(h = "top", major.at = breaks, labels = breaks, 
                           labels.cex = 1,labels.pos.adjust = FALSE)
             })


write_xlsx(sheets_RTs_TRs_DLGs,paste(datdir,'RTs_TRs_DLGs.xlsx',sep='/'))

# loading RNA data
MCF7vsTR_UpR = read_excel('/Users/KunFang/Documents/lab/Jin_lab/Labmates/Yini/RNA-seq/xlsx_files/SupplfileS7_DEGs_MCF7vsMCF7TR.xlsx',sheet = 'Up-Regulated')
MCF7vsTR_DownR = read_excel('/Users/KunFang/Documents/lab/Jin_lab/Labmates/Yini/RNA-seq/xlsx_files/SupplfileS7_DEGs_MCF7vsMCF7TR.xlsx',sheet = 'Down-Regulated')
T47DvsTR_UpR = read_excel('/Users/KunFang/Documents/lab/Jin_lab/Labmates/Yini/RNA-seq/xlsx_files/SupplfileS9_DEGs_T47DvsT47DTR.xlsx',sheet = 'Up-Regulated')
T47DvsTR_DownR = read_excel('/Users/KunFang/Documents/lab/Jin_lab/Labmates/Yini/RNA-seq/xlsx_files/SupplfileS9_DEGs_T47DvsT47DTR.xlsx',sheet = 'Down-Regulated')

# RTs-TRs-Gained-UpR
sheets_RTs_TRs_DLGs_DEGs = list()
for ( i in new_gain_cols){
    sample = "MCF7TR"
    current_key = paste(i,sample,"Gained",sep="_")
    current_deg_df_select = as.data.frame(MCF7vsTR_UpR[MCF7vsTR_UpR$symbol %in% rownames(sheets_RTs_TRs_DLGs[[current_key]]),])
    current_dlg_df = as.data.frame(cellTR_list[paste(sample,"Gained",sep="_")])
    current_dlg_gene_col = paste(paste(sample,"Gained",sep="_"),"Genes",sep='.')
    current_dlg_df_select = current_dlg_df[current_dlg_df[,current_dlg_gene_col] %in% current_deg_df_select$symbol, ]
    current_merge = merge(current_deg_df_select,current_dlg_df_select,by.x='symbol',by.y=current_dlg_gene_col,all=TRUE)
    print(paste(i,dim(current_df)[1],sep=" "))
    if (dim(current_merge)[1]!=0){
      out_key = paste(i,sample,"Gained","UpR",sep="_")
      sheets_RTs_TRs_DLGs_DEGs[[out_key]] = current_merge
    }
}
  
for ( i in new_gain_cols){
  sample = "T47DTR"
  current_key = paste(i,sample,"Gained",sep="_")
  current_deg_df_select = as.data.frame(T47DvsTR_UpR[T47DvsTR_UpR$symbol %in% rownames(sheets_RTs_TRs_DLGs[[current_key]]),])
  current_dlg_df = as.data.frame(cellTR_list[paste(sample,"Gained",sep="_")])
  current_dlg_gene_col = paste(paste(sample,"Gained",sep="_"),"Genes",sep='.')
  current_dlg_df_select = current_dlg_df[current_dlg_df[,current_dlg_gene_col] %in% current_deg_df_select$symbol, ]
  current_merge = merge(current_deg_df_select,current_dlg_df_select,by.x='symbol',by.y=current_dlg_gene_col,all=TRUE)
  print(paste(i,dim(current_df)[1],sep=" "))
  if (dim(current_merge)[1]!=0){
    out_key = paste(i,sample,"Gained","UpR",sep="_")
    sheets_RTs_TRs_DLGs_DEGs[[out_key]] = current_merge
  }
}

# RTs-TRs-Lost-DownR
for ( i in new_loss_cols){
  sample = "MCF7TR"
  current_key = paste(i,sample,"Lost",sep="_")
  current_deg_df_select = as.data.frame(MCF7vsTR_DownR[MCF7vsTR_DownR$symbol %in% rownames(sheets_RTs_TRs_DLGs[[current_key]]),])
  current_dlg_df = as.data.frame(cellTR_list[paste(sample,"Lost",sep="_")])
  current_dlg_gene_col = paste(paste(sample,"Lost",sep="_"),"Genes",sep='.')
  current_dlg_df_select = current_dlg_df[current_dlg_df[,current_dlg_gene_col] %in% current_deg_df_select$symbol, ]
  current_merge = merge(current_deg_df_select,current_dlg_df_select,by.x='symbol',by.y=current_dlg_gene_col,all=TRUE)
  print(paste(i,dim(current_df)[1],sep=" "))
  if (dim(current_merge)[1]!=0){
    out_key = paste(i,sample,"Lost","DownR",sep="_")
    sheets_RTs_TRs_DLGs_DEGs[[out_key]] = current_merge
  }
}

for ( i in new_loss_cols){
  sample = "T47DTR"
  current_key = paste(i,sample,"Lost",sep="_")
  current_deg_df_select = as.data.frame(T47DvsTR_DownR[T47DvsTR_DownR$symbol %in% rownames(sheets_RTs_TRs_DLGs[[current_key]]),])
  current_dlg_df = as.data.frame(cellTR_list[paste(sample,"Lost",sep="_")])
  current_dlg_gene_col = paste(paste(sample,"Lost",sep="_"),"Genes",sep='.')
  current_dlg_df_select = current_dlg_df[current_dlg_df[,current_dlg_gene_col] %in% current_deg_df_select$symbol, ]
  current_merge = merge(current_deg_df_select,current_dlg_df_select,by.x='symbol',by.y=current_dlg_gene_col,all=TRUE)
  print(paste(i,dim(current_df)[1],sep=" "))
  if (dim(current_merge)[1]!=0){
    out_key = paste(i,sample,"Lost","DownR",sep="_")
    sheets_RTs_TRs_DLGs_DEGs[[out_key]] = current_merge
  }
}

# RTs-TRs-Gain-DownR
for ( i in new_gain_cols){
  sample = "MCF7TR"
  current_key = paste(i,sample,"Gained",sep="_")
  current_deg_df_select = as.data.frame(MCF7vsTR_DownR[MCF7vsTR_DownR$symbol %in% rownames(sheets_RTs_TRs_DLGs[[current_key]]),])
  current_dlg_df = as.data.frame(cellTR_list[paste(sample,"Gained",sep="_")])
  current_dlg_gene_col = paste(paste(sample,"Gained",sep="_"),"Genes",sep='.')
  current_dlg_df_select = current_dlg_df[current_dlg_df[,current_dlg_gene_col] %in% current_deg_df_select$symbol, ]
  current_merge = merge(current_deg_df_select,current_dlg_df_select,by.x='symbol',by.y=current_dlg_gene_col,all=TRUE)
  print(paste(i,dim(current_df)[1],sep=" "))
  if (dim(current_merge)[1]!=0){
    out_key = paste(i,sample,"Gained","DownR",sep="_")
    sheets_RTs_TRs_DLGs_DEGs[[out_key]] = current_merge
  }
}

for ( i in new_gain_cols){
  sample = "T47DTR"
  current_key = paste(i,sample,"Gained",sep="_")
  current_deg_df_select = as.data.frame(T47DvsTR_DownR[T47DvsTR_DownR$symbol %in% rownames(sheets_RTs_TRs_DLGs[[current_key]]),])
  current_dlg_df = as.data.frame(cellTR_list[paste(sample,"Gained",sep="_")])
  current_dlg_gene_col = paste(paste(sample,"Gained",sep="_"),"Genes",sep='.')
  current_dlg_df_select = current_dlg_df[current_dlg_df[,current_dlg_gene_col] %in% current_deg_df_select$symbol, ]
  current_merge = merge(current_deg_df_select,current_dlg_df_select,by.x='symbol',by.y=current_dlg_gene_col,all=TRUE)
  print(paste(i,dim(current_df)[1],sep=" "))
  if (dim(current_merge)[1]!=0){
    out_key = paste(i,sample,"Gained","DownR",sep="_")
    sheets_RTs_TRs_DLGs_DEGs[[out_key]] = current_merge
  }
}

# RTs-TRs-Loss-UpR
for ( i in new_loss_cols){
  sample = "MCF7TR"
  current_key = paste(i,sample,"Lost",sep="_")
  current_deg_df_select = as.data.frame(MCF7vsTR_UpR[MCF7vsTR_UpR$symbol %in% rownames(sheets_RTs_TRs_DLGs[[current_key]]),])
  current_dlg_df = as.data.frame(cellTR_list[paste(sample,"Lost",sep="_")])
  current_dlg_gene_col = paste(paste(sample,"Lost",sep="_"),"Genes",sep='.')
  current_dlg_df_select = current_dlg_df[current_dlg_df[,current_dlg_gene_col] %in% current_deg_df_select$symbol, ]
  current_merge = merge(current_deg_df_select,current_dlg_df_select,by.x='symbol',by.y=current_dlg_gene_col,all=TRUE)
  print(paste(i,dim(current_df)[1],sep=" "))
  if (dim(current_merge)[1]!=0){
    out_key = paste(i,sample,"Lost","UpR",sep="_")
    sheets_RTs_TRs_DLGs_DEGs[[out_key]] = current_merge
  }
}

for ( i in new_loss_cols){
  sample = "T47DTR"
  current_key = paste(i,sample,"Lost",sep="_")
  current_deg_df_select = as.data.frame(T47DvsTR_UpR[T47DvsTR_UpR$symbol %in% rownames(sheets_RTs_TRs_DLGs[[current_key]]),])
  current_dlg_df = as.data.frame(cellTR_list[paste(sample,"Lost",sep="_")])
  current_dlg_gene_col = paste(paste(sample,"Lost",sep="_"),"Genes",sep='.')
  current_dlg_df_select = current_dlg_df[current_dlg_df[,current_dlg_gene_col] %in% current_deg_df_select$symbol, ]
  current_merge = merge(current_deg_df_select,current_dlg_df_select,by.x='symbol',by.y=current_dlg_gene_col,all=TRUE)
  print(paste(i,dim(current_df)[1],sep=" "))
  if (dim(current_merge)[1]!=0){
    out_key = paste(i,sample,"Lost","UpR",sep="_")
    sheets_RTs_TRs_DLGs_DEGs[[out_key]] = current_merge
  }
}

write_xlsx(sheets_RTs_TRs_DLGs_DEGs,paste(datdir,"RTs_TRs_DLGs_DEGs.xlsx",sep='/'))

# circlize heatmap
plot_dlg_deg <- function(samplename,dlg_deg_type,sheets_RTs_TRs_DLGs_DEGs,groupname){
  split <- c()
  deg_mat_list = list()
  dlg_mat_list = list()
  for (current_type in dlg_deg_type){
    sample = samplename
    # deg info
    current_key = paste(">=2RTs",sample,current_type,sep="_")
    current_df = sheets_RTs_TRs_DLGs_DEGs[[current_key]]
    rownames(current_df) = current_df$symbol
    split <- c(split,rep(current_type,dim(current_df)[1]))
    deg_mat_list[[current_type]] = current_df[,7:12]
    
    # dlg_info
    gain_loss_type = str_split(current_type,"_")[[1]][1]
    current_dlg_df = current_df[,(dim(current_df)[2]-1):dim(current_df)[2]]
    colnames(current_dlg_df) = c("SubDiff","DivDiff")
    current_common = as.data.frame(RTsvsPTs_common_list[[paste(">=2RTs",gain_loss_type,sep="_")]])
    current_common_select =current_common[current_common$Genes %in% rownames(current_df),]
    rownames(current_common_select) = current_common_select$Genes
    current_common_select = current_common_select[rownames(current_dlg_df),]
    current_dlg_df = cbind(current_dlg_df,"Common_RTs"=as.character(current_common_select[,"Commons"]))
    dlg_mat_list[[current_type]] = current_dlg_df
  }
  # factor split vector
  split <- factor(split, levels = dlg_deg_type)
  # deg
  deg_mat = dplyr::bind_rows(deg_mat_list)
  deg_mat_z = t(scale(t(deg_mat)))
  colnames(deg_mat_z) = c(paste0(rep(paste(groupname[1],"rep",sep='_')),1:3),paste0(rep(paste(groupname[1],"rep",sep='_')),1:3))
  # dlg
  dlg_mat = dplyr::bind_rows(dlg_mat_list)
  dlg_mat[,"Scaled(DivDiff)"] = log2(dlg_mat$DivDiff)
  dlg_mat[,"Scaled(SubDiff)"] = sign(dlg_mat$SubDiff) * log2(abs(dlg_mat$SubDiff))
  dlg_mat[,"CexSubDiff"] = 2*abs(dlg_mat$`Scaled(SubDiff)`)/max(dlg_mat$`Scaled(SubDiff)`)
  
  circos.clear()
  num_sectors = length(dlg_deg_type)
  # set gaps for each sectors
  circos.par(gap.after = c(rep(2,num_sectors-1), 20))
  # set heatmap colors 
  col_fun1 = colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red"))
  # set dend colors
  dend_col = structure(c("red","orange","yellow3","cyan"), names = dlg_deg_type)
  # plot heatmap
  circos.heatmap(deg_mat_z, split = split, col = col_fun1, 
                 dend.side = "outside",
                 dend.track.height = 0.1,
                 dend.callback = function(dend, m, si) {
                   # when k = 1, it renders one same color for the whole dendrogram
                   color_branches(dend, k = 1, col = dend_col[si])
                 }
  )
  # add group rectangel
  circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
    if(CELL_META$sector.numeric.index == num_sectors) { # the last sector
      circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), 0,
                  CELL_META$cell.xlim[2] + convert_x(5, "mm"), dim(deg_mat_z)[2]/2,
                  col = "orange", border = NA)
      circos.text(CELL_META$cell.xlim[2] + convert_x(3, "mm"), dim(deg_mat_z)[2]/4,
                  groupname[2], cex = 0.5, facing = "clockwise")
      
      circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), dim(deg_mat_z)[2]/2,
                  CELL_META$cell.xlim[2] + convert_x(5, "mm"), dim(deg_mat_z)[2],
                  col = "pink", border = NA)
      circos.text(CELL_META$cell.xlim[2] + convert_x(3, "mm"), 3*dim(deg_mat_z)[2]/4,
                  groupname[1], cex = 0.5, facing = "clockwise")
    }
  }, bg.border = NA)
  # add dlg point plot
  max_axs_value = ceiling(max(abs(dlg_mat$`Scaled(DivDiff)`)))
  circos.track(ylim = c(-max_axs_value,max_axs_value), panel.fun = function(x, y) {
    y = dlg_mat[CELL_META$subset,"Scaled(DivDiff)"]
    y = y[CELL_META$row_order]
    circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "grey")
    circos.points(seq_along(y) - 0.5, y, col = ifelse(y > 0, "red", "blue") ,cex=dlg_mat[CELL_META$subset,"CexSubDiff"])
  }, cell.padding = c(0.02, 0, 0.02, 0))
  # add y-axis label for point plot
  circos.track(track.index = get.current.track.index(), 
               panel.fun = function(x, y) {
                 if(CELL_META$sector.numeric.index == num_sectors) { # the last sector
                   circos.text(x =rep(CELL_META$cell.xlim[2], 2) + convert_x(1, "mm"),
                               y=c(-max_axs_value,max_axs_value)+c(0,0),
                               c(-max_axs_value,max_axs_value),
                               cex = 1, adj = c(0, 0.5), facing = "inside")
                 }
               }, bg.border = NA)
  # add RTs information
  col_common_type = structure(c("purple", "blue", "green"), names = unique(dlg_mat$Common_RTs))
  common_RTs_list = as.data.frame(dlg_mat[,"Common_RTs"])
  circos.heatmap(common_RTs_list, col = col_common_type, track.height = 0.01,rownames.side="inside")
  # add gene names inside
  circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
    circos.text(1:CELL_META$cell.xlim[2], CELL_META$cell.ylim[1] - convert_y(2.5, "mm"), 
                rownames(dlg_mat[CELL_META$subset,][CELL_META$row_order,]),
                facing = "clockwise", cex = 0.6,
                adj = c(0.5, 0), niceFacing = TRUE,col = dend_col[CELL_META$sector.index])
  }, bg.border = NA)
  
  # add legend
  library(ComplexHeatmap)
  lgd_exprs = Legend(title = "Expression", col_fun = col_fun1)
  lgd_type = Legend(title = "Type", at = names(dend_col), 
                    legend_gp = gpar(fill = dend_col))
  lgd_tr = Legend(title = "RTs Common", at = names(col_common_type), 
                  legend_gp = gpar(fill = col_common_type))
  lgd_sub = Legend(title="Scaled(SubDiff)",labels = c("0.2","0.6","1.2","2"),
                   graphics = list(
                     function(x, y, w, h) {grid.points(x, y, gp = gpar(col = "black"), pch = 1,size=unit(0.1,'npc') )},
                     function(x, y, w, h) {grid.points(x, y, gp = gpar(col = "black"), pch = 1,size=unit(0.2,'npc') )},
                     function(x, y, w, h) {grid.points(x, y, gp = gpar(col = "black"), pch = 1,size=unit(0.3,'npc'))},
                     function(x, y, w,h) {grid.points(x, y, gp = gpar(col = "black"), pch = 1, size=unit(0.4,'npc'))}
                   ))
  
  
  h = dev.size()[2]
  lgd_list = packLegend(lgd_exprs, lgd_type, lgd_tr, lgd_sub,max_height = unit(0.5*h, "inch"))
  grid.draw(lgd_list)
  return(list("deg"=deg_mat,"dlg"=dlg_mat))
}
# MCF7TR
mcf7tr_deg_dlg_mat_list <- plot_dlg_deg("MCF7TR",dlg_deg_type,sheets_RTs_TRs_DLGs_DEGs,c("MCF7","MCF7TR"))
# T47DTR
t47dtr_deg_dlg_mat_list <- plot_dlg_deg("T47DTR",dlg_deg_type,sheets_RTs_TRs_DLGs_DEGs,c("T47D","T47DTR"))


# # RTs common TRs, alternative way to find RTs_TRs_DLGs
# RTTR_Gain_common = binary_Gain[rowSums(binary_Gain[rt_samples])>=1&rowSums(binary_Gain[cell_samples])>=1,]
# RTTR_Loss_common = binary_Loss[rowSums(binary_Loss[rt_samples])>=1&rowSums(binary_Loss[cell_samples])>=1,]
# RTTR_Gain_common["RTs_MCF7TR"] = rowSums(RTTR_Gain_common[,rt_samples]) * RTTR_Gain_common['MCF7TR']
# RTTR_Gain_common["RTs_T47DTR"] = rowSums(RTTR_Gain_common[,rt_samples]) * RTTR_Gain_common['T47DTR']
# RTTR_Gain_common["RTs_ZR75.1TR"] = rowSums(RTTR_Gain_common[,rt_samples]) * RTTR_Gain_common['ZR75.1TR']
# RTTR_Loss_common["RTs_MCF7TR"] = rowSums(RTTR_Loss_common[,rt_samples]) * RTTR_Loss_common['MCF7TR']
# RTTR_Loss_common["RTs_T47DTR"] = rowSums(RTTR_Loss_common[,rt_samples]) * RTTR_Loss_common['T47DTR']
# RTTR_Loss_common["RTs_ZR75.1TR"] = rowSums(RTTR_Loss_common[,rt_samples]) * RTTR_Loss_common['ZR75.1TR']
# 
# 
# for (sample in cell_samples){
#   for ( i in names(table(RTTR_Gain_common[paste("RTs",sample,sep="_")]))){
#     current_df = RTTR_Gain_common[RTTR_Gain_common[paste("RTs",sample,sep="_")]>=i,]
#     current_df = cbind("Genes"=rownames(current_df), current_df)
#     current_dim = dim(current_df)[1]
#     if (i != 0 & current_dim !=0){
#       print(paste(sample,"Gained",i,current_dim,sep=" "))
#       sheets_RTs_TRs_DLGs[[paste(i,paste("RTs",sample,"Gained",sep="_"),sep=">=")]] = current_df
#     }
#   }
#   for ( i in names(table(RTTR_Loss_common[paste("RTs",sample,sep="_")]))){
#     current_df = RTTR_Loss_common[RTTR_Loss_common[paste("RTs",sample,sep="_")]>=i,]
#     current_df = cbind("Genes"=rownames(current_df), current_df)
#     current_dim = dim(current_df)[1]
#     if (i != 0 & current_dim !=0){
#       print(paste(sample,"Lost",i,current_dim,sep=" "))
#       sheets_RTs_TRs_DLGs[[paste(i,paste("RTs",sample,"Lost",sep="_"),sep=">=")]] = current_df
#     }
#   }
# }

