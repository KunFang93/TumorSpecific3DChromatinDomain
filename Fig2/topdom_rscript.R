source("/data/kfang/software/TopDom/TopDom_v0.0.2.R")
#library(TopDom)
library(purrr)
# args 1 is full path of files, args[2] is full path of result, args[3] is the window size
args <- commandArgs(TRUE)
file_dir <- args[1]
result_dir <- args[2]
win_size<-args[3]
dir.create(file.path(dirname(result_dir)), showWarnings = FALSE)
dir.create(file.path(result_dir),showWarnings = FALSE)
files<-list.files(path=file_dir,pattern="*dense*")
for (file in files)
{
        print(file)
	print(strsplit(file,"_")[[1]][4])
#       print(paste(file_dir,file,sep="/"))
#	print(result_dir)
	res<- try(TopDom(matrix.file=paste(file_dir,file,sep="/"),window.size=win_size,outFile=paste(result_dir,file,sep="/")))
        if(inherits(res,"try-error")){
            next
        }
}
