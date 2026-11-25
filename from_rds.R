# Export native image, coordinates, measured library depth and target counts.
Sys.setlocale("LC_CTYPE","English_United States.utf8")
suppressPackageStartupMessages(library(Matrix))
args <- commandArgs(trailingOnly=TRUE)
decode <- function(s) {i <- seq(1,nchar(s),2); value <- rawToChar(as.raw(strtoi(substring(s,i,i+1L),16L))); Encoding(value) <- "UTF-8"; value}
input <- decode(args[1]); dest <- decode(args[2])
dir.create(dest,showWarnings=FALSE,recursive=TRUE)
genes <- c("Acan","Col2a1","Col1a1","Col9a1","Comp","Dcn","Bgn","Fmod","Sox9","Krt19","Mmp3","Mmp13","Adamts5")
files <- list.files(input,"^new_(Control|Day1|Day3|Day28)(_st|B)[.]rds$",full.names=TRUE)
stopifnot(length(files)==8)
for (f in files) {
 obj <- readRDS(f); assay <- obj@assays[["Spatial"]]; mat <- assay@counts
 img <- obj@images[[1]];bc <- colnames(mat);co <- img@coordinates[bc,,drop=FALSE]
 sec <- unique(as.character(obj@meta.data[bc,"new_sampleid"])); symbols <- toupper(rownames(assay@meta.features))
 stopifnot(length(sec)==1,identical(rownames(co),bc),!anyDuplicated(symbols),!anyDuplicated(bc))
 ids <- match(toupper(genes),symbols);stopifnot(!anyNA(ids))
 tab <- data.frame(section=sec,barcode=bc,full_x=co$imagecol,full_y=co$imagerow,lowres_scale=img@scale.factors$lowres,total=as.numeric(Matrix::colSums(mat)),t(as.matrix(mat[ids,,drop=FALSE])),check.names=FALSE)
 names(tab)[7:ncol(tab)] <- genes
 write.csv(tab,file.path(dest,paste0(sec,".csv")),row.names=FALSE)
 png::writePNG(img@image,file.path(dest,paste0(sec,".png")))
 cat(sec,nrow(tab),"native spots\n");flush.console()
 rm(obj,assay,mat,img,tab);gc()
}
