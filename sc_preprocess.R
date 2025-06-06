#!/usr/bin/env Rscript


print('sc_pipeline version 6.0')
#########################################################################
# 0) prep
# 0A. input

suppressMessages(library(optparse))
option_list <- list( 
  make_option(c("-i", "--lib"), type='character', default = NULL,
              help="library ID"),
  #  make_option(c("-s", "--species"), type='character', default = 'human', 
  #              help="human or mouse"),
  make_option(c("-d", "--inpath"), type='character', default = './input/', 
              help="cellranger count path"),
  make_option(c("-o", "--outpath"), type='character', default = './sc_pipeline/', 
              help="output path")
)
opt <- parse_args(OptionParser(option_list=option_list))


suppressMessages(library(ComplexHeatmap))
suppressMessages(library(dplyr))
suppressMessages(library(Seurat))
suppressMessages(library(Matrix))
suppressMessages(library(DropletUtils))
suppressMessages(library(SoupX))


run.name <- opt$lib
#run.species <-  opt$species 
count.path <- opt$inpath
setwd(opt$outpath)

# 0B. load custom functions
source("./functions.R")


# 0C. create the standard directories for each sample
if (!dir.exists(run.name)) {system(paste('mkdir', run.name))}
setwd(paste0(opt$outpath, run.name))
system('mkdir -p data plots')

run.species <- 'mouse'
#########################################################################
# 1) raw count
h5path <- paste0(count.path, '/', run.name, '/outs/filtered_feature_bc_matrix.h5')
count_data <-  Read10X_h5(h5path)
if (class(count_data)=='list') {
  count_data=count_data[[1]]
}
###########################################################
# 2) QC

# 2A. Create Seurat Object without filtering, and calculate fraction of Mito reads
srt <-  CreateSeuratObject(counts = count_data[rowSums(count_data)>0, ], project = "10x", min.cells = 0, min.features = 0)


if (run.species=="human"){
  mt.pattern = "^MT-"
} else if (run.species=="mouse"){
  mt.pattern = "^mt-"
} else {
  print("enter run.species")
}

srt$percent.mt = PercentageFeatureSet(srt, pattern = mt.pattern)


################################################################
# . Doublet detection
srt <-  sc.proc(srt)
sce <-  as.SingleCellExperiment(DietSeurat(srt,  graphs = "pca"))

set.seed(100)
dbl.scores1 <-  scDblFinder::computeDoubletDensity(sce, subset.row= VariableFeatures(srt), 
                                               d=ncol(reducedDim(sce)))
dbl.calls1 <-  scDblFinder::doubletThresholding(data.frame(score=dbl.scores1),
                                 method="griffiths", returnType="call")

set.seed(10010101)
dbl.calls2 <-  scDblFinder::scDblFinder(sce, clusters= Idents(srt))

dbl.df <-  data.frame(cell = colnames(srt), doublet.score1 = dbl.scores1, doublet.call1 = dbl.calls1,
                    doublet.score2 = dbl.calls2$scDblFinder.score, doublet.call2 = dbl.calls2$scDblFinder.class, stringsAsFactors = F)


#######################################################
# 2B. identify outliers in Libsize, Feature, Mito content 
libsize.drop <-  scater::isOutlier(srt$nCount_RNA, nmads = 3, type = "both", log = TRUE)
feature.drop <-  scater::isOutlier(srt$nFeature_RNA, nmads = 3, type = "both", log = TRUE)
mito.drop <-  scater::isOutlier(srt$percent.mt, nmads = 3, type = "higher") | srt$percent.mt>=80

cellfilt.df <- c()
cellfilt.df$mito <- mito.drop
cellfilt.df$libsize <- libsize.drop
cellfilt.df$feature <- feature.drop

cellfilt.df <- as.data.frame(cellfilt.df)
rownames(cellfilt.df) <- colnames(count_data)

# write out cellfilter 
cellfilt.df <- cbind.data.frame(cellfilt.df, srt@meta.data[,-1])
cellfilt.df <-  cbind.data.frame(cellfilt.df, dbl.df)
write.csv(cellfilt.df, "data/cellannot_full.csv")

###################################################################
# 3. Clustering 
srt <-  sc.proc(srt) # by default, 2000 genes selected


#################################################################
# 5. ambinent RNA
sc <- load10X(paste0(count.path, run.name, '/outs/')) # for multiome libraries, need to copy clustering folder under gex to a level up

f <- 1
n <- 0
pdf('plots/soupX.pdf')
sc.try <- try(autoEstCont(sc, tfidfMin=f)) 

while ('try-error' %in% class(sc.try) & n<=10) {
  f=f-0.1
  n=n+1
  sc.try=try(autoEstCont(sc, tfidfMin=f))
}

if ('try-error' %in% class(sc.try)) {
  out='error'
} else {
  out=adjustCounts(sc.try)
  sink( 'data/ambient.SoupX.txt')
  print(paste('estimated contamination by SoupX:', sc.try$fit$rhoEst))
  sink()
  saveRDS(list(sc=sc.try, out=out), 'data/ambient.soupX.rds')
  
}
dev.off()
rm(sc)

raw <- Read10X_h5(paste0(count.path, run.name, '/outs/raw_feature_bc_matrix.h5'))
if (class(raw)=='list') {
  raw=raw[[1]]
}
amb.est <- estimateAmbience(raw, round=FALSE, good.turing=FALSE, lower = min(srt$nCount_RNA)) # 
amb.res <- sort(amb.est, decreasing = T)
amb.res <- data.frame(gene=names(amb.res), count=amb.res, logNormUMI=log2(1e4*amb.res/sum(amb.res)+1))
write.csv(amb.res,  'data/ambient.DropletUtils.csv')





