
library(Signac)
library(EnsDb.Mmusculus.v79)
library(data.table)
library(Matrix)
library(GenomicRanges)
library(JASPAR2020)
library(TFBSTools)
library(patchwork)

annotations <-  GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevels(annotations) = paste0('chr', seqlevels(annotations))
genome(annotations) = "mm10"

ids = readRDS("class2.young.old.ids.rds")

###############################################

atac.proc = function(x) {
  p = paste0("./cellranger_arc/count/", x) #"/outs/filtered_feature_bc_matrix.h5"
  
  d = Read10X_h5(paste0(p, "/outs/filtered_feature_bc_matrix.h5"))[[2]]
  m = read.csv(paste0(p, "/outs/per_barcode_metrics.csv"), 
               header = T, row.names = 1)
  f = paste0(p, "/outs/atac_fragments.tsv.gz")
  srt = CreateChromatinAssay(counts = d,  sep = c(":", "-"), genome = "mm10",
                             fragments = f, min.cells = 1)
  srt = CreateSeuratObject(
    counts = srt,
    assay = 'peaks',
    project = 'ATAC',
    meta.data = m[colnames(srt),]
  )
  Annotation(srt) = annotations
  srt = NucleosomeSignal(srt)
  srt = TSSEnrichment(srt, fast = TRUE)
  srt$pct_frag_in_peaks <- srt$atac_peak_region_fragments/srt$atac_fragments 
  srt$log_frag_in_peaks <- log10(srt$atac_peak_region_fragments)
  srt$library = x
  write.csv(srt@meta.data , paste0("qc/",x,".csv"))
  return(srt)
}
atac.list <-  lapply(ids, function(x) atac.proc(x))  

# read annotation of cells based on gene expression
gex.meta <-  readRDS("pool.epi.meta.rds")


peaklist <-  lapply(ids, function(x) {
  p = paste0("./cellranger_arc/count/", x) 
  
  # read in peak sets
  peaks <- read.table(
    file = paste0(p,"/outs/atac_peaks.bed"),
    col.names = c("chr", "start", "end")
  )
})

# convert to genomic ranges
peaklist <-  lapply(peaklist, makeGRangesFromDataFrame )

# Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = c(peaklist[[1]], peaklist[[2]], peaklist[[3]],
                               peaklist[[4]], peaklist[[5]], peaklist[[6]],
                               peaklist[[7]], peaklist[[8]], peaklist[[9]],
                               peaklist[[10]], peaklist[[11]], peaklist[[12]]
                               ))

# check peaks width
peakwidths <- width(combined.peaks)
summary(peakwidths)
lapply(peaklist, function(x) summary(width(x)))

file.path <-  data.frame(lib = ids, stringsAsFactors = F)
file.path$path = sapply(ids, function(x) {
  p = paste0("./cellranger_arc/count/", x) #"/outs/filtered_feature_bc_matrix.h5"
  return(p)
})
# create fragment objects
qflt <- read.csv("qc.csv") 
qflt <-  q[q$atac_peak_region_fragments>=1000,]
fraglist <-  lapply(1:length(ids), function(i) {
  frags <- CreateFragmentObject(
    path = paste0(file.path$path[i],"/outs/atac_fragments.tsv.gz"),
    cells = qflt$V1[qflt$library==ids[i]]
  )
})
names(fraglist) <- ids
ctlist <-  lapply(ids, function(x) {
  counts <- FeatureMatrix(
    fragments = fraglist[[x]],
    features = combined.peaks,
    cells = qflt$V1[qflt$library==x]
  )
})

#######################

names(ctlist) <-  ids
assaylist <-  lapply(ids, function(x) {
  assayobj <- CreateChromatinAssay(ctlist[[x]], fragments = fraglist[[x]])
  assayobj <- CreateSeuratObject(assayobj, assay = "ATAC")
})

# merge all datasets, adding a cell ID to make sure cell names are unique
combined <- merge(
  x = assaylist[[1]],
  y = assaylist[-1],
  add.cell.ids = ids
)
combined[["ATAC"]]
rownames(qflt) <-  qflt$id
qflt <-  qflt[colnames(combined),]
combined$library <- substr(colnames(combined), 1, 17)
gex.meta.all <-  readRDS("pool.meta.rds")
combined$gex.anno.epi <-  gex.meta[match(colnames(combined), rownames(gex.meta)),]$pred.crowley.new
combined$gex.anno <-  gex.meta.all[match(colnames(combined), rownames(gex.meta.all)),]$pred.mctp
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 50)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, dims = 2:30, reduction = 'lsi')

DimPlot(combined, group.by = 'library', pt.size = 0.1)
DimPlot(combined,group.by = "gex.anno", label = T )
tmp <-  gex.meta.all[,c("library","trt","age","lobe","group")]
tmp <-  tmp[!duplicated(tmp),]
combined$trt <-  tmp$trt[match(combined$library, tmp$library)]
combined$age <-  tmp$age[match(combined$library, tmp$library)]
DimPlot(combined, group.by = "trt")+DimPlot(combined, group.by = "age")
combined$anno <-  combined$gex.anno
combined$anno <-  sapply(combined$anno, function(x) strsplit(x, "-")[[1]][1])
combined$anno[combined$gex.anno %in% c("Luminal-3","Luminal-5")] <- combined$gex.anno[combined$gex.anno %in% c("Luminal-3","Luminal-5")]


combined <- FindNeighbors(
  object = combined,
  reduction = 'lsi',
  dims = 2:30
)
combined <- FindClusters(
  object = combined,
  algorithm = 3,
  resolution = 1,
  verbose = FALSE
)
DimPlot(combined, label = T)
aggregate(nCount_ATAC~seurat_clusters, FUN=median, data= combined@meta.data)
saveRDS(combined, "combined.rds")

g1 <-  DimPlot(combined, label = T)+NoLegend()
g2 <-  DimPlot(combined, group.by = "anno", label = T)+NoLegend()
g1+g2
m <-  table(combined$anno, combined$seurat_clusters)
a <-  rownames(m)[apply(m, 2, which.max)]
names(a) <-  colnames(m)
combined$anno <-  a[as.character(combined$seurat_clusters)]
combined$anno[combined$anno=="Luminal-3"] <-  "SV"
pdf("umap.anno.pdf", width = 6, height = 5, useDingbats = F)
DimPlot(combined, group.by = "anno", label = T)+NoLegend()
DimPlot(combined, group.by = "gex.anno", label = T)+NoLegend()
dev.off()

trt.cls <-  setNames(c("cyan4","coral"), c("ctrl","case"))
lib.cls <-  setNames(brewer.pal(12, "Paired"), unique(combined$library))
pdf("umap.by.library.age.trt.pdf", width = 6, height = 4.5, useDingbats = F)
DimPlot(combined, group.by = "library", cols = lib.cls)
DimPlot(combined, group.by = "age")+scale_color_nejm()
DimPlot(combined, group.by = "trt", cols = trt.cls)
dev.off()


# keep lum cells
combined.lum <-  combined[, grepl("Lum", combined$anno)]
DefaultAssay(combined.lum) <-  "ATAC"
combined.lum <- RunTFIDF(combined.lum)
combined.lum <- FindTopFeatures(combined.lum, min.cutoff = 20)
combined.lum <- RunSVD(combined.lum)
combined.lum <- RunUMAP(combined.lum, dims = 2:50, reduction = 'lsi')
combined.lum <- FindNeighbors(
  object = combined.lum,
  reduction = 'lsi',
  dims = 2:50
)
combined.lum <- FindClusters(
  object = combined.lum,
  algorithm = 3,
  resolution = 0.5,
  verbose = FALSE
)

# add motif information
DefaultAssay(combined.lum) = "ATAC"
rn = rownames(combined.lum)
combined.lum = combined.lum[grepl("^chr",rn), ]

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)
pfm.names <-  sapply(pfm, function(x) x@name)
names(pfm) <-  pfm.names

combined.lum <- AddMotifs(
  object = combined.lum,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)
combined.lum$gex.cluster <-  gex.meta$cluster[match(colnames(combined.lum), rownames(gex.meta))]
DimPlot(combined.lum, group.by = "gex.cluster", label = T)
Idents(combined.lum) = as.factor(combined.lum$trt)
jco.cls <-  pal_jco()(7)[c(1,3,2,4:7)]
names(jco.cls) <-  1:7
g1 <-  DimPlot(combined.lum, label = T)+NoLegend()
g2 <-  DimPlot(combined.lum[,!is.na(combined.lum$gex.cluster)], group.by = "gex.cluster", label = T, cols = jco.cls)+NoLegend()
g1+g2

# manual rename clusters based on cluster assignment from GEX 
combined.lum$cluster = 2
combined.lum$cluster[combined.lum$seurat_clusters %in% c(0,6,12)] = 3
combined.lum$cluster[combined.lum$seurat_clusters %in% c(1)] = 4
combined.lum$cluster[combined.lum$seurat_clusters %in% c(10)] = 7
combined.lum$cluster[combined.lum$seurat_clusters %in% c(4)] = 6
combined.lum$cluster[combined.lum$seurat_clusters %in% c(8)] = 5
g3 <-  DimPlot(combined.lum, group.by = "cluster", cols = jco.cls, label = T)
g2+g3
table(combined.lum$cluster, combined.lum$trt)
combined.lum$cluster.group <-  "ctrl.clusters"
combined.lum$cluster.group[combined.lum$cluster %in% c(3,4)] <-  "case.clusters"
combined.lum$cluster.group[combined.lum$cluster==5] <-  NA
Idents(combined.lum) <-  as.factor(combined.lum$cluster.group)

pdf("lum.umap.anno.cluster.pdf", width = 13, height = 5, useDingbats = F)
DimPlot(combined.lum[,!is.na(combined.lum$gex.anno.epi)], group.by = "gex.anno.epi", label = T, split.by = "trt")+NoLegend()
DimPlot(combined.lum, group.by = "cluster", label = T, cols = jco.cls,  split.by = "trt")+NoLegend()
dev.off()
pdf("lum.umap.denocluster.pdf", width = 13, height = 5, useDingbats = F)
DimPlot(combined.lum, group.by = "seurat_clusters", label = T,  split.by = "trt")+NoLegend()
dev.off()

pdf("lum.umap.by.library.age.trt.pdf", width = 6, height = 4.5, useDingbats = F)
DimPlot(combined.lum, group.by = "library", cols = lib.cls)
DimPlot(combined.lum, group.by = "age")+scale_color_nejm()
DimPlot(combined.lum, group.by = "trt", cols = trt.cls)
dev.off()

################################################
# identify differential peaks: cluster 3/4 vs 2/6/7
combined.lum$cluster.flt <-  combined.lum$cluster.group
combined.lum$cluster.flt[combined.lum$cluster %in% c(3,4) & combined.lum$trt=="ctrl"] <-  NA
combined.lum$cluster.flt[combined.lum$cluster %in% c(2,6,7) & combined.lum$trt=="case"] <-  NA
table(combined.lum$cluster, combined.lum$cluster.flt)

DefaultAssay(combined.lum) <- "ATAC"
da_peaks.Lum.flt <- FindMarkers(
  object = combined.lum,
  ident.1 = 'case.clusters',
  ident.2 = 'ctrl.clusters',group.by = "cluster.flt",
  only.pos = F, logfc.threshold = 0.1,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_ATAC'
)
da_peaks.Lum.flt0 <- FindMarkers( # to get results from all peaks
  object = combined.lum,
  ident.1 = 'case.clusters',
  ident.2 = 'ctrl.clusters',group.by = "cluster.flt",
  only.pos = F, logfc.threshold = 0,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_ATAC'
)
da_peaks.Lum.flt0$sig <-  NA
da_peaks.Lum.flt0$sig[da_peaks.Lum.flt0$avg_log2FC>0.1 & da_peaks.Lum.flt0$p_val_adj<0.01] <-  "UP"
da_peaks.Lum.flt0$sig[da_peaks.Lum.flt0$avg_log2FC< -0.1 & da_peaks.Lum.flt0$p_val_adj<0.01] <- "DN"
saveRDS(da_peaks.Lum.flt0, "da_peaks.lum.flt0.rds")

pdf("case.vs.ctrl.peaks.volcano.pdf", width = 6, height = 5, useDingbats = F)
ggplot( data = da_peaks.Lum.flt0[is.na(da_peaks.Lum.flt0$sig),],aes(avg_log2FC, -log10(p_val_adj)))+
  geom_point(size = 0.2, color ="grey")+theme_classic(base_size = 12)+
  geom_point(data = da_peaks.Lum.flt0[!is.na(da_peaks.Lum.flt0$sig),], aes(color = sig), size =0.2)+
  scale_color_manual(values = c("UP"="red","DN"="Blue"))+
  geom_vline(xintercept = c(-0.1, 0.1), lty=3)+
  geom_hline(yintercept = 2, lty=3)+
  annotate("text", x = -0.4, y = c(220, 210,200,190), label = c("number of peaks:"," UP:8884"," DN:2161"," insignif:34509"), hjust=0)+
  annotate("text", x = -0.25, y =7, label ="p_val_adj = 0.01")
dev.off()



up.peaks.flt <-   da_peaks.Lum.flt[da_peaks.Lum.flt$avg_log2FC>0  & da_peaks.Lum.flt$p_val_adj<0.01,]
enriched.motifs.case.flt.top2k <- FindMotifs(
  object = combined.lum,
  features = rownames(up.peaks.flt)[1:2000]
)
enriched.motifs.case.flt <- FindMotifs(
  object = combined.lum,
  features = rownames(up.peaks.flt) # all ~8k peaks
)

dn.peaks.flt <-   da_peaks.Lum.flt[da_peaks.Lum.flt$avg_log2FC<0  & da_peaks.Lum.flt$p_val_adj<0.01,]
enriched.motifs.ctrl.flt.top2k <- FindMotifs(
  object = combined.lum,
  features = rownames(dn.peaks.flt)[1:2000]
)
enriched.motifs.ctrl.flt <- FindMotifs(
  object = combined.lum,
  features = rownames(dn.peaks.flt)
)

gex.meta$cluster_trt <-  paste(gex.meta$cluster, gex.meta$trt, sep = "_")
combined.lum$cluster_trt <-  paste(combined.lum$cluster, combined.lum$trt, sep = "_")


tfs <-  c("Tcf12","Tcf21","Tcf3","Tcf4","Tcfl2","Tcf7l1","Tcf7l2","Tcf7","Tcfcp211","Lef1","Sox9","Klf5","Gata2","Gata3","Hoxb13","Nkx3-1","Ar","Foxa1")
enriched.motifs.case.flt.top2k$gene_name <-  sapply(enriched.motifs.case.flt.top2k$motif, function(x) strsplit(x, "[(]")[[1]][1])
enriched.motifs.case.flt.top2k$gene_name <-  toupper(enriched.motifs.case.flt.top2k$gene_name)
tmp <-  enriched.motifs.case.flt.top2k[1:100,]
tmp <-  enriched.motifs.case.flt.top2k[enriched.motifs.case.flt.top2k$gene_name %in% union(tmp$gene_name, toupper(tfs)),]
tmp <-  tmp[!duplicated(tmp$gene_name),]


gex.atac.dotplot(g = tmp, filename = "case-peaks-top2k.motifs.more-tfs-added.pdf", split = T, col_flt = TRUE, zscore = FALSE, color.limits = TRUE)
gex.atac.dotplot(g = tmp, filename = "case-peaks-top2k.motifs.more-tfs-added.scale.pdf", split = T, col_flt = TRUE, zscore = TRUE)

DefaultAssay(combined.lum) <-  "chromvar"
pdf("featureplot.chromvar.selected.genes.pdf", width = 7, height = 3, useDingbats = F)
FeaturePlot(combined.lum, features = "Ar", split.by = "trt", max.cutoff = "q90",min.cutoff = 'q10')
FeaturePlot(combined.lum, features = "FOS", split.by = "trt", max.cutoff = "q90",min.cutoff = 'q10')
FeaturePlot(combined.lum, features = "JUN", split.by = "trt", max.cutoff = "q90",min.cutoff = 'q10')
dev.off()

saveRDS(combined.lum, "combined.lum.rds")


motif.volcano(x = enriched.motifs.case.flt.top2k, add.genes = c("Ar","Lef1","Foxa1"), w=7, h=2.5,
               filename = "volcano.motifs.enriched.in.case.top2k.peaks.v1.top.pdf", ylimits = c(30,100))
motif.volcano(x = enriched.motifs.case.flt.top2k, add.genes = c("Ar","Lef1","Foxa1","Foxa2","Foxa3"), 
               filename = "volcano.motifs.enriched.in.case.top2k.peaks.v2.bottom.pdf", ylimits = c(0,50))
