
# read data
ids <-  paste0("SI_", c(37326:37328, 37335))

sc.srt <-  lapply(ids, function(x) {
  print(x)
  s = readRDS(paste0("./sc_pipeline/", x, "/data/ambient.soupX.rds"))[[2]]
  m = read.csv(paste0("./sc_pipeline/", x, "/data/cellannot_filt.csv"),
               stringsAsFactors = F, row.names = 1)
  s = s[,rownames(m)]
  s = CreateSeuratObject(counts = s)
})
names(sc.srt) <-  ids 

for (i in 1:length(sc.srt)) {
  sc.srt[[i]] <-  sc.proc(sc.srt[[i]])
}

for (i in 1:length(sc.srt)) {
  sc.srt[[i]][["percent.mt"]] <- PercentageFeatureSet(sc.srt[[i]], pattern = "^mt-")
}

##################
# further qc
sc.qc <-  do.call(rbind, lapply(sc.srt, function(x) x@meta.data))
sc.qc$library <-  substr(rownames(sc.qc), 1, 8)
sc.qc$nCount_RNA.log <-  log10(sc.qc$nCount_RNA)
qc.density <-  lapply(c("nCount_RNA.log","nFeature_RNA","percent.mt"), function(x) {
  ggplot(data = sc.qc)+geom_density(aes_string(x, fill="library"), alpha =0.4)+theme_classic()+ggtitle(x)
})

sc.qc.pipeline <-  lapply(ids, function(x) {
  m = read.csv(paste0("./sc_pipeline/", x, "/data/cellannot_filt.csv"),
               stringsAsFactors = F, row.names = 1)
  m$library = x
  return(m)
})
sc.qc$doublet <-  unlist(lapply(sc.qc.pipeline, function(x) x$doublet.call2))

for (i in 1:length(sc.srt)) {
  sc.srt[[i]]$doublet <-  sc.qc.pipeline[[i]]$doublet.call2
}
DimPlot(sc.srt[[1]], group.by = "doublet")
mt.genes <-  grep("^mt-",rownames(sc.srt[[1]]), value = T)

sc.srt.flt <-  lapply(1:length(sc.srt), function(i) {
  srt = sc.srt[[i]][,sc.srt[[i]]$nFeature_RNA>=200 & sc.srt[[i]]$percent.mt<10 & sc.srt[[i]]$doublet=="singlet" ]
  srt = srt[!(rownames(srt) %in% mt.genes), ]
})

for (i in 1:length(sc.srt.flt)) {
  sc.srt.flt[[i]] <-  sc.proc(sc.srt.flt[[i]])
}


# annotate with inhouse reference (PMID: 39763898)
mctp.ref <-  readRDS("mctp.ref.rds")
mctp.ref <-  sc.proc(mctp.ref)

predictions <-  lapply(sc.srt.flt, function(x) {
  anchors <- FindTransferAnchors(reference = mctp.ref, query = x,
                                 dims = 1:30, reference.reduction = "pca")
  pred <- TransferData(anchorset = anchors, refdata = mctp.ref$anno.adj3,
                       dims = 1:30)
})
for (i in 1:length(sc.srt.flt)) {
  sc.srt.flt[[i]]$pred.new <-  predictions[[i]]$predicted.id
}


sc.pool <-  merge(sc.srt.flt[[1]], sc.srt.flt[-1], add.cell.ids = ids)
sc.pool$library <-  substr(colnames(sc.pool), 1, 8)
smp.anno <-  setNames(c("ctl","case1","case2","case2"),ids)
sc.pool$sample <-  smp.anno[sc.pool$library]
sc.pool.psk <-  get.psk(obj = sc.pool, anno = "pred.new", sample = "sample")
sc.pool$major <- sc.pool$pred.new
sc.pool$major[grepl("Luminal", sc.pool$major)] <-  "Luminal"
sc.pool.psk.major <-  get.psk(obj = sc.pool, anno = "major", sample = "sample")


trt.cls <- setNames(c("blue","orange"), c("ctrl","case"))
trt.cls1 <-  adjustcolor(trt.cls, 0.2)
pdf("fig.s3ab.class1.umap.anno.allcells.pdf", width = 5.4, height = 4, useDingbats = F)
DimPlot(sc.pool, group.by = "anno", label = T, repel = T,  label.size = 3.6)+NoLegend()
DimPlot(sc.pool, group.by = "trt", cols = trt.cls1)
dev.off()

df <-  Embeddings(sc.pool, "umap")
df <- cbind(df, sc.pool@meta.data[,c("pred","pred.new","library","sample","trt","anno")])
write.csv(df, "fig.s3ab.class1.umap.anno.allcells.csv")


#####################################################################

# subset basal and luminal 
# exclude Luminal-7,8 and 3 which are ventral/lateral/SV luminal cells because the samples were taken from anterior and dorsal
sc.epi <-  sc.pool[, grepl("Basal|Lum", sc.pool$anno) & !(sc.pool$pred.new %in% c("Luminal-8","Luminal-7","Luminal-3"))]
sc.epi <-  sc.proc(sc.epi, n = 500)

# read imputed data
imputed <- readRDS("imputed.rds")
sc.epi[['imputed']] <-  CreateAssayObject(counts = imputed)
DefaultAssay(sc.epi) <-  "imputed"

# annotate cells based on puroR and hFOXA1
d <-  t(as.matrix(sc.epi@assays$imputed@data[c("puroR","hFOXA1"),]))
d <-  cbind(d, sc.epi@meta.data)
d$puroR <-  ifelse(d$puroR>0, "puroR+","puroR-")
d$hFOXA1 <-  ifelse(d$hFOXA1>0, 'hFOXA1+','hFOXA1-')
table(d$puroR, d$hFOXA1)
d$group <-  paste(d$puroR, d$hFOXA1, sep = ".")
d$cell <-  1
dtab <-  aggregate(cell~trt+anno+group, FUN=sum, data=d)

gp.cls <-  setNames(c("blue","red","grey","grey"), sort(unique(d$group)))
dtab$group <-  factor(dtab$group, levels = rev(names(gp.cls)))
ggplot(data = dtab, aes(x = trt, y = cell, fill = group ))+geom_bar(stat = "identity")+facet_wrap(~anno, scales = "free")+theme_bw()+
  scale_fill_manual(values = gp.cls)

sc.epi$group <-  d$group
DimPlot(sc.epi, group.by = "group")


sc.epi$anno.adj <-  NA
sc.epi$anno.adj[which(sc.epi$group %in% c("puroR-.hFOXA1-","puroR-.hFOXA1+"))] <-  "Lum"
sc.epi$anno.adj[grepl("Basal",sc.epi$anno) & grepl("puroR[+]", sc.epi$group) ] <-  "Basal"
DimPlot(sc.epi, group.by = "ann1", split.by = "trt")

DefaultAssay(sc.epi) <-  "RNA"
selgenes <-  c("puroR","hFOXA1","Foxa1","Nkx3-1","Ar")
pdf("fig.s3c.lum.vlnplots.hFOXA1.etc.pdf", width = 3.6, height = 4, useDingbats = F)
VlnPlot(sc.epi[, which(sc.epi$anno.adj=="Lum")],group.by = "anno.adj", features = selgenes, split.by = "trt", cols = trt.cls, stack = T, flip = T)
dev.off()

df <-  t(sc.epi@assays$RNA@data[selgenes, ])
df <-  cbind(df, sc.epi@meta.data[, c("anno.adj","trt")])
write.csv(df, "fig.s3c.lum.vlnplots.hFOXA1.etc.csv")


#################################################
# denovo clustering
sc.epi <-  FindClusters(sc.epi, resolution = 0.2)
sc.epi$cluster <-  as.integer(sc.epi$RNA_snn_res.0.2)

sc.epi$cluster.reorder <-  sc.epi$cluster
sc.epi$cluster.reorder[sc.epi$cluster==3] <-  5
sc.epi$cluster.reorder[sc.epi$cluster==4] <-  3
sc.epi$cluster.reorder[sc.epi$cluster==5] <-  4

jco.cols <-  pal_jco()(5)[c(1,2,4,5,3)]
names(jco.cols) <-  1:5
pdf("fig.1i.epi.louvain.clusters.pdf", width = 4, height = 3.2, useDingbats = F)
DimPlot(sc.epi, label = T, group.by = "cluster.reorder", label.size = 5, repel = T, pt.size = 0.2, cols = jco.cols)
dev.off()
pdf("fig.1i.epi.louvain.clusters.split.pdf", width = 7, height = 3.2, useDingbats = F)
DimPlot(sc.epi, label = T, group.by = "cluster.reorder", label.size = 5, repel = T, pt.size = 0.2, cols = jco.cols, split.by = "trt")
dev.off()

df <- Embeddings(sc.epi, "umap")
df <- cbind(df, sc.epi@meta.data[,c("cluster.reorder","trt")])
write.csv(df, "fig.1i.csv")

####################################
vbls <-  c("Ar","hFOXA1","Foxa1","Nkx3-1","Myc") 
sc.epi.meta <-  cbind(sc.epi@meta.data, t(as.matrix(sc.epi@assays$imputed@data[vbls,])))
cluster.mean <-  aggregate(.~cluster.reorder, data = sc.epi.meta[, c("cluster.reorder", vbls)], FUN=median)
vln.cols <-  brewer.pal(5, "RdBu")
custom.vln <-  function(x, ylab = "score") {
  d = cluster.mean[order(-cluster.mean[,x]),]
  col.val = setNames(vln.cols, d$cluster.reorder)
  dt = sc.epi1.meta
  dt$cluster = as.factor(dt$cluster.reorder)
  colnames(dt)[colnames(dt)==x]= "score"
  ggplot(data = dt, aes(cluster, score, fill = cluster))+
    geom_violin(scale = "width")+ylab(ylab)+
    scale_fill_manual(values = col.val)+ggtitle(x)+
    theme_bw()+theme()
}
ps <-  lapply(vbls,function(x) custom.vln(x, ylab = "log2UMI"))
pdf("fig.1j.epi.violin.Ar.etc.imputed.pdf", width = 3.5, height = 16, useDingbats = F)
gridExtra::grid.arrange(grobs = ps, ncol=1)
dev.off()

write.csv(sc.epi.meta[, c("cluster.reorder", vbls)], "fig.1j.csv")
