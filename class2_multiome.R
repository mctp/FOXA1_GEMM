# copy this file to class2_multiome.R
# R4.1.3 Seurat 4.1.0
setwd("/mctp/share/users/yupingz/projects/abhijit_multiome/class2_young_old/")

# library ids
ids <- readRDS("class2.young.old.ids.rds")

# read floating-RNA corrected (by SoupX) counts
rna.srt <- lapply(ids, function(x) {
  print(x)
  s = readRDS(paste0("/mctp/share/users/yupingz/projects/sc_pipeline/", x, "/data/ambient.soupX.rds"))[[2]]
  m = read.csv(paste0("/mctp/share/users/yupingz/projects/sc_pipeline/", x, "/data/cellannot_filt.csv"),
               stringsAsFactors = F, row.names = 1)
  s = s[,rownames(m)]
  s = CreateSeuratObject(counts = s)
})
names(rna.srt) <-  ids 

for (i in 1:length(rna.srt)) {
  rna.srt[[i]] <-  sc.proc(rna.srt[[i]])
}

for (i in 1:length(rna.srt)) {
  rna.srt[[i]][["percent.mt"]] <- PercentageFeatureSet(rna.srt[[i]], pattern = "^mt-")
}


##################
# mitochondrial genes
mt.genes <-  grep("^mt-",rownames(rna.srt[[1]]), value = T)

# filter genes and cells
rna.srt.flt <-  lapply(1:length(rna.srt), function(i) {
  srt = rna.srt[[i]][,rna.srt[[i]]$nFeature_RNA>=400 & rna.srt[[i]]$percent.mt<40 ]
  srt = srt[!(rownames(srt) %in% mt.genes), ]
})

# sample information
smp <-  data.frame(library = ids, trt = c(na.omit(libs$`Case/Control`[11:30]),"Control","Control"), stringsAsFactors = F)
smp$age <-  c(rep(c("14wk","66wk"),c(2,8)),"15wk","21wk")
smp$filter1 <-  sapply(rna.srt, ncol)
smp$filter2 <-  sapply(rna.srt.flt, ncol)
smp$lobe <-  c("whole","whole",rep(c("DVL","AL"), 4),"DVL","DVL")

# merge libraries
pool <-  merge(rna.srt.flt[[1]], rna.srt.flt[-1], add.cell.ids = ids)
pool$library <-  substr(colnames(pool), 1, 17)
pool$trt <-  smp$trt[match(pool$library, smp$library)]
pool$age <-  smp$age[match(pool$library, smp$library)]
pool <-  sc.proc(pool)

# annotate using reference dataset from Crowley et al
crowley <-  readRDS("/mctp/share/users/yupingz/projects/sc_prostate/mouse/elife_paper/elif.pool.rds")
i <-  which(grepl("Dying|stress|trans|prolif",crowley$anno))
crowley <-  crowley[, -i]
crowley <-  sc.proc(crowley, cluster = FALSE)

anchors <- FindTransferAnchors(reference = crowley, query = pool,
                               dims = 1:30, reference.reduction = "pca")
pred <- TransferData(anchorset = anchors, refdata = crowley$anno,
                     dims = 1:30)
pool$pred.crowley <-  pred$predicted.id
pool$trt <-  ifelse(pool$trt=="Case","case","ctrl")
pool$trt <-  factor(pool$trt, levels = c("ctrl","case"))

# annotate using inhouse reference data 
mctp.ref <-  readRDS("/mctp/share/users/hanbyul/projects/sc_prostate/for_paper/figure1/whole.rds")
mctp.ref <-  sc.proc(mctp.ref)
anchors <-  FindTransferAnchors(reference = mctp.ref, query = pool,
                               dims = 1:30, reference.reduction = "pca")
pred1 <- TransferData(anchorset = anchors, refdata = mctp.ref$anno.adj3,
                     dims = 1:30)
pool$pred.mctp <-  pred1$predicted.id

########################################
# rename and finalize annotation
pool$anno <-  pool$pred.mctp
pool$anno[grep("Fibro", pool$anno)] <-  "Fibroblast"
pool$anno[grep("Schwann", pool$anno)] <-  "Schwann"
pool$anno[grep("Macro", pool$anno)] <-  "Macrophage"

pool$anno[grepl("Lum|Basal", pool$pred.crowley)] = pool$pred.crowley[grepl("Lum|Basal",pool$pred.crowley)]
pool$anno[pool$seurat_clusters ==17] <-  "SV"
pool$anno[pool$pred.mctp=="NE"] <-  "NE"
# remove cells with conflicting prediction
i <-  which(grepl("Luminal", pool$pred.mctp) & !grepl("Lum|Basal", pool$pred.crowley))
pool <- pool[,-i]

pdf("figS8a.class2.umap.all.cells.pdf", width = 8, height = 3.5, useDingbats = F)
DimPlot(pool, group.by = "anno", label = T, repel = T, split.by = "trt",  label.size = 3)+NoLegend()
dev.off()

df <-  Embeddings(pool, "umap")
df <-  cbind(df, pool@meta.data)
write.csv(df, "figS8a.class2.umap.all.cells.csv")

##############################
# subset epithelial cells
pool.epi <-  pool[, grepl("Lum|Basal", pool$anno) ]
pool.epi <-  sc.proc(pool.epi, n = 1000, pc = 20)

# remove cluster 1 and 12 which had lower umi and cluster 1 has some Vim expression
pool.epi <-  pool.epi[, !(pool.epi$seurat_clusters %in% c(1,12))]
pool.epi <-  RunUMAP(pool.epi, dims = 1:5, reduction.name = "umap3")

##########################
# denovo clusters
pool.epi <-  FindClusters(pool.epi, resolution = 0.2)
pool.epi <-  FindClusters(pool.epi, resolution = 0.1)
pool.epi <-  FindClusters(pool.epi, resolution = 0.05)

pool.epi$cluster <-  pool.epi$RNA_snn_res.0.2
# merge cluster 7
pool.epi$cluster[pool.epi$cluster==7] <-  1
# recorder clusters
m <-  setNames(c(1:4, 6,5,7), 0:6)
jco.cls <-  pal_jco()(7)[c(1,3,2,4:7)]
pool.epi$cluster <-  m[as.character(pool.epi$cluster)]

pdf("fig3b.class2.epi.umap.by.cluster.pdf", width = 4, height = 3)
DimPlot(pool.epi, group.by = "cluster", label = T, reduction = "umap3", label.size = 3, cols = jco.cls) # label.box = TRUE,
dev.off()



# read signatures
ctrl.mks.list <-  readRDS("ctrl.mks.list.rds") # basal, lum and lumP markers from wildtype mouse prostate
cast.genes2 <-  readRDS("cast.signatures.refine.rds")
lump.markers <-  readRDS("lump.markers.rds")
l1 <-  openxlsx::read.xlsx("kirk2024.mmc2.xlsx", sheet = 4)
l1.genes <- l1$gene[1:100]
ar.gset = readRDS("ar.signatures.rds")
all.signatures <- c(ctrl.mks.list, cast.genes2, lump.markers)
all.signatures$Kirk.L1 <- l1.genes


# signature score
pool.epi$score.crowley.LumP = get.score(lump.markers$crowley, pool.epi@assays$RNA@data)
pool.epi$score.guo.LumC <-  get.score(lump.markers$guo, pool.epi@assays$RNA@data)
pool.epi$score.karthaus.Lum2 <-  get.score(lump.markers$karthaus, pool.epi@assays$RNA@data)
pool.epi$LumP.overlap <-  get.score(genelist = lump.markers$overlap, mat = pool.epi@assays$RNA@data )
pool.epi$Kirk.L1.top100 <- get.score(l1.genes, pool.epi@assays$RNA@data )
pool.epi$AR.UP <- get.score(ar.gset$AR.UP, pool.epi@assays$RNA@data )

pdf("fig.s8f.class2-epi.featureplot.lumP.score.pdf", width = 7, height = 12, useDingbats = F)
grd.cls = viridis::rocket(9, direction = -1)
FeaturePlot(pool.epi, features = c("score.crowley.LumP","score.guo.LumC","score.karthaus.Lum2","LumP.overlap"), split.by = "trt",
            cols = grd.cls, reduction = "umap3", max.cutoff = "q90")
dev.off()


pool.epi$Castration20.refine <-  get.score(cast.genes2$castrated20, pool.epi@assays$RNA@data)
pool.epi$score.Lum <-  get.score(ctrl.mks.list$Lum, pool.epi@assays$RNA@data)
pool.epi$score.Basal <-  get.score(ctrl.mks.list$Basal, pool.epi@assays$RNA@data)
pool.epi$LumP.score <-  get.score(ctrl.mks.list$LumP, pool.epi@assays$RNA@data)



pdf("fig.s9b.class2-epi.featureplot.castration.score.pdf", width = 7, height = 3, useDingbats = F)
FeaturePlot(pool.epi, features = c("Castration20.refine"), 
            split.by = "trt", cols = grd.cls1, max.cutoff = "q90", reduction = "umap3")
dev.off()



##################################################################################
###################################################################################


#########################################

pool.epi$cluster = factor(pool.epi$cluster, levels = c(1,2,6,7,3,4,5))


###########################################
# color gradient, with cell filtering
i <-  which(pool.epi$cluster %in% c(3,4) & pool.epi$trt=="ctrl")
j <-  which(pool.epi$cluster %in% c(2,6,7,5) & pool.epi$trt == "case")
table(pool.epi$cluster[-c(i,j)], pool.epi$trt[-c(i,j)])
tmp.scores <-  lapply(all.signatures, function(x) {
  get.score(x, pool.epi@assays$RNA@data[, -c(i,j)])
})
tmp.scores <-  do.call(cbind, tmp.scores)
tmp.scores <-  as.data.frame(tmp.scores)
tmp.scores$cluster <-  pool.epi$cluster[-c(i,j)]

vbls <-  c("Castration20.refine","LumP.overlap","score.crowley.LumP","score.guo.LumC","score.karthaus.Lum2","Kirk.L1.top100","score.Basal","score.Lum","score.LumP")
cluster.mean <-  aggregate(.~cluster, data = tmp.scores, FUN=median)
ps <-  lapply(vbls, function(x) custom.vln(x = x, dt = tmp.scores))
pdf("fig3cdh.class2-epi.vln.signatures.pdf", width = 3.5, height = 2, useDingbats = F)
ps
dev.off()


######################################
# violin plots of individual genes, imputation, cell filtering
# read imputed matrix from magic
pool.epi[["imputed"]] <- CreateAssayObject(counts = impute)
DefaultAssay(pool.epi) <-  "imputed"

g <-  c("Ppp1r1b","Nkx3-1","Trp63","Hoxb13","Tacstd2","Psca","Krt4")
tmp <-  cbind(pool.epi@meta.data[-c(i,j),], t(impute[g, -c(i,j)]))
cluster.mean <-  aggregate(.~cluster, data = tmp[, c("cluster", g)], FUN=median)
ps <-  lapply(g, function(x) custom.vln(x, dt = tmp, ylim = TRUE))
pdf("fig3f.s8c.class2-epi.vln.markers.imputed.pdf", width = 3.5, height = 2, useDingbats = F)
ps
dev.off()

write.csv(tmp[, c("cluster",g)], "fig3f.s8c.class2-epi.vln.markers.imputed.csv")


#######################################
pdf("fig.s8g.class2-epi.vln.LumP.score.stat.pdf", width = 15, height = 3, useDingbats = F)
ggplot(data = pool.epi@meta.data, aes(trt, score.crowley.LumP))+geom_violin(aes(fill = trt), scale = "width" )+
  facet_wrap(~pred.crowley.new, scales = "free_y", ncol =7)+theme_classic()+stat_compare_means()+
  scale_fill_manual(values = trt.cls)+theme(legend.position = "none")
ggplot(data = pool.epi@meta.data, aes(trt, score.guo.LumC))+geom_violin(aes(fill = trt), scale = "width" )+
  facet_wrap(~pred.guo, scales = "free_y", ncol=5)+theme_classic()+stat_compare_means()+
  scale_fill_manual(values = trt.cls)+theme(legend.position = "none")
ggplot(data = pool.epi@meta.data, aes(trt, score.karthaus.Lum2))+geom_violin(aes(fill = trt), scale = "width" )+
  facet_wrap(~pred.karthaus, scales = "free_y", ncol=5)+theme_classic()+stat_compare_means()+
  scale_fill_manual(values = trt.cls)+theme(legend.position = "none")
dev.off()


######################
# differential analysis between clusters 3+4 from case vs clusters 2+6+7 from control, at pseudo-bulk level
pool.epi$group <-  paste(pool.epi$cluster, pool.epi$trt, sep = "_")
pool.epi.psk <-  get.psk(obj = pool.epi, anno = "group", sample = "library")
pool.epi.psk$dgelist$samples$trt <-  substr(pool.epi.psk$dgelist$samples$anno, 3, 6)
pool.epi.psk$dgelist$samples$cluster <- substr(pool.epi.psk$dgelist$samples$anno, 1, 1)
pool.epi.psk$dgelist$samples$include <- pool.epi.psk$dgelist$samples$anno %in% c("3_case","4_case","2_ctrl","6_ctrl","7_ctrl")


i <- which(pool.epi.psk$dgelist$samples$include)
pool.epi.psk$dgelist$samples$trt <-  factor(pool.epi.psk$dgelist$samples$trt, levels = c("ctrl","case"))
des <-  model.matrix(~trt, data = pool.epi.psk$dgelist$samples[i,])
flt <-  rowMeans(pool.epi.psk$logtpm[,i])>1

v <-  voom(pool.epi.psk$dgelist[flt,i], design = des, plot = T)
v <-  lmFit(v, design = des)
v <-  eBayes(v)
case.tab <-  topTable(v, coef = 2, number = Inf)
case.tab$LumP.overlap <-  rownames(case.tab) %in% lump.markers$lump.overlap
case.tab$gene <-  rownames(case.tab)

pdf("fig3d.class2-epi.volcano.pdf", width = 8, height = 5, useDingbats = F)
ggplot(data = case.tab[!case.tab$LumP.overlap,], aes(logFC, -log10(adj.P.Val)))+geom_point(color ="grey")+
  geom_point(data = case.tab[(case.tab$LumP.overlap), ],color ="purple")+
  geom_text_repel(data = case.tab[(case.tab$LumP.overlap), ], aes(label = gene), max.overlaps = Inf)+
  theme_classic()
dev.off()

write.csv(case.tab, "fig3d.class2-epi.volcano.csv")


################################################
# fig3j
vbls = c("AR.UP","LumP.overlap","Castration20.refine")
df <- pool.epi@meta.data[pool.epi$cluster!=1, c("trt", vbls)] # exclude basal cells and only focus on luminal cells

# center by median of control
df.med <-  aggregate(.~trt, data = df, FUN=median)
df$AR.UP = df$AR.UP - df.med$AR.UP[1]
df$Castration20.refine <-  df$Castration20.refine - df.med$Castration20.refine[1]
df$LumP.overlap <-  df$LumP.overlap - df.med$LumP.overlap[1]
write.csv(df, "fig2j.class2-epi.boxplot.signature.by.trt.csv")

df <-  reshape2::melt(df)
colnames(df)[c(2:3)] <-  c("signature","score")

trt.cls <-  setNames(c("cyan4","coral"), c("ctrl","case"))
pdf("fig3j.class2-epi.boxplot.signature.by.trt.pdf", width = 6, height = 5, useDingbats = F)
ggplot(data = df, aes(trt, score, fill = trt))+geom_boxplot(outlier.size = 0, outlier.color = "white")+
  facet_wrap(~signature, ncol=4)+theme_classic()+#ylim(c(-0.2,1))+
  theme(legend.position = "none")+
  stat_compare_means(label.y = 0.9, cex = 2)+
  scale_fill_manual(values = trt.cls)
dev.off()


#########################
df <-  Embeddings(pool.epi, "umap3")
df <-  cbind(df , pool.epi@meta.data[, c("library",  "trt"  ,"age","cluster",
                                         "pred.crowley.new", "pred.mctp","score.crowley.LumP", 
                                         "score.guo.LumC", "score.karthaus.Lum2", "LumP.overlap")])

write.csv(df, "fig.s8.class2.epitheial.cell.signatures.csv")
