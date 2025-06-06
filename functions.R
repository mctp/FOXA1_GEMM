
get.score = function(genelist, mat) {
  m = mat[rownames(mat) %in% genelist,]
  m = m[rowSums(m)>0,]
  return(colMeans(t(scale(t(m)))))
}

get.psk = function(obj, anno, sample, keep.all.genes = FALSE, min.cell = 5, assay = "RNA") {
  require(edgeR)
  obj$anno_sample = paste(as.character(obj@meta.data[,anno]), obj@meta.data[,sample], sep = ".")
  a = table(obj$anno_sample)
  a = a[a>min.cell]
  obj = obj[, obj$anno_sample %in% names(a)]
  whole.psk = Seurat::AggregateExpression(object = obj, assays = assay, return.seurat = F, group.by = "anno_sample")
  whole.psk = whole.psk[[1]]
  colnames(whole.psk) = names(a)
  if (keep.all.genes) {
    whole.psk = DGEList(counts = whole.psk)
  } else {
    whole.psk = DGEList(counts = whole.psk[rowSums(whole.psk)>0, ])
  }
  
  whole.psk = calcNormFactors(whole.psk)
  logtpm = log2(edgeR::cpm(whole.psk)+1)
  s = sapply(colnames(whole.psk), function(x) strsplit(x, "[.]")[[1]])
  s = as.data.frame(t(s))
  colnames(s) = c("anno", "sample")
  whole.psk$samples = cbind.data.frame(whole.psk$samples, s)
  return(list(dgelist = whole.psk, logtpm = logtpm))
}

mouse.human.conv = function(mouse.genes = NULL, human.genes = NULL) {
  require(biomaRt)
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
  if (!is.null(mouse.genes)) {
    m2h=getLDS(attributes = c("mgi_symbol"), 
               filters = "mgi_symbol", values = mouse.genes, mart = mouse, 
               attributesL = c("hgnc_symbol"), martL = human , uniqueRows = TRUE) 
  } else {
    m2h=getLDS(attributes = c("hgnc_symbol"), 
               filters = "hgnc_symbol", values = human.genes, mart = human, 
               attributesL = c("mgi_symbol"), martL = mouse , uniqueRows = TRUE)
  }
  return(m2h)
}

stackbar=function(dat, anno, smp, col=NULL) {
  d=table(dat[, anno], dat[, smp])
  d=as.data.frame.matrix(d)
  d=100*t(t(d)/colSums(d))
  tmp=reshape2::melt(d)
  colnames(tmp)=c('Cell','Sample','Percentage')
  if (is.null(col)) {
    g=ggplot(tmp, aes(Sample, Percentage, fill=Cell))+geom_bar(stat = 'identity', position = 'stack', color='black')+
      theme_classic()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
  } else {
    g=ggplot(tmp, aes(Sample, Percentage, fill=Cell))+geom_bar(stat = 'identity', position = 'stack', color='black')+
      scale_fill_manual(values = col)+
      theme_classic()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
  }
  return(g)
}

sc.proc = function(sc, n = 2000, pc =30, res =0.5, normalize = TRUE, hvg = TRUE, pca = TRUE, umap = TRUE, cluster = TRUE, assay = "RNA") {
  DefaultAssay(sc) = assay
  if (normalize) {
    sc <- NormalizeData(sc, verbose = FALSE)
  } 
  if (hvg) {
    sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = n, verbose = FALSE)
  }
  if (pca) {
    sc <- ScaleData(sc, verbose = FALSE)
    sc <- RunPCA(sc, npcs = pc, verbose = FALSE)
  }
  if (umap) {
    sc <- RunUMAP(sc, reduction = "pca", dims = 1:pc)
  }
  #sc <- RunTSNE(sc, reduction = "pca", dims = 1:pc)
  if (cluster) {
    sc <- FindNeighbors(sc, dims=1:pc)
    sc <- FindClusters(sc, resolution = res)
    sc <- cluster.reorder(sc)
  }
  return(sc)
}


fmt.enrichplot=function(geneset, setname, rnk, gsea.res, title, padj=TRUE) {
  require(fgsea)
  require(gridExtra)
  d=gsea.res[gsea.res$pathway==setname]
  if (d$ES<0) {p=d$ES+0.2} else {p=d$ES-0.1}
  label.pos=floor(0.8*length(rnk))
  
  if (!padj) {
    plotEnrichment(geneset[[setname]], rnk)+labs(title= title)+
      annotate("text",
               label=paste0("NES = ", round(d$NES,2)), 
               x=label.pos, y=p )+
      annotate("text",
               label=paste0("pval = ",formatC(d$pval, format = "e", digits = 2)  ), #round(d$pval,3)
               x=label.pos, y=p-0.1 )+
      theme_minimal(base_size = 10)+ theme(plot.title = element_text(size=10))
  } else {
    plotEnrichment(geneset[[setname]], rnk)+labs(title= title)+
      annotate("text",
               label=paste0("NES = ", round(d$NES,2)), 
               x=label.pos, y=p )+
      annotate("text",
               label=paste0("p.adj = ", formatC(d$pval, format = "e", digits = 2) ), #round(d$padj,3)
               x=label.pos, y=p-0.1 )+
      theme_minimal(base_size = 10)+ theme(plot.title = element_text(size=10))
  }
}


custom.vln = function(x, dt = pool.epi@meta.data, ylim = FALSE) {
  vln.cols = brewer.pal(7, "RdBu")
  d = cluster.mean[order(-cluster.mean[,x]),]
  col.val = setNames(vln.cols, d$cluster)
  #dt = pool.epi@meta.data
  dt$cluster = as.factor(dt$cluster)
  colnames(dt)[colnames(dt)==x]= "score"
  if (!(ylim)) {
    ggplot(data = dt, aes(cluster, score, fill = cluster))+
      geom_violin(scale = "width")+
      scale_fill_manual(values = col.val)+ggtitle(x)+
      theme_bw()+theme()
  } else {
    ggplot(data = dt, aes(cluster, score, fill = cluster))+
      geom_violin(scale = "width")+ylim(c(0, max(dt$score)+1))+
      scale_fill_manual(values = col.val)+ggtitle(x)+
      theme_bw()+theme()
  }
}
