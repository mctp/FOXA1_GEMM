
get.score <-  function(genelist, mat) {
  m = mat[rownames(mat) %in% genelist,]
  m = m[rowSums(m)>0,]
  return(colMeans(t(scale(t(m)))))
}

get.psk <-  function(obj, anno, sample, keep.all.genes = FALSE, min.cell = 5, assay = "RNA") {
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

mouse.human.conv <-  function(mouse.genes = NULL, human.genes = NULL) {
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

stackbar <- function(dat, anno, smp, col=NULL) {
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

sc.proc <-  function(sc, n = 2000, pc =30, res =0.5, normalize = TRUE, hvg = TRUE, pca = TRUE, umap = TRUE, cluster = TRUE, assay = "RNA") {
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


fmt.enrichplot <- function(geneset, setname, rnk, gsea.res, title, padj=TRUE) {
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


custom.vln <-  function(x, dt = pool.epi@meta.data, ylim = FALSE) {
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


volplot <-  function(tab,fc.col='logFC', fdr.col='adj.P.Val', name, add.genes=NULL, lfc=0, circle.genes=NULL, label.genes=NULL) {
  if (!('gene_name' %in% colnames(tab))) {
    tab$gene_name=rownames(tab)
  }
  tab$logFDR= -log10(tab[, fdr.col])
  sig.cls=c(brewer.pal(6,'Reds')[c(3,6)], brewer.pal(6, 'Blues')[c(3,6)])
  names(sig.cls)=c('UP1','UP2','DN1','DN2')
  tab$dir=ifelse(tab[, fc.col]>0, 'UP','DN')
  tab$sig=''
  tab$sig[tab[,fdr.col]<0.05 & abs(tab[,fc.col])>lfc ]=1
  
  tab$dir=paste0(tab$dir, tab$sig)
  tab$dir=factor(tab$dir, levels = c('UP','DN','UP1','DN1'))
  
  if (!is.null(label.genes)) {
    tab2label=label.genes
  } else {
    if (sum(tab[, fc.col][1:30]<0)>=20) {
      tab2label=rbind.data.frame(tab[1:20, ], tab[tab[, fc.col]>lfc,][1:10,])
    } else {
      tab2label=tab[1:30, ]
    }
    
    if (!is.null(add.genes)) {
      
      tab2label=rbind.data.frame(tab2label, tab[tab$gene_name %in% setdiff(add.genes, tab2label$gene_name),])
    }
    tab2label=tab2label[!duplicated(tab2label$gene_name),]
  }
  
  
  a=table(tab$dir)
  lbs=paste0('UP/DN FDR<0.05, logFC=',lfc, ' :', a['UP1'],'/',a['DN1'])
  if (is.null(circle.genes)) {
    g=ggplot(tab)+geom_point(aes_string(fc.col, 'logFDR'), col='grey')+
      geom_point(data=tab[tab$dir %in% c('UP1','DN1'),] , aes_string(fc.col, 'logFDR', color='dir') )+
      #geom_point(data=tab[tab$dir %in% c('UP2','DN2'),] , aes(logFC, -log10(P.Value), color=dir) )+
      scale_color_manual(values=sig.cls)+
      geom_vline(xintercept = c(-lfc, lfc), lty=3)+
      geom_point(data=tab2label ,aes_string(fc.col, 'logFDR'), shape=1)+
      geom_text_repel(data=tab2label ,aes_string(fc.col, 'logFDR', label='gene_name'), max.overlaps=Inf)+
      annotate('text', x=0, y=max(tab$logFDR)+1, label=lbs)+
      theme_classic()+ggtitle(name)
  } else {
    g=ggplot(tab)+geom_point(aes_string(fc.col, 'logFDR'), col='grey')+
      geom_point(data=tab[tab$dir %in% c('UP1','DN1'),] , aes_string(fc.col, 'logFDR', color='dir') )+
      #geom_point(data=tab[tab$dir %in% c('UP2','DN2'),] , aes(logFC, -log10(P.Value), color=dir) )+
      scale_color_manual(values=sig.cls)+
      geom_point(data=tab[tab$gene_name %in% circle.genes,] , aes_string(fc.col, 'logFDR'), shape=1, color='red' )+
      geom_vline(xintercept = c(-lfc, lfc), lty=3)+
      geom_point(data=tab2label ,aes_string(fc.col, 'logFDR'), shape=1)+
      geom_text_repel(data=tab2label ,aes_string(fc.col, 'logFDR', label='gene_name'), max.overlaps=Inf)+
      annotate('text', x=0, y=max(tab$logFDR)+1, label=lbs)+
      theme_classic()+ggtitle(name)
  }
  return(g)
}

motif.volcano <-  function(x, add.genes, filename, ylimits= NULL, w = 7, h =5) {
  x$motif.gene = sapply(x$motif.name, function(x) strsplit(x, "[(]")[[1]][1])
  x$gene = full.g[match(toupper(x$motif.gene), toupper(full.g))]
  x$tf.group = "others"
  x$tf.group[grepl("^Klf", x$motif.gene, ignore.case = T)] = "KLFs"
  x$tf.group[grepl("FOS|JUN|BATF|BACH|AP1", x$motif.gene, ignore.case = T)] = "FOS/JUN/AP1"
  x$tf.group[grepl("FOX", x$motif.gene, ignore.case = T)] = "FOXs" # added later
  
  x$size = log10(x$observed)
  x$fdr = p.adjust(x$pvalue)
  x$logfdr = -log10(x$fdr)
  x$logfdr[x$fdr==0] = max(x$logfdr[x$fdr>0])
  x$sig = x$fold.enrichment>0.5 & x$fdr<0.05
  x.label = x[x$gene %in% union(na.omit(x$gene)[1:30], add.genes), ]
  x.label = x.label[!duplicated(x.label$gene),]
  cls = setNames(c("red","blue","orange","grey"),c("FOS/JUN/AP1","KLFs","FOXs","others"))
  pdf(filename, height = h, width = w, useDingbats = F)
  p = ggplot(data = x[x$tf.group!="FOXs",], aes(fold.enrichment, logfdr))+
    geom_point(aes(color = tf.group, size = percent.observed))+
    scale_color_manual(values = cls )+
    scale_size_continuous(range  = c(0.5,5))+
    geom_point(data = x[x$tf.group=="FOXs",],aes(color = tf.group, size = percent.observed))+
    #geom_point(data = x.label, shape=1)+
    geom_text_repel(data = x.label, aes(label = gene), max.overlaps = Inf)+
    theme_classic()
  if (!is.null(ylimits)) {
    p = p + ylim(ylimits)
  }
  print(p)
  dev.off()
}

gex.atac.dotplot <-  function(g, filename, split = FALSE, col_flt = FALSE, zscore = FALSE, color.limits = FALSE) {
  # raw chromar values have large value outliers, here i set the upper limit and lower limit and 
  #   do average for each cluster 
  #   for plotting, one version is do zscore first then average, the other version is no zscore
  int = intersect(toupper(g$gene_name), toupper(rownames(pool.epi@assays$RNA@data)))
  
  full.g = rownames(pool.epi@assays$RNA@data)
  g$gene_name.mouse = full.g[match(g$gene_name, toupper(full.g))]
  g = g[!is.na(g$gene_name.mouse), ]
  g.expr = pool.epi@assays$RNA@data[g$gene_name.mouse,]
  g.motif = combined.lum@assays$chromvar@data[g$motif, ]
  #g.motif.q10 = apply(g.motif, 1, function(x) quantile(x, 0.1))
  #g.motif.q90 = apply(g.motif, 1, function(x) quantile(x, 0.9))
  g.motif[g.motif>6] = 6
  g.motif[g.motif< -6] = -6
  if (zscore) {
    g.motif = t(scale(t(g.motif)))
  }
  
  if (split) {
    g.expr.ave = get.ave(g.expr, gex.meta$cluster_trt)
    g.motif.ave = get.ave(g.motif, combined.lum$cluster_trt)
    w = 5
  } else {
    g.expr.ave = get.ave(g.expr, gex.meta$cluster)
    g.motif.ave = get.ave(g.motif, combined.lum$cluster)
    w = 2.7
  }
  
  flt = rowSums(g.expr.ave) >0
  g.expr.ave = g.expr.ave[flt,]
  g.motif.ave = g.motif.ave[flt, ]
  rownames(g.motif.ave) = rownames(g.expr.ave) 
  if (split) {
    g.expr.ave = g.expr.ave[, -c(1:2)]
  } else {
    g.expr.ave = g.expr.ave[, -c(1)]
  }
  if (col_flt) {
    cls = c("2_ctrl","3_case","4_case","5_ctrl","6_ctrl","7_ctrl")
    g.expr.ave = g.expr.ave[, cls]
    g.motif.ave = g.motif.ave[, cls]
    w = 2.7
  }
  pdf(filename, width = w, height = 12, useDingbats = F)
  print(myplot2(ma = g.expr.ave, mt = g.motif.ave, motif.scale = FALSE, color.limits = color.limits))
  dev.off()
}