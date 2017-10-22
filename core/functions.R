# functions for scRNA-seq secondary data analysis 
# vim: tabstop=9 expandtab shiftwidth=3 softtabstop=3
# Chang Xu, 11OCT2017

# mode function
Mode <- function(x) {
   ux <- unique(x)
   ux[which.max(tabulate(match(x, ux)))]
}

# mean expression function
meanExp <- function(x){
  posExp <- sum(x > 0)
  if(posExp == 0){
    meanExpre <- 0
  } else{
    meanExpre <- sum(x) / posExp 
  }
  return(meanExpre)
}

# quality check for genes
gene.check <- function(counts.table){
  totalCountsPerGene <- apply(counts.table, 1, sum)
  cellsWithExpression <- apply(counts.table, 1, function(x) sum(x > 0))
  AvCountsPerCellsWithExpression <- apply(counts.table, 1, meanExp)
  geneInclude <- ifelse(totalCountsPerGene >= 5 & cellsWithExpression >= 5 & AvCountsPerCellsWithExpression >= 5, TRUE, FALSE)
  if(all(geneInclude)){
    qc.pass <- TRUE
    geneDropped <- c()
    new.counts <- counts.table
  } else{
    qc.pass <- FALSE
    geneDropped <- rownames(counts.table)[!geneInclude]
    new.counts <- counts.table[geneInclude,]
  }
  out <- list(qc.pass = qc.pass, gene.drop = geneDropped, new.counts = new.counts)
  return(out)
}

# quality check for cells
cell.check <- function(counts.table){
  totalCountsPerCell <- apply(counts.table, 2, sum)
  cellInclude <- ifelse(totalCountsPerCell >= 5, TRUE, FALSE)
  if(all(cellInclude)){
    qc.pass <- TRUE
    cellDropped <- c()
    new.counts <- counts.table
  } else{
    qc.pass <- FALSE
    cellDropped <- colnames(counts.table)[!cellInclude]
    new.counts <- counts.table[,!cellInclude,]
  }
  out <- list(qc.pass = qc.pass, cell.drop = cellDropped, new.counts = new.counts)
  return(out)
}

# overall quality check - alternatively checks row and col until both meet minimum requirement
overall.check <- function(counts.table){
  gene.pass <- FALSE
  cell.pass <- FALSE
  gene.drop <- c()
  cell.drop <- c()
  new.table <- counts.table
  while(!gene.pass | !cell.pass){
    # check genes first
    geneCheck <- gene.check(new.table)
    new.table <- geneCheck$new.counts
    gene.drop <- c(gene.drop, geneCheck$gene.drop)
    gene.pass <- geneCheck$qc.pass
    
    # check cells second
    cellCheck <- cell.check(new.table)
    new.table <- cellCheck$new.counts
    cell.drop <- c(cell.drop, cellCheck$cell.drop)
    cell.pass <- cellCheck$qc.pass
  }
  out <- list(new.table = new.table, gene.drop = gene.drop, cell.drop = cell.drop)
  return(out)
}

# cell clustering
kmeans.cluster <- function(top.pc, tsne, nclust, subdir, run.id, col){
   # k-means clustering
   cluster.kmeans <- kmeans(top.pc, centers=nclust, nstart=5, iter.max=100)
   tsne$kmeans <- cluster.kmeans$cluster
   # t-sne plot with k-means clustering
   p.tsne.kmeans <- ggplot(tsne, aes(x=X1, y=X2, colour=factor(kmeans))) + geom_point() + scale_colour_manual(values = col) + 
      xlab('TSNE-1') + ylab('TSNE-2') + theme(legend.position="bottom", legend.title = element_blank()) 
   ggsave(paste0(subdir, "/", run.id, ".kmeans.clustering.png"), dpi=300, height=7, width=7)
   # save clustering results to disk
   cell <- rownames(top.pc)
   cluster <- factor(cluster.kmeans$cluster)
   TSNE_1 <- tsne$X1
   TSNE_2 <- tsne$X2
   cluster.res <- data.frame(cell, cluster, TSNE_1, TSNE_2)
   write.csv(cluster.res, paste0(subdir, "/", run.id, ".kmeans.clustering.csv"), row.names = FALSE, quote = FALSE)
   # return cluster results (group membership)
   return(cluster)
}

# select the optimal K using Gap statistics
find.k <- function(top.pc, kmax, n.trial, n.boot) {
   k.gap <- rep(NA, n.trial)
   for(i in 1:n.trial){
      gap.res <- clusGap(top.pc, FUNcluster=kmeans, K.max=kmax, B=n.boot, verbose=F)$Tab
      gap <- gap.res[,3]
      se <- gap.res[,4]
      k.gap[i] <- maxSE(gap, se, method='Tibs2001SEmax')
   }
   k.optimal <- Mode(k.gap)
   return(k.optimal)
}
   

# differential expression analysis; norm.table: normalized umi counts, genes in rows and cells in cols. counts must be in integer; clusters: vector of clusters in factor
pair.de <- function(o.ifm, o.prior, norm.table, cluster, pair, nclust, n.cpu){
   ediff <- scde.expression.difference(o.ifm, norm.table, o.prior, groups  =  cluster, n.randomizations  =  100, n.cores  =  n.cpu, verbose  =  0)
   ediff %>% mutate(K = nclust, Group_A = pair[1], Group_B = pair[2], p = 2 * pnorm(-abs(Z)), adj_p = 2 * pnorm(-abs(cZ)), abs.cz = abs(cZ), Gene = rownames(ediff)) %>%
             rename(Fold_change = mle, Lower_limit = lb, Upper_limit = ub, Z_score = Z, Adjusted_z_score =cZ, P_value = p, Adjusted_p_value = adj_p) %>%
             select(c(K, Group_A, Group_B, Gene, Fold_change, Lower_limit, Upper_limit, Z_score, P_value, Adjusted_z_score, Adjusted_p_value)) %>% 
             arrange(Adjusted_p_value) -> ediff
   return(ediff)
}

# plot heatmap 
heatmap_wrapper <- function(counts, cluster, filename, col){
   # sort cells by cluster; within each cluster, sort by median expression
   med <- apply(counts, 1, median)
   group <- cluster$cluster
   counts.sorted <- rownames_to_column(counts, "cell") %>% mutate(med = med, group = group) %>% arrange(group, med) %>% column_to_rownames("cell") %>% select(-c(group, med))
   # annotation colors
   names(col) <- levels(group)
   ann_colors <- list(cluster = col)
   # make plot
   pheatmap(mat = counts.sorted, cluster_cols=F, cluster_rows=F, show_rownames=T, show_colnames=T, annotation_row=cluster, filename = filename, fontsize_row=3, fontsize_col=3, 
            annotation_colors = ann_colors)
}



