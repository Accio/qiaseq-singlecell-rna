# scRNA-seq secondary data analysis 
# vim: tabstop=9 expandtab shiftwidth=3 softtabstop=3
# Chang Xu, 03OCT2017

# clear all R objects
rm(list=ls())

#list of arguments
args <- commandArgs(TRUE)
wd <- args[1]
umi.counts <- args[2]
ercc.input <- args[3]
qc.metrics <- args[4]
run.id <- args[5]
n.iter <- args[6]
n.cpu <- args[7]
k <- args[8]
perplexity <- args[9]
hvg.thres <- args[10]

# create log file and record starting time
sink(paste0("./", run.id, ".log"))
writeLines("Parameters:")
writeLines(paste0("wd = ", wd))
writeLines(paste0("umi.counts = ", umi.counts))
writeLines(paste0("ercc.input = ", ercc.input))
writeLines(paste0("qc.metrics = ", qc.metrics))
writeLines(paste0("run.id = ", run.id))
writeLines(paste0("n.iter = ", n.iter))
writeLines(paste0("n.cpu = ", n.cpu))
writeLines(paste0("k = ", k))
writeLines(paste0("perplexity = ", perplexity))
writeLines(paste0("hvg.thres = ", hvg.thres, "\n"))

n.iter <- as.numeric(n.iter)
n.cpu <- as.numeric(n.cpu)
k <- as.numeric(k)
perplexity <- as.numeric(perplexity)
hvg.thres <- as.numeric(hvg.thres)

t0 <- Sys.time()
writeLines(paste0("Start time: ", t0))

# MCMC iteration must be multiply of 10 for thining purpose
if(n.iter %% 10 != 0) stop("n.iter must be multiply of 10!")

# load packages
suppressMessages(library(ggplot2))
suppressMessages(library(BASiCS))
suppressMessages(library(dplyr))
suppressMessages(library(Rtsne))
suppressMessages(library(MAST))
suppressMessages(library(scde))
suppressMessages(library(stringr))
suppressMessages(library(ggrepel))
suppressMessages(library(cluster))

options(stringsAsFactors=F)
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename

# mode function
Mode <- function(x) {
   ux <- unique(x)
   ux[which.max(tabulate(match(x, ux)))]
}

# constants 
seed <- 8062017
n.tsne <- 20
writeLines(paste0("Random number seed: ", as.character(seed)))

if(file.exists(wd)){
  setwd(wd)
} else{
  stop(paste0(wd, " not found! Program stopped."))
}
set.seed(seed)

# make output directory
if(!file.exists('output')) dir.create('output')
if(!file.exists('misc')) dir.create('misc')

##############################################################################
# import and process raw data
##############################################################################
# read in raw data
if(file.exists(ercc.input)){
  ercc.orig <- read.csv(ercc.input, header=T, check.names=F)
} else{
  stop(paste0(ercc.input, " not found! Program stopped."))
}

if(file.exists(umi.counts)){
  counts.orig <- read.csv(umi.counts, header=T, check.names=F)
} else{
  stop(paste0(umi.counts, " not found! Program stopped."))
}

if(file.exists(qc.metrics)){
  qc.orig <- read.table(qc.metrics, header=T, sep="\t", check.names=F)
} else{
  stop(paste0(qc.metrics, " not found! Program stopped."))
}

# process umi counts; sum up primer level counts to get gene level counts
counts.orig %>% select(-c(chromosome, start, stop,strand,gene_type)) %>%
            group_by(gene) %>%
            summarise_all(sum) -> dat
n.cell <- ncol(dat) - 1
writeLines(paste0("Number of valid cell indices: ", as.character(n.cell)))

##############################################################################
# prepare datasets required by BASiCS
##############################################################################
# ERCC spike-in input information; drop spike-in's with < 1 expected copy
ercc.orig %>% select(ERCC_ID, input_1ul) %>%
   rename(SpikeID = ERCC_ID, SpikeInput = input_1ul) %>% 
   #dplyr::filter(SpikeInput >= 1) %>% 
   mutate(SpikeInput = round(SpikeInput)) -> SpikeInfo
SpikeInput <- SpikeInfo$SpikeInput

# gene-cell UMI counts matrix
Counts <- as.matrix(dat[,2:ncol(dat)])
class(Counts) <- 'numeric'
rownames(Counts) <- dat$gene

# vector to indicate endogenous genes vs ERCC spike-in
Tech <- ifelse(grepl('ERCC', rownames(Counts)), TRUE, FALSE)

##############################################################################
# quality check - % of UMIs on endogenous genes; RLE plot (for <= 96 cells); observed vs. expected ERCC 
##############################################################################
# plot total UMIs and % of biological gene UMIs for quality control
totalUMIs <- apply(Counts, 2, sum)
totalEndoUMIs <- apply(Counts[!Tech, ], 2, sum)
pctEndoUMIs <- totalEndoUMIs / totalUMIs
pctEndoUMIs.df <- data.frame(pctEndoUMIs)
p.pctEndoUMIs <- ggplot(pctEndoUMIs.df, aes(x=pctEndoUMIs)) + geom_density() + xlab('% of UMIs mapped to endogenous genes')
ggsave(paste0('output/', run.id, '.pct_endoGene_UMI.png'), dpi=300)

# relative log expression (RLE) plot for ERCC; make plot if <= 96 cells. TODO: what to do when n>96? 
if(n.cell <= 96){
   ercc <- Counts[Tech, ]
   medCnt <- apply(ercc, 1, median)
   f <- function(x) log(x[1:length(x)] / x[length(x)])
   logr <- t(apply(ercc[medCnt > 0,], 1, f))
   cnt <- as.vector(logr)
   cell_ID <- rep(colnames(logr), each=nrow(logr))
   data.frame(cnt, cell_ID) %>% mutate(cell_ID = factor(cell_ID)) -> df
   # sort by median RLE
   df.med <- df %>% group_by(cell_ID) %>% filter(cnt != Inf) %>% summarise(med.cnt = median(cnt, na.rm=T)) %>% arrange(desc(med.cnt)) 
   df %>% mutate(cell_ID = factor(cell_ID, levels=df.med$cell_ID)) %>% 
      ggplot(aes(x=cell_ID, y=cnt)) + geom_boxplot() + geom_hline(yintercept=0, linetype='dashed') +
      theme(axis.text.x = element_text(angle = 270, hjust = 1)) + ylab('Relative Log Expression (ERCC)') + xlab('Cell ID') -> p.rle
   ggsave(paste0('output/', run.id, '.RLE_ERCC.png'), dpi=400, height=8.94, width=15)
}

# observed vs expected ERCC; averaged across all cells
#Counts.ercc <- Counts[Tech, ]
#mean.obs.ercc <- apply(Counts.ercc, 1, sum) / n.cell
#df.ercc <- data.frame(mean.obs.ercc, SpikeInput)
#p.scatter.ercc <- ggplot(df.ercc, aes(x=mean.obs.ercc,  y=SpikeInput)) + geom_point() + xlab('observed ERCC UMI per cell (average)') + ylab('expected ERCC UMI per cell')
#ggsave(paste0('output/', run.id, '.observed_vs_expected_UMI_ERCC.png'), dpi=300)

##############################################################################
# filter low quality cells 
##############################################################################
# exclude cells with 1) endogenous gene reads below 5th percentile OR 2) % of endogenous reads below 5th percentile OR 3) detected genes below 5th percentile
qc.orig %>% mutate(Map_to_genes = Mapped_reads - Map_to_ERCC, 
                   pct_Map_to_genes = ifelse(Mapped_reads > 0, Map_to_genes / Mapped_reads, 0), 
                   keep = ifelse(Map_to_genes >= max(100, quantile(Map_to_genes, 0.05)) & 
                                 pct_Map_to_genes >= max(0.1, quantile(pct_Map_to_genes, 0.05)) & 
                                 Detected_genes >= max(1, quantile(Detected_genes, 0.05)), TRUE, FALSE)) -> qc

# write dropped cells to file
qc.drop <- filter(qc, keep == FALSE)
cellsToDrop <- qc.drop$Cell
cellsToKeep <- filter(qc, keep == TRUE)$Cell
write.csv(qc.drop, paste0('output/', run.id, '.cell_dropped.csv'), row.names=F, quote=F)
writeLines(paste0("Low-quality cells dropped: ", as.character(length(cellsToDrop))))
writeLines(paste0("High-quality cells for downstream analyses: ", as.character(length(cellsToKeep))))

# filter very lowly expressed transcripts
totalCountsPerGene <- apply(Counts, 1, sum)
cellsWithExpression <- apply(Counts, 1, function(x) sum(x > 0))
meanExp <- function(x){
   posExp <- sum(x > 0)
   if(posExp == 0){
      meanExpre <- 0
   } else{
      meanExpre <- sum(x) / posExp 
   }
   return(meanExpre)
}

AvCountsPerCellsWithExpression <- apply(Counts, 1, meanExp)

geneInclude <- ifelse(totalCountsPerGene >= 2 & cellsWithExpression >= 2 & AvCountsPerCellsWithExpression >= 2, TRUE, FALSE)
geneDropped <- data.frame(rownames(Counts)[!geneInclude])
colnames(geneDropped) <- 'Genes_dropped'
write.csv(geneDropped, paste0('output/', run.id, '.gene_dropped.csv'), row.names=F, quote=F)
writeLines(paste0("Low-expression genes dropped: ", as.character(nrow(geneDropped))))
writeLines(paste0("Genes for downstream analyses: ", as.character(sum(geneInclude))))

# update data
newCounts <- Counts[geneInclude, cellsToKeep]
newTech <- ifelse(grepl('ERCC', rownames(newCounts)), TRUE, FALSE)
spikeInclude <- rownames(newCounts)[newTech]
newSpikeInfo <- SpikeInfo[SpikeInfo$SpikeID %in% spikeInclude,]
newSpikeInput <- newSpikeInfo$SpikeInput
FilterData <- newBASiCS_Data(Counts=newCounts, Tech=newTech, SpikeInfo=newSpikeInfo)
#save(FilterData, file=paste0('misc/', run.id, '.FilterData.RData'))

##############################################################################
# MCMC sampling to estimate model parameters
##############################################################################
#n.iter <- 400
MCMC_Output <- BASiCS_MCMC(FilterData, N=n.iter, Thin=5, Burn=n.iter/10, StoreChains=T, StoreDir='./misc', RunName=run.id, PrintProgress=F)
save(MCMC_Output, file=paste0('misc/', run.id, '.MCMCoutput.RData'))
#load(paste0('misc/', run.id, '.MCMCoutput.RData'))

# TODO: verify convergence
# TODO: back-up plan if MCMC fails or poor convergence

# visually check convergence afterwards. Need to do it on-the-fly
png(paste0('misc/', run.id, '.traceplot_gene1_cell1.png'), res=300, height=1500, width=1200)
par(mfrow=c(3,2))
plot(MCMC_Output, Param = "mu", Gene = 1, log = "y", cex.lab=1)
plot(MCMC_Output, Param = "delta", Gene = 1, cex.lab=1)
plot(MCMC_Output, Param = "phi", Cell = 1, cex.lab=1)
plot(MCMC_Output, Param = "s", Cell= 1, cex.lab=1)
plot(MCMC_Output, Param = "nu", Cell = 1, cex.lab=1)
plot(MCMC_Output, Param = "theta", Batch = 1, cex.lab=1)
dev.off()

# summarize model fit
MCMC_Summary <- Summary(MCMC_Output)

# variance decomposition
VarDecomp <- BASiCS_VarianceDecomp(FilterData, MCMC_Output)
write.csv(VarDecomp, paste0('output/', run.id, '.Variance_decomposition.csv'), row.names=F, quote=F)

# highly and lowly variable genes
DetectHVG <- BASiCS_DetectHVG(FilterData, MCMC_Output, VarThreshold = hvg.thres, Plot = F)
HVG <- DetectHVG$Table[DetectHVG$Table$HVG, 'GeneNames']
df.HVG <- data.frame(HVG)
write.csv(df.HVG, paste0('output/', run.id, '.highly_variable_genes.csv'), row.names=F, quote=F)

# log-log plot of inter-cell CV vs. mean transcripts per cell
newCounts <- data.frame(counts(FilterData), row.names=rownames(counts(FilterData)))
meanUMI <- apply(newCounts, 1, mean)
CV <- apply(newCounts, 1, function(x) sd(x) / mean(x))
newCounts %>% mutate(meanUMI = meanUMI, CV = CV, 
                     geneType = factor(ifelse(grepl('ERCC', rownames(newCounts)), 'Spike-in', ifelse(rownames(newCounts) %in% HVG, 'HVG', 'Others')), levels=c('HVG', 'Spike-in', 'Others'))) %>% 
            ggplot(aes(x=log(meanUMI), y=log(CV), colour=geneType)) + geom_point() + theme_bw() + xlab('log(mean expression)') + ylab('log(CV)') +
                  geom_abline(slope=-0.5, intercept=0, linetype='dashed') + scale_color_manual(values=c('red', 'blue', 'gray')) -> p.cv_mean
ggsave(paste0('output/', run.id, '.CV_vs_mean_expression.png'), dpi=300)

# normalized UMI counts
DenoisedCounts <- BASiCS_DenoisedCounts(Data=FilterData, Chain=MCMC_Output)
write.csv(DenoisedCounts, paste0('output/', run.id, '.normalized_UMI.csv'), row.names=T, quote=F)

##############################################################################
# cell clustering and visualization; TODO: all genes or HVG only?
##############################################################################
# T-SNE based on highly variable genes
hvg.DenoisedCounts <- t(DenoisedCounts[rownames(DenoisedCounts) %in% HVG, ])
n.pc <- min(c(50,dim(hvg.DenoisedCounts)))
tsne.hvg <- Rtsne(hvg.DenoisedCounts, dims=2, pca=TRUE, perplexity=perplexity, theta=0.0, initial_dims=n.pc)
tsne.hvg.y <- data.frame(tsne.hvg$Y)
p.tsne.hvggenes <- ggplot(tsne.hvg.y, aes(x=X1, y=X2, label=str_sub(rownames(hvg.DenoisedCounts), start=5))) + 
    geom_point() + geom_label_repel() + xlab('TSNE-1') + ylab('TSNE-2') + theme(legend.position="none")
ggsave(paste0('misc/', run.id, '.tsne.highly_variable_genes.png'), dpi=300)

# PCA + Kmeans classification
pca.res <- prcomp(hvg.DenoisedCounts, scale. = T, center=T)
pca.hvg <- pca.res$x[,1:n.pc]

# determine k (number of clusters) if not set by user (i.e. k=0)
if(k == 0){
   # set maximum number of clusters
   k.max <- min(round(nrow(hvg.DenoisedCounts)/10), 15)
   writeLines(paste0("Maximum number of cell populations: ", as.character(k.max)))

   if(k.max >= 2){
      # select optimal number of clusters based on Gap statistic using cluster package
      # Important NOTE: k is chosen based on the first two T-SNE dimensions (multiple T-SNE's with different seeds). This is to ensure consistency between clustering and T-SNE plots  
      k.gap <- rep(NA, n.tsne)
      k.sil <- rep(NA, n.tsne)
      for(j in 1:n.tsne){
         set.seed(j)
         tsne.tmp <- Rtsne(hvg.DenoisedCounts, dims=2, pca=TRUE, perplexity=perplexity, theta=0.0, initial_dims=n.pc)
         gap <- clusGap(tsne.tmp$Y, FUNcluster=kmeans, K.max=k.max, B=500, verbose=F)
         gaps <- gap$Tab[,3]
         se <- gap$Tab[,4]
         k.gap[j] <- maxSE(gaps, se, method='Tibs2001SEmax')

         # select optimal number of clusters based on Silhoutte statistic 
         dists <- dist(tsne.tmp$Y)
         sil <- rep(0, k.max)
         for(i in 2:k.max){
          km.res <- kmeans(tsne.tmp$Y, centers = i, nstart = 5, iter.max=100)
          ss <- silhouette(km.res$cluster, dist = dists)
          sil[i] <- mean(ss[,3])
         }
         k.sil[j] <- which.max(sil)
      }
      k.gap.mode <- Mode(k.gap)
      k.sil.mode <- Mode(k.sil)
      k <- min(k.gap.mode, k.sil.mode)
      writeLines(paste0("Cell populations based on Gap statistics: ", as.character(k.gap.mode)))
      writeLines(paste0("Cell populations based on Silhoutte method: ", as.character(k.sil.mode)))
      writeLines(paste0("final cell populations (minimum of K_gap and K_silhoutte): ", as.character(k)))    
      set.seed(seed)
  } 
}

# if k = 1, 
if(k == 1){
   writeLines("Clustering and differential expression are not performed because there is only one cluster.")
}

if(k >= 2){
   # k-means clustering
   cluster.kmeans.hvg <- kmeans(pca.hvg, centers=k, nstart=5, iter.max=100)
   tsne.hvg.y$kmeans <- cluster.kmeans.hvg$cluster
   
   # t-sne plot with k-means clustering
   p.tsne.hvg.kmeans <- ggplot(tsne.hvg.y, aes(x=X1, y=X2, colour=factor(kmeans))) + 
      geom_point() + xlab('TSNE-1') + ylab('TSNE-2') + theme(legend.position="bottom", legend.title = element_blank()) + guides(color = FALSE)
   ggsave(paste0('output/', run.id, '.kmeans.hvg.png'), dpi=300)
   
   # save clustering results to disk
   cell <- rownames(hvg.DenoisedCounts)
   cluster <- factor(cluster.kmeans.hvg$cluster)
   TSNE_1 <- tsne.hvg.y$X1
   TSNE_2 <- tsne.hvg.y$X2
   kmeans.clustering <- data.frame(cell, cluster, TSNE_1, TSNE_2)
   write.csv(kmeans.clustering, paste0('output/', run.id, '.kmeans.clustering.csv'), row.names = FALSE, quote = FALSE)
   
   ##############################################################################
   # Differential expression between cell clusters (based on HVG) using scde
   ##############################################################################
   df.hvg.DenoisedCounts <- data.frame(t(hvg.DenoisedCounts), check.names=F)
   df.hvg.DenoisedCounts <- apply(df.hvg.DenoisedCounts,2,function(x) {storage.mode(x) <- 'integer'; x}) 

   # run pair-wise differential expression analysis on HVG
   final.de <- NULL
   pairs <- combn(sort(unique(cluster)), 2)  # all pairs of clusters; each col is a pair
   for(i in 1:ncol(pairs)){
      pair <- pairs[,i]
      writeLines(paste0("Differential expression analysis of group ", as.character(pair[1]), " and ", as.character(pair[2])))
      # select cells of current 2 clusters being compared
      df.hvg.pair <- df.hvg.DenoisedCounts[, names(cluster[cluster %in% pair])]
      # fit individual error model for each cell. 
      o.ifm.hvg <- scde.error.models(counts = df.hvg.pair, n.cores = n.cpu, min.size.entries = 50, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 0)
      # remove cells with abnormal fits (corr.a <= 0); all valid in this data
      valid.cells.hvg <- o.ifm.hvg$corr.a > 0
      o.ifm.hvg <- o.ifm.hvg[valid.cells.hvg,]
      # estimate gene expression prior
      o.prior.hvg <- scde.expression.prior(models = o.ifm.hvg, counts = df.hvg.pair, length.out = 400, show.plot = FALSE)
      # re-define clusters on o.ifm.hvg because some cells may be invalid
      cluster.new.hvg <- factor(cluster[names(cluster) %in% rownames(o.ifm.hvg) & cluster %in% pair])
      ediff.hvg <- scde.expression.difference(o.ifm.hvg, df.hvg.pair, o.prior.hvg, groups  =  cluster.new.hvg, n.randomizations  =  100, n.cores  =  n.cpu, verbose  =  0)
      # add to final differencial expression results
      ediff.hvg$Group_A <- pair[1]
      ediff.hvg$Group_B <- pair[2]
      final.de <- rbind(final.de, ediff.hvg)
   }
   # save all results to disk
   final.de %>% mutate(absZ = abs(Z)) %>% arrange(Group_A, Group_B, desc(absZ)) %>% select(c(Group_A, Group_B, lb, mle, ub, ce, Z, cZ)) -> final.de
   write.csv(final.de, file = paste0('output/', run.id, '.diff.exp.hvg.csv'), row.names = FALSE, quote = FALSE)
}

# record end time 
writeLines(paste0("End time: ", Sys.time()))
# print session info
sessionInfo()
sink()


