# this R code has been condensed from the tutorial here: 
	# github.com/CBC-UCONN/RNA-seq-with-reference-genome-and-annotation

# Load the libraries we'll need in the following code:
library("DESeq2")
library("apeglm")
library("pheatmap")
library("tidyverse")
library("dplyr")
library("ashr")
library("goseq")
library("biomaRt")
library("ggplot2")
library("ggrepel")
library("clusterProfiler")
library("enrichplot")
library("ggupset")

######################################################
# Point R to count data, set an output file prefix 
######################################################

# create an object with the directory containing your counts:
	# !!edit this to point to your own count file directory!!
directory <- "06_count/counts"

# ensure the count files are where you think they are
list.files(directory)

sampleFiles <- list.files(directory, pattern = ".*counts$")


######################################################
# Read the count data into R along with treatment information
######################################################

# load in metadata table
meta <- read.csv("01_raw_data/metadata.txt") %>%
	mutate(population=str_replace(population, "Eliza.*", "ER")) %>%
	mutate(population=str_replace(population, "King.*", "KC")) %>%
	mutate(pcb_dosage = case_when(pcb_dosage == 0 ~ "control", pcb_dosage > 0 ~ "exposed"))

# ensure that sampleFiles and metadata table are in the same order

all( str_remove(sampleFiles, ".counts") == meta[,1] )

# now create a data frame with sample names, file names and treatment information. 
sampleTable <- data.frame(
	sampleName = meta$Run,
	fileName = sampleFiles,
	population = meta$population,
	dose = meta$pcb_dosage
	)

# look at the data frame to ensure it is what you expect:
sampleTable

# create the DESeq data object
ddsHTSeq <- DESeqDataSetFromHTSeqCount(
		sampleTable = sampleTable, 
		directory = directory, 
		design = ~ population*dose
		)

######################################################
# Reset treatment factors
######################################################

# factor order determines reference level (i.e. control vs exposed)
	# factors are automatically ordered alphabetically
	# we want KC as the reference level

# To see the levels as they are now:
ddsHTSeq$population

# To replace the order with one of your choosing, create a vector with the order you want:
treatments <- c("KC","ER")

# Then reset the factor levels:
ddsHTSeq$population <- factor(ddsHTSeq$population, levels = treatments)

# verify the order
ddsHTSeq$population

######################################################
# Filter out genes with very low expression
######################################################

# what does expression look like across genes?

# sum counts for each gene across samples
sumcounts <- rowSums(counts(ddsHTSeq))
# take the log
logsumcounts <- log(sumcounts,base=10)
# plot a histogram of the log scaled counts
hist(logsumcounts,breaks=100)

# you can see the typically high dynamic range of RNA-Seq, with a mode in the distribution around 1000 fragments per gene, but some genes up over 1 million fragments. 

# get genes with summed counts greater than 20
keep <- sumcounts > 50

# keep only the genes for which the vector "keep" is TRUE
ddsHTSeq <- ddsHTSeq[keep,]

######################################################
# Run the statistical analysis
######################################################

dds <- DESeq(ddsHTSeq)

######################################################
# Plot dispersion estimates
######################################################
plotDispEsts(dds)


######################################################
# Which results do we want?
######################################################

# we are interested in:
	# 1) genes responding to exposure in KC, the sensitive population
	# 2) genes responding to exposure in ER, the tolerant population
	# 3) genes that respond differently to exposure in ER and KC

# how do we extract that?

# in this experiment there are four coefficents. you can see them here:

resultsNames(dds)

# we can get (1) by the two below methods. 
	# they extract ONLY the effect of treatment in the reference population. 

res1 <- results(dds, contrast=c("dose","exposed","control"))

# or

res1b <- results(dds, name="dose_exposed_vs_control")

# we can get (2) by:

res2 <- results(dds, contrast=list(c("dose_exposed_vs_control","populationER.doseexposed")))

# we can get (3) by:

res3 <- results(dds, name="populationER.doseexposed")

# just for fun, we can also ask about baseline differences between ER and KC

res4 <- results(dds, name="population_ER_vs_KC")


######################################################
# Quickly summarize results
######################################################

# get a quick summary of the table
summary(res1)
summary(res1b)
summary(res2) # note that independent filtering goes haywire here b/c very few significant genes
summary(res3)
summary(res4)

# check out the first few lines
head(res1)

# or sort the table by pvalue and get the first n lines

as.data.frame(res1) %>% arrange(pvalue) %>% head(n=10)


######################################################
# Get a table of shrunken log2 fold changes
######################################################

# get shrunken log fold changes
res_shrink1 <- lfcShrink(dds,type="ashr",coef="dose_exposed_vs_control")
res_shrink2 <- lfcShrink(dds,type="ashr",contrast=list(c("dose_exposed_vs_control","populationER.doseexposed")))
res_shrink3 <- lfcShrink(dds,type="ashr",coef="populationER.doseexposed")
res_shrink4 <- lfcShrink(dds,type="ashr",coef="population_ER_vs_KC")

# plot the shrunken log2 fold changes against the raw changes:

	data.frame(l2fc=res1$log2FoldChange, l2fc_shrink=res_shrink1$log2FoldChange, padj=res1$padj) %>%
		filter(l2fc > -5 & l2fc < 5 & l2fc_shrink > -5 & l2fc_shrink < 5) %>%
		ggplot(aes(x=l2fc, y=l2fc_shrink,color=padj > 0.1)) +
			geom_point(size=.25) + 
			geom_abline(intercept=0,slope=1, color="gray")


# get the top 20 genes by shrunken log2 fold change
	# this will include genes with outliers
arrange(data.frame(res_shrink1), -abs(log2FoldChange)) %>% head(., n=20)


######################################################
# Data visualization
######################################################

# MA plot
plotMA(res1, ylim=c(-4,4))
plotMA(res_shrink1, ylim=c(-4,4))

##############

#Volcano plot

# negative log-scaled adjusted p-values
log_padj <- -log(res_shrink1$padj,10)
log_padj[log_padj > 100] <- 100

# plot
plot(x=res_shrink1$log2FoldChange,
     y=log_padj,
     pch=20,
     cex=.5,
     col=(log_padj > 10)+1, # color padj < 0.1 red
     ylab="negative log-scaled adjusted p-value",
     xlab="shrunken log2 fold changes")

#############

# PCA plot

# normalized, variance-stabilized transformed counts for visualization
vsd <- vst(dds, blind=FALSE)

plotPCA(vsd, intgroup=c("population","dose"))

# alternatively, using ggplot

dat <- plotPCA(vsd,returnData=TRUE,intgroup=c("population","dose"))

p <- ggplot(dat,aes(x=PC1,y=PC2,col=paste(population, dose)))
p <- p + geom_point() + 
	xlab(paste("PC1: ", round(attr(dat,"percentVar")[1],2)*100, "% variation explained", sep="")) + 
	ylab(paste("PC2: ", round(attr(dat,"percentVar")[2],2)*100, "% variation explained", sep="")) +
	geom_label_repel(aes(label=name))
p

##############

##Just KC dose response
vsd_kcdose <- vsd[which(res1$padj < 0.1),]

dat <- plotPCA(vsd_kcdose,returnData=TRUE,intgroup=c("population","dose"))

p <- ggplot(dat,aes(x=PC1,y=PC2,col=population, shape=dose))
p <- p + geom_point() + 
  xlab(paste("PC1: ", round(attr(dat,"percentVar")[1],2)*100, "% variation explained", sep="")) + 
  ylab(paste("PC2: ", round(attr(dat,"percentVar")[2],2)*100, "% variation explained", sep="")) +
  geom_label_repel(aes(label=name))
p

# heatmap of DE genes in Results 1 - 
# 1) genes responding to exposure in KC, the sensitive population

# create a metadata data frame to add to the heatmap
df <- data.frame(colData(dds)[,c("population","dose")])
rownames(df) <- colnames(dds)
colnames(df) <- c("population","dose")

# pull the top 30 genes by shrunken log2 fold change
top30_KCDose <- res_shrink1 %>% 
  data.frame() %>%
  filter(padj < 0.1) %>%
  arrange(-abs(log2FoldChange)) %>%
  rownames %>%
  head(n=30)

# make the heatmap
pheatmap(
  assay(vsd)[top30_KCDose,], 
  cluster_rows=TRUE, 
  show_rownames=TRUE,
  cluster_cols=TRUE,
  annotation_col=df
)

##heatmap with pearson correlation
pheatmap(
  assay(vsd)[top30_KCDose,], 
  cluster_rows=TRUE, 
  show_rownames=TRUE,
  cluster_cols=TRUE,
  annotation_col=df,
  clustering_distance_rows="correlation"
)
##heatmap with values rescaled by base mean
rescaled <- assay(vsd) - rowMeans(assay(vsd))

pheatmap(
  rescaled[top30_KCDose,], 
  cluster_rows=TRUE, 
  show_rownames=TRUE,
  cluster_cols=TRUE,
  annotation_col=df,
  clustering_distance_rows="correlation"
)

###heatmap for results 3: using regularized log scaled counts
  #genes that respond differently to exposure in ER and KC

# regularized log transformation of counts
rld <- rlog(dds, blind=FALSE)

# order gene names by absolute value of shrunken log2 fold change (excluding cook's cutoff outliers)
lfcorder <- data.frame(res_shrink3) %>%
  filter(!is.na(padj)) %>% 
  arrange(-abs(log2FoldChange)) %>% 
  rownames() 

# use regularized log-scaled counts
pheatmap(
  assay(rld)[lfcorder[1:30],], 
  cluster_rows=TRUE, 
  show_rownames=TRUE,
  cluster_cols=TRUE,
  annotation_col=df
)

# re-scale regularized log-scaled counts by baseMean (estimated mean across all samples)
pheatmap(
  assay(rld)[lfcorder[1:30],] - log(res1[lfcorder[1:30],"baseMean"],2), 
  cluster_rows=TRUE, 
  show_rownames=TRUE,
  cluster_cols=TRUE,
  annotation_col=df
)

# re-scale regularized log-scaled counts by reference level (KC control samples)
pheatmap(
  assay(rld)[lfcorder[1:30],] - rowMeans(assay(rld)[lfcorder[1:30],dds$population=="KC" & dds$dose=="control"]), 
  cluster_rows=TRUE, 
  show_rownames=TRUE,
  cluster_cols=TRUE,
  annotation_col=df
)

##############

# plot counts for individual genes

plotCounts(dds, gene=lfcorder[9], intgroup=c("population","dose"))

##or take a look at results 1:
as.data.frame(res1) %>% arrange(pvalue) %>% head(n=10)

plotCounts(dds, gene="ENSFHEG00000008855", intgroup=c("population","dose"))
plotCounts(dds, gene="ENSFHEG00000021655", intgroup=c("population","dose"))

######################################################
# Get gene annotations using biomaRt
######################################################


##############################
# Select a mart and dataset
##############################

# ensembl host:
  # most recent is "https://ensembl.org"
  # to list archived version hosts: listEnsemblArchives()

ensemblhost <- "https://oct2024.archive.ensembl.org"

listMarts(host=ensemblhost)

# create an object for the Ensembl Genes v100 mart
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host=ensemblhost)

# occasionally ensembl will have connectivity issues. we can try an alternative function:
	# select a mirror: 'www', 'uswest', 'useast', 'asia'
	# mart <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", mirror = "useast")

# see a list of datasets within the mart
	# at the time of writing, there were 203
listDatasets(mart)

# figure out which dataset is the killifish
	# be careful using grep like this. verify the match is what you want
searchDatasets(mart,pattern="fheteroclitus_gene_ensembl")

# there's only one match, get the name
killidata <- searchDatasets(mart,pattern="Mummichog")[,1]

# create an object for the killifish dataset
killi_mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = ensemblhost, dataset = killidata)

# if above there were connectivity issues and you used the alternative function then:
	# select a mirror: 'www', 'uswest', 'useast', 'asia'
	# killi_mart <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = killidata, mirror = "useast")


#########################
# Query the mart/dataset
#########################

# filters, attributes and values

# see a list of all "filters" available for the mummichog dataset.
	# at the time of writing, over 300
listFilters(killi_mart)

# see a list of all "attributes" available
	# 130 available at the time of writing
listAttributes(mart = killi_mart, page="feature_page")

# we can also search the attributes and filters
searchAttributes(mart = killi_mart, pattern = "ensembl_gene_id")

searchFilters(mart = killi_mart, pattern="ensembl")

# get gene names and transcript lengths when they exist
ann <- getBM(filter="ensembl_gene_id",value=rownames(res1),attributes=c("ensembl_gene_id","description","transcript_length","external_gene_name"),mart=killi_mart)

# pick only the longest transcript for each gene ID
ann <- group_by(ann, ensembl_gene_id) %>% 
  summarize(.,description=unique(description),external_gene_name=unique(external_gene_name),transcript_length=max(transcript_length)) %>%
  as.data.frame()

# get GO term info
  # each row is a single gene ID to GO ID mapping, so the table has many more rows than genes in the analysis
go_ann <- getBM(filter="ensembl_gene_id",value=rownames(res1),attributes=c("ensembl_gene_id","description","go_id","name_1006","definition_1006","namespace_1003"),mart=killi_mart)

#let's look at the annotation:
filter(ann, description!="") %>% head()
head(ann)
# and Go term:
head(go_ann)
# put results and annotation in the same table
as.data.frame(res1)

res1_data <- res1 %>% 
  as.data.frame() %>%
  rownames_to_column(var = "ensembl_gene_id") 

res1_ann <- merge(x=res1_data,y=ann, by.x="ensembl_gene_id", by.y="ensembl_gene_id", all.x=TRUE)

## Over-representation and GSEA with clusterProfiler

# get ENSEMBL gene IDs for genes with padj < 0.1
genes <- rownames(res1[which(res1$padj < 0.1),])
# get ENSEMBL gene IDs for universe (all genes with non-NA padj passed independent filtering)
univ <- rownames(res1[!is.na(res1$padj),])
# pull out the columns of go_ann containing GO IDs and descriptions, keep only unique entries. 
go_ann[,c(3,4)] %>% head()
gonames <- unique(go_ann[,c(3,4)]) %>% unique()

res1_enrich <- enricher(
  gene=genes,
  universe=univ,
  TERM2GENE=go_ann[,c(3,1)],
  TERM2NAME=gonames
)

as.data.frame(res1_enrich) %>% head()

# GeneRatio: number of DE genes in category / total DE genes
# BgRation: total genes in category / total genes considered
# pvalue: raw p-value
# p.adjust fdr adjusted p-value (BH method)
# qvalue: alternative fdr adjusted p-value

#### Gene-set enrichment analysis
  # determine if functional categories tend to be enriched for higher or lower expression levels
  #shrunken or raw log2fc? depends on preference for smoothing noise but reducing signal

# extract log2 fold changes, only for genes that passed independent filtering. 
l2fcs <- as.data.frame(res1) %>% 
  filter(!is.na(padj)) %>%
  select(log2FoldChange)
View(l2fcs)

# put log2FCs in a vector, add gene IDs as names, sort 
l2fcvec <- l2fcs[,1]
names(l2fcvec) <- rownames(l2fcs)
View(l2fcvec)
l2fcvec <- sort(l2fcvec, decreasing=TRUE)
View(l2fcvec)

res1_gsea <- GSEA(
  geneList=l2fcvec,
  TERM2GENE=go_ann[,c(3,1)],
  TERM2NAME=gonames
) 


##Upregualted term: 
gseaplot(res1_gsea, by = "all", title = res1_enrich$Description[2], geneSetID = 2)
gseaplot(res1_gsea, by = "all", title = res1_enrich$Description[1], geneSetID = 1)

## GeneRatio in GSEA results: fraction of genes found in leading edge subset

###visuals
dotplot(res1_enrich)
dotplot(res1_gsea)
upsetplot(res1_enrich)
upsetplot(res1_gsea)
ridgeplot(res1_gsea,label_format=30)



