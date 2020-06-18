# this R code has been condensed from the tutorial here: 
	# https://github.com/CBC-UCONN/RNA-seq-with-reference-genome-and-annotation


# Load DESeq2 library
library("DESeq2")


######################################################
# Point R to count data, set an output file prefix 
######################################################

# create an object with the directory containing your counts:
	# !!edit this to point to your own count file directory!!
directory <- "../count"

# ensure the count files are where you think they are
list.files(directory)

sampleFiles <- list.files(directory, pattern = ".*counts_2", full.names = TRUE)


######################################################
# Read the count data into R along with treatment information
######################################################

# create a vector of sample names. ensure these are in the same order as the "sampleFiles" object!
sampleNames <- c("LB2A_1","LB2A_2","LC2A_1","LC2A_2")

# create a vector of conditions. again, mind that they are ordered correctly!
sampleCondition <- c("control","control","treated","treated")

# now create a data frame from these three vectors. 
sampleTable <- data.frame(
		sampleName = sampleNames,
		fileName = sampleFiles,
		condition = sampleCondition
		)

# look at the data frame to ensure it is what you expect:
sampleTable

# create the DESeq data object
ddsHTSeq <- DESeqDataSetFromHTSeqCount(
		sampleTable = sampleTable, 
		directory = directory, 
		design = ~ condition
		)

######################################################
# Reset treatment factors
######################################################

# To see the levels as they are now:
ddsHTSeq$condition

# To replace the order with one of your choosing, create a vector with the order you want:
treatments <- c("control","treated")

# Then reset the factor levels:
ddsHTSeq$condition <- factor(ddsHTSeq$condition, levels = treatments)

# verify the order
ddsHTSeq$condition

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
keep <- sumcounts > 20

# keep only the genes for which the vector "keep" is TRUE
ddsHTSeq <- ddsHTSeq[keep,]

######################################################
# Run the statistical analysis
######################################################

dds <- DESeq(ddsHTSeq)

######################################################
# Get a table of results
######################################################

# get results table
res <- results(dds)

# get a quick summary of the table
summary(res)

# check out the first few lines
head(res)

######################################################
# Get a table of shrunken log2 fold changes
######################################################

# get shrunken log fold changes
res_shrink <- lfcShrink(dds,coef="condition_treated_vs_control")

# plot the shrunken log2 fold changes against the raw changes:
plot(
	x=res$log2FoldChange,
	y=res_shrink$log2FoldChange,pch=20,
	cex=.2,
	col=1+(res$padj < 0.05),
	xlab="raw log2 fold change",
	ylab="shrunken log2 fold change"
	)
abline(0,1)

# get the top 20 genes by shrunken log2 fold change
res_shrink[order(-abs(res_shrink$log2FoldChange)),][1:20,]


######################################################
# Data visualization
######################################################

# MA plot
plotMA(res_shrink, ylim=c(-4,4))

##############

# PCA plot

# normalized, variance-stabilized transformed counts for visualization
vsd <- vst(dds, blind=FALSE)

plotPCA(vsd, intgroup="condition")

# alternatively, using ggplot
library(ggplot2)

dat <- plotPCA(vsd, intgroup="condition",returnData=TRUE)

p <- ggplot(dat,aes(x=PC1,y=PC2,col=group))
p + geom_point()

##############

# heatmap of DE genes

# load pheatmap library (install.packages(pheatmap) if you don't have it.)
library(pheatmap)

# regularized log transformation of counts
rld <- rlog(dds, blind=FALSE)

# get top 50 log fold change genes
top50 <- order(-abs(res_shrink$log2FoldChange))[1:50]
df <- data.frame(colData(dds)[,"condition"])
	rownames(df) <- colnames(dds)
	colnames(df) <- "condition"
pheatmap(
	assay(rld)[top50,], 
	cluster_rows=TRUE, 
	show_rownames=TRUE,
	cluster_cols=FALSE,
	annotation_col=df
	)

##############

# plot counts for individual genes

plotCounts(dds, gene=order(-abs(res_shrink$log2FoldChange))[1], intgroup="condition")




######################################################
# Get gene annotations using biomaRt
######################################################


# load biomaRt
library(biomaRt) 


##############################
# Select a mart and dataset
##############################

# see a list of "marts" available at host "ensembl.org"
listMarts(host="ensembl.org")

# create an object for the Ensembl Genes v100 mart
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="ensembl.org")

# see a list of datasets within the mart
	# at the time of writing, there were 203
listDatasets(mart)

# figure out which dataset is the croaker
	# be careful using grep like this. verify the match is what you want
searchDatasets(mart,pattern="lcrocea")

# there's only one match, get the name
croakerdata <- searchDatasets(mart,pattern="lcrocea")[,1]


# create an object for the croaker dataset
croaker_mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "ensembl.org", dataset = croakerdata)

#########################
# Query the mart/dataset
#########################

# filters, attributes and values

# see a list of all "filters" available for the lcrocea dataset.
	# at the time of writing, over 300
listFilters(croaker_mart)

# see a list of all "attributes" available
	# 129 available at the time of writing
listAttributes(mart = croaker_mart, page="feature_page")

# we can also search the attributes and filters
searchAttributes(mart = croaker_mart, pattern = "ensembl_gene_id")

searchFilters(mart = croaker_mart, pattern="ensembl")

# get gene names, when they exist
ann <- getBM(filter="ensembl_gene_id",value=rownames(res),attributes=c("ensembl_gene_id","description"),mart=croaker_mart)

# get GO term info
go_ann <- getBM(filter="ensembl_gene_id",value=rownames(res),attributes=c("ensembl_gene_id","description","go_id","name_1006","namespace_1003"),mart=croaker_mart)

######################################################
# Save some of the results to a table
######################################################

# set a prefix for output file names
outputPrefix <- "Croaker_DESeq2"

# save data results and normalized reads to csv
resdata <- merge(
  as.data.frame(res_shrink), 
  as.data.frame(counts(dds,normalized =TRUE)), 
  by = 'row.names', sort = FALSE
)
names(resdata)[1] <- 'gene'

write.csv(resdata, file = paste0(outputPrefix, "-results-with-normalized.csv"))

# send normalized counts to tab delimited file for GSEA, etc.
write.table(
  as.data.frame(counts(dds),normalized=T), 
  file = paste0(outputPrefix, "_normalized_counts.txt"), 
  sep = '\t'
)

