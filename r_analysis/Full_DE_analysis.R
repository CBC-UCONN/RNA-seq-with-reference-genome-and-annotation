# this R code has been condensed from the tutorial here: 
	# https://github.com/CBC-UCONN/RNA-seq-with-reference-genome-and-annotation


# Load the libraries we'll need in the following code:
library("DESeq2")
library("apeglm")
library("pheatmap")
library("tidyverse")


######################################################
# Point R to count data, set an output file prefix 
######################################################

# create an object with the directory containing your counts:
	# !!edit this to point to your own count file directory!!
directory <- "../count"

# ensure the count files are where you think they are
list.files(directory)

sampleFiles <- list.files(directory, pattern = ".*counts")


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
	ylab="shrunken log2 fold change",
	xlim=c(-5,5),
	ylim=c(-5,5)
	)
abline(0,1)

# get the top 20 genes by shrunken log2 fold change
top20 <- order(-abs(res_shrink$log2FoldChange))[1:20]
res_shrink[top20,]


######################################################
# Data visualization
######################################################

# MA plot
plotMA(res_shrink, ylim=c(-4,4))

##############

#Volcano plot

# negative log-scaled adjusted p-values
log_padj <- -log(res_shrink$padj,10)
log_padj[log_padj > 100] <- 100

# plot
plot(x=res_shrink$log2FoldChange,
     y=log_padj,
     pch=20,
     cex=.2,
     col=(log_padj > 10)+1, # color padj < 0.1 red
     ylab="negative log-scaled adjusted p-value",
     xlab="shrunken log2 fold changes")

#############

# PCA plot

# normalized, variance-stabilized transformed counts for visualization
vsd <- vst(dds, blind=FALSE)

plotPCA(vsd, intgroup="condition")

# alternatively, using ggplot

dat <- plotPCA(vsd, intgroup="condition",returnData=TRUE)

p <- ggplot(dat,aes(x=PC1,y=PC2,col=group))
p <- p + geom_point()
p

##############

# heatmap of DE genes

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

l2fc_ord <- order(-abs(res_shrink$log2FoldChange))
plotCounts(dds, gene=l2fc_ord[1], intgroup="condition")




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

# get gene names and transcript lengths when they exist
ann <- getBM(filter="ensembl_gene_id",value=rownames(res),attributes=c("ensembl_gene_id","description","transcript_length"),mart=croaker_mart)

# pick only the longest transcript for each gene ID
ann <- group_by(ann, ensembl_gene_id) %>% 
  summarize(.,description=unique(description),transcript_length=max(transcript_length))

# get GO term info
  # each row is a single gene ID to GO ID mapping, so the table has many more rows than genes in the analysis
go_ann <- getBM(filter="ensembl_gene_id",value=rownames(res),attributes=c("ensembl_gene_id","description","go_id","name_1006","namespace_1003"),mart=croaker_mart)

# get KEGG info
# kegg_ann <- getBM(filter="ensembl_gene_id",value=rownames(res),attributes=c("ensembl_gene_id","description","kegg_enzyme"),mart=croaker_mart)

# put results and annotation in the same table
res_ann <- cbind(res_shrink,ann)



######################################################
# Do a couple GO enrichment analyses
######################################################

# first use 'goseq', a bioconductor package
  # goseq can automatically pull annotations for some organisms, but not L crocea. 
  # we need to put together our own input data:
    # a vector of 1/0 for each gene, indicating DE/not DE
    # a vector of transcript lengths (the method tries to account for this source of bias)
    # a table of gene ID to category IDs (in this case GO term IDs) 


# load library `goseq`
library(goseq)

# 0/1 vector for DE/not DE
de <- as.numeric(res$padj < 0.1)
names(de) <- rownames(res)
# length of each gene
len <- ann[[3]]

# first try to account for transcript length bias by calculating the
# probability of being DE based purely on gene length
pwf <- nullp(DEgenes=de,bias.data=len)

# use the Wallenius approximation to calculate enrichment p-values
GO.wall <- goseq(pwf=pwf,gene2cat=go_ann[,c(1,3)])

# do FDR correction on p-values using Benjamini-Hochberg, add to output object
GO.wall <- cbind(
  GO.wall,
  padj_overrepresented=p.adjust(GO.wall$over_represented_pvalue, method="BH"),
  padj_underrepresented=p.adjust(GO.wall$under_represented_pvalue, method="BH")
  )

# explore the results

head(GO.wall)

# identify ensembl gene IDs annotated with to 2nd from top enriched GO term
  
g <- go_ann$go_id==GO.wall[2,1]
gids <- go_ann[g,1]

# inspect DE results for those genes
res_ann[gids,]


# plot log2 fold changes for those genes, sorted
ord <- order(res_ann[gids,]$log2FoldChange)
plot(res_ann[gids,]$log2FoldChange[ord],
     ylab="l2fc of genes in top enriched GO term",
     col=(res_ann[gids,]$padj[ord] < 0.1) + 1,
     pch=20,cex=.5)
abline(h=0,lwd=2,lty=2,col="gray")


########################

# we can also run a similar analysis using the web platform gProfiler

# execute the following code, then copy the all ensembl gene IDs printed to the screen to your clipboard
cat(rownames(res)[res$padj < 0.1])

# then visit https://biit.cs.ut.ee/gprofiler/gost
# select Laramichtys crocea as the organism
# paste the output into the query window and press "Run Query"
# if you explore the results, you'll see they are very similar. 

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

