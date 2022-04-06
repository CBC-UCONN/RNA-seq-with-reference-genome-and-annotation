# this R code has been condensed from the tutorial here: 
	# github.com/CBC-UCONN/RNA-seq-with-reference-genome-and-annotation


# Load the libraries we'll need in the following code:
library("DESeq2")
library("apeglm")
library("pheatmap")
library("tidyverse")
library("ggrepel")

######################################################
# Point R to count data, set an output file prefix 
######################################################

# create an object with the directory containing your counts:
	# !!edit this to point to your own count file directory!!
directory <- "../06_count/counts"

# ensure the count files are where you think they are
list.files(directory)

sampleFiles <- list.files(directory, pattern = ".*counts")


######################################################
# Read the count data into R along with treatment information
######################################################

# load in metadata table
meta <- read.csv("../01_raw_data/metadata.txt") %>%
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
summary(res2)
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
arrange(data.frame(res_shrink1), log2FoldChange) %>% head(., n=20)


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

# heatmap of DE genes

# regularized log transformation of counts
rld <- rlog(dds, blind=FALSE)

# order gene names by absolute value of shrunken log2 fold change (excluding cook's cutoff outliers)
lfcorder <- data.frame(res_shrink3) %>%
  filter(!is.na(padj)) %>% 
  arrange(-abs(log2FoldChange)) %>% 
  rownames() 

# create a metadata data frame to add to the heatmaps
df <- data.frame(colData(dds)[,c("population","dose")])
  rownames(df) <- colnames(dds)
  colnames(df) <- c("population","dose")

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

# occasionally ensembl will have connectivity issues. we can try an alternative function:
	# select a mirror: 'www', 'uswest', 'useast', 'asia'
	# mart <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", mirror = "useast")

# see a list of datasets within the mart
	# at the time of writing, there were 203
listDatasets(mart)

# figure out which dataset is the killifish
	# be careful using grep like this. verify the match is what you want
searchDatasets(mart,pattern="Mummichog")

# there's only one match, get the name
killidata <- searchDatasets(mart,pattern="Mummichog")[,1]

# create an object for the killifish dataset
killi_mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "ensembl.org", dataset = killidata)

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
	# 129 available at the time of writing
listAttributes(mart = killi_mart, page="feature_page")

# we can also search the attributes and filters
searchAttributes(mart = killi_mart, pattern = "ensembl_gene_id")

searchFilters(mart = killi_mart, pattern="ensembl")

# get gene names and transcript lengths when they exist
ann <- getBM(filter="ensembl_gene_id",value=rownames(res1),attributes=c("ensembl_gene_id","description","transcript_length"),mart=killi_mart)

# pick only the longest transcript for each gene ID
ann <- group_by(ann, ensembl_gene_id) %>% 
  summarize(.,description=unique(description),transcript_length=max(transcript_length)) %>%
  as.data.frame()

# get GO term info
  # each row is a single gene ID to GO ID mapping, so the table has many more rows than genes in the analysis
go_ann <- getBM(filter="ensembl_gene_id",value=rownames(res1),attributes=c("ensembl_gene_id","description","go_id","name_1006","namespace_1003"),mart=killi_mart)

# get KEGG info
# kegg_ann <- getBM(filter="ensembl_gene_id",value=rownames(res),attributes=c("ensembl_gene_id","description","kegg_enzyme"),mart=killi_mart)

# put results and annotation in the same table
res_ann1 <- cbind(res_shrink1,ann)
res_ann2 <- cbind(res_shrink2,ann)
res_ann3 <- cbind(res_shrink3,ann)
res_ann4 <- cbind(res_shrink4,ann)



######################################################
# Do a couple GO enrichment analyses
######################################################

# first use 'goseq', a bioconductor package
  # goseq can automatically pull annotations for some organisms, but not F heteroclitus. 
  # we need to put together our own input data:
    # a vector of 1/0 for each gene, indicating DE/not DE
    # a vector of transcript lengths (the method tries to account for this source of bias)
    # a table of gene ID to category IDs (in this case GO term IDs) 


# load library `goseq`
library(goseq)

# 0/1 vector for DE/not DE
de <- as.numeric(res3$padj[!is.na(res3$padj)] < 0.1)
names(de) <- rownames(res3)[!is.na(res3$padj)]
# length of each gene
len <- ann[[3]][!is.na(res3$padj)]

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
res_ann3[gids,] %>% data.frame()


# plot log2 fold changes for those genes, sorted
ord <- order(res_ann3[gids,]$log2FoldChange)
plot(res_ann3[gids,]$log2FoldChange[ord],
     ylab="l2fc of genes in top enriched GO term",
     col=(res_ann3[gids,]$padj[ord] < 0.1) + 1,
     pch=20,cex=.5)
abline(h=0,lwd=2,lty=2,col="gray")


########################

# we can also run a similar analysis using the web platform gProfiler

# execute the following code, then copy the all ensembl gene IDs printed to the screen to your clipboard
cat(rownames(res3)[which(res3$padj < 0.1)])

# then visit biit.cs.ut.ee/gprofiler/gost
# select Laramichtys crocea as the organism
# paste the output into the query window and press "Run Query"
# if you explore the results, you'll see they are very similar. 

######################################################
# Save some of the results to a table
######################################################

# set a prefix for output file names
outputPrefix <- "killifish_DESeq2"

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

