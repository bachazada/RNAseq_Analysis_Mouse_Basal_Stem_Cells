#from https://combine-australia.github.io/RNAseq-R/06-rnaseq-day1.html 
#This study examines the expression profiles of basal stem-cell enriched cells (B) and committed luminal cells (L) in the mammary gland of virgin, pregnant and lactating mice. Six groups are present, with one for each combination of cell type and mouse status. Each group contains two biological replicates. We will first use the counts file as a starting point for our analysis. This data has already been aligned to the mouse genome. The command line tool featureCounts (Liao, Smyth, and Shi 2014) was used to count reads mapped to mouse genes from Refseq annotation (see the paper for details).


################################################
################################################
################################################
rm(list=ls())
# install edgeR & limma
y <- readRDS(file = "./data/data_mat.rds")
sampleinfo <- readRDS(file = "./data/sample_info.rds")
logcounts <- cpm(y,log=TRUE)
y <- calcNormFactors(y)

str(y)
y$samples
#The normalization factors multiply to unity across all libraries. A normalization factor below one indicates that the library size will be scaled down, as there is more suppression (i.e., composition bias) in that library relative to the other libraries. This is also equivalent to scaling the counts upwards in that sample. Conversely, a factor above one scales up the library size and is equivalent to downscaling the counts.

#The last two samples have much smaller normalisation factors, and MCL1.LA and MCL1.LB have the largest. If we plot mean difference plots using the plotMD function for these samples, we should be able to see the composition bias problem. We will use the logcounts, which have been normalised for library size, but not for composition bias.
which(colnames(y$counts)=="MCL1.LA")
which(colnames(y$counts)=="MCL1.LB")

pdf("M_versus_A_cpm.pdf", width=10, height=5)
par(mfrow=c(1,2))
plotMD(logcounts,column = 7)
abline(h=0,col="red")
plotMD(logcounts,column = 11)
abline(h=0,col="red")
dev.off()

pdf("M_versus_A_TMM.pdf", width=10, height=5)
par(mfrow=c(1,2))
plotMD(y,column = 7)
abline(h=0,col="red")
plotMD(y,column = 11)
abline(h=0,col="red")
dev.off()




# define the factors with names 
combinations <- paste(sampleinfo$CellType, sampleinfo$Status, sep=".")
# create factors from combinations
f <- factor(combinations, levels=unique(combinations))
design <- model.matrix(~0+f)
# beautify column names
colnames(design) <- levels(f)
design

pdf("voom.pdf",width=8,height=5)
par(mfrow=c(1,1))
v <- voom(y,design,plot = TRUE)
dev.off()


pdf(file = "log_cpm_vs_voom.pdf",width=10, height=5)
par(mfrow=c(1,2))
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2,main="Unnormalised logCPM")
## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
boxplot(v$E, xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")
## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(v$E),col="blue")
dev.off()

#Testing for differential expression
fit <- lmFit(v)
names(fit)
#Note that the group names must exactly match the column names of the design matrix.
cont.matrix <- makeContrasts(B.PregVsLac=basal.pregnant - basal.lactate,levels=design)
cont.matrix
#Now we can apply the contrasts matrix to the fit object to get the statistics and estimated parameters of our comparison that we are interested in. Here we call the contrasts.fit function in limma.



fit.cont <- contrasts.fit(fit, cont.matrix)

#The final step is to call the eBayes function, which performs empirical Bayes shrinkage on the variances, and estimates moderated t-statistics and the associated p-values.

fit.cont <- eBayes(fit.cont)
topTable(fit.cont)


summa.fit <- decideTests(fit.cont)
summary(summa.fit)


summa.fit <- decideTests(fit.cont, lfc = 1)
summary(summa.fit)

###Adding annotation and saving the results
library(org.Mm.eg.db)
columns(org.Mm.eg.db)
IDs = gsub("ID_","",rownames(fit.cont))
ann <- select(org.Mm.eg.db,keys=IDs,columns=c("ENTREZID","SYMBOL","GENENAME"))
head(ann)

table(paste("ID_",ann$ENTREZID,sep="")==rownames(fit.cont))

fit.cont$genes <- ann
topTable(fit.cont,coef="B.PregVsLac",sort.by="p")
limma.res <- topTable(fit.cont,coef="B.PregVsLac",sort.by="p",n="Inf")

write.csv(limma.res,file="limma_result.csv",row.names=FALSE)

pdf("significance_plots",width=10, height=5)
par(mfrow=c(1,2))
plotMD(fit.cont,coef=1,status=summa.fit[,"B.PregVsLac"], values = c(-1, 1))
# For the volcano plot we have to specify how many of the top genes to highlight.
# We can also specify that we want to plot the gene symbol for the highlighted genes.
# let's highlight the top 100 most DE genes
volcanoplot(fit.cont,coef=1,highlight=100,names=fit.cont$genes$SYMBOL)
dev.off()

sink("session_info.txt")
sessionInfo()
sink()
