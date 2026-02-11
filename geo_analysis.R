# Step 0: Set working directory (optional)
getwd()   # check current directory

# Step 1: Load GEO dataset
library(GEOquery)
gse <- getGEO("GSE70970", GSEMatrix = TRUE)
exprSet <- exprs(gse[[1]])   # expression matrix
exprSet <- as.matrix(exprSet) # ensure numeric
head(exprSet[,1:5])

# Step 2: Basic analysis
summary(rowMeans(exprSet))   # average expression per gene

# Visualization of distribution
hist(rowMeans(exprSet),
     main = "Gene Expression Distribution",
     xlab = "Average Expression",
     col = "lightblue")

# Step 3: Differential Expression
pheno <- pData(gse[[1]])
head(pheno)

library(limma)

# IMPORTANT: build design matrix using a valid phenotype column
# Adjust column name based on your dataset
group <- factor(pheno$characteristics_ch1.1)   # replace with correct column
design <- model.matrix(~ group)

# Fit model
fit <- lmFit(exprSet, design)
fit <- eBayes(fit)
results <- topTable(fit, coef=2, number=Inf)

# Step 4: Volcano Plot
with(results, plot(logFC, -log10(P.Value),
                   pch=20, main="Volcano Plot",
                   xlab="Log2 Fold Change", ylab="-Log10 P-value"))

# Highlight significant genes
with(subset(results, adj.P.Val < 0.05 & abs(logFC) > 1),
     points(logFC, -log10(P.Value), pch=20, col="red"))

# Step 5: Heatmap
#install.packages("pheatmap")   # install if not already
library(pheatmap)

# Select top 50 DE genes
topGenes <- rownames(results)[1:50]
exprSubset <- exprSet[topGenes, ]
exprSubset <- as.matrix(exprSubset)

# Align phenotype metadata with samples
pheno <- pheno[colnames(exprSubset), ]

# Heatmap
pheatmap(exprSubset,
         scale="row",                # normalize each gene
         show_rownames=FALSE,
         annotation_col=pheno,
         main="Top 50 DE Genes Heatmap")

#Fixing heatmap
# Make sure pheatmap is installed
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}
library(pheatmap)

# Select top 50 DE genes
topGenes <- rownames(results)[1:50]
exprSubset <- exprSet[topGenes, ]
exprSubset <- as.matrix(exprSubset)   # ensure numeric

# Align phenotype metadata with samples
rownames(pheno) <- colnames(exprSet)   # set rownames explicitly

# Create a simple annotation dataframe
# Pick one useful column (adjust based on your dataset)
annotation <- data.frame(Group = pheno$characteristics_ch1.1)
rownames(annotation) <- rownames(pheno)

# Heatmap
pheatmap(exprSubset,
         scale="row",                # normalize each gene
         show_rownames=FALSE,
         annotation_col=annotation,  # simplified metadata
         main="Top 50 DE Genes Heatmap")
#Saving plots
# Save volcano plot
png("volcano_plot.png", width=800, height=600)
with(results, plot(logFC, -log10(P.Value),
                   pch=20, main="Volcano Plot",
                   xlab="Log2 Fold Change", ylab="-Log10 P-value"))
with(subset(results, adj.P.Val < 0.05 & abs(logFC) > 1),
     points(logFC, -log10(P.Value), pch=20, col="red"))
dev.off()

# Save heatmap
pdf("heatmap.pdf")
pheatmap(exprSubset,
         scale="row",
         show_rownames=FALSE,
         annotation_col=annotation,
         main="Top 50 DE Genes Heatmap")
dev.off()
#Saving results table
write.csv(results, "top_genes.csv")
