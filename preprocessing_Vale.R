options(stringsAsFactors = FALSE)

# DATA IMPORT 
setwd("~/Desktop/Data Mining")
expr.matrix <- read.table("colorectal_FullExprMatrix_FPKM.tsv",
                          sep='\t', header=TRUE, row.names = 1, check.names = F)


# -------------------
# We check the range of the data to see if a log transformation is required

range(expr.matrix)

# [1]      0.0 151505.6

# -------------------

# Log transformation of the data
# log(data + 1) to avoid -inf

expr.matrix <- log2(expr.matrix + 1)

# NORMALIZE TO Z-SCORES
# each column is normalized using z-score
# Z = ( x - mean(x) ) / std(x)

mean.cols <- apply(expr.matrix, 2, mean)
median.cols <- apply(expr.matrix, 2, median)

expr.matrix.mmedian <- sweep(expr.matrix, 2, median.cols, '-')

std.cols <- apply(expr.matrix, 2, sd)

zscore <- function(x) {
  
}
# mean centering + standard deviation scaling
expr.matrix.z <- scale(expr.matrix)

# median centered data -> standard deviation scaling
expr.matrix.m <- scale(expr.matrix.mmedian, center=F)

# plot of # of zeros vs IQR size for each gene (for both normalizations)

draw.iqr.dist <- function(dataset, threshold) {
  # IQR for every gene
  iqrs <- apply(dataset, 1, IQR)
  hist(iqrs,col = "mistyrose", main = "IQRs distribution", breaks = 80)
  abline(v=quantile(iqrs, threshold), lty=2, lwd=2, col="red")
}

draw.iqr.zeros <- function(dataset, threshold.iqr, threshold.zeros, ...) {
  # % of zero plotted against IQR values
  iqrs <- apply(dataset, 1, IQR)
  zeros <- apply(dataset==0, 1, sum)* 100/dim(dataset)[2]
  plot(iqrs, 100 - zeros, xlab = "IQR", ylab = "100 - % of zeros", cex=1.5, ...)
  text(quantile(iqrs, threshold.iqr),10, label = paste(as.character(threshold.iqr*100),"%"," min IQR", sep = ''), pos=4, col="red", cex=1.6)
  text(6,100-threshold.zeros*100, label = paste(as.character((1-threshold.zeros)*100),"%", " max zeros", sep = ''), pos=1, col="red", cex=1.6)
  abline(v=quantile(iqrs, threshold.iqr), h=100-threshold.zeros*100, lty=2, col="red")
}

iqr.zeros <- draw.iqr.zeros(expr.matrix, threshold.iqr = 0.75, threshold.zeros = 0.75, ylim=c(0,100), main="IQR size and percentage of zeros")

length(unique(colnames(expr.matrix)))


# Remove replicates (barcode with same vial)
library(Biobase)
no_vials <- isUnique(substring(colnames(expr.matrix), 1, 12))
table(no_vials)
expr.matrix <- expr.matrix[, no_vials]
dim(expr.matrix)

# Remove genes that have 0 expression across all samples 
table(rowSums(expr.matrix) == 0)
to_keep <- rowSums(expr.matrix) == 0
expr.matrix.nozero <- expr.matrix[!to_keep, ]
dim(expr.matrix.nozero)

# Load metadata 
metadata <- read.table("colorectal_clinical.tsv",
                       sep='\t', header=TRUE, row.names = 1, check.names = F)
# Replace row names with submitter_id column to match the sample names in the expression matrix
library(tidyverse)
metadata <- metadata %>% remove_rownames %>% column_to_rownames(var = "submitter_id")

# Change the sample names to remove the information about the vial (not available in the metadata)
new_names <- substring(colnames(expr.matrix), 1, 12)
colnames(expr.matrix) <- new_names
# Check how many of the samples in the expression matrix are contained in the metadata 
table(colnames(expr.matrix) %in% rownames(metadata)) # ALL 
# Check the opposite 
table(rownames(metadata) %in% colnames(expr.matrix)) # 9 samples in metadata are not in expression data (probably the replicates)

# Transpose the matrix for subsequent dimensionality reduction 
t.expr.matrix.nozero <- t(expr.matrix.nozero)
dim(t.expr.matrix.nozero)