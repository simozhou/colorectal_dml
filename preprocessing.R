options(stringsAsFactors = FALSE)

# DATA IMPORT 

expr.matrix <- read.table("Documents/master/lab_data_mining/data/colorectal_FullExprMatrix_FPKM.tsv",
                          sep='\t', header=TRUE, row.names = 1, check.names = F)


## -------------------
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

