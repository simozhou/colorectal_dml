library(sva)

expr.matrix <- read.table("filtered_genes_FPKM.csv",
                          sep=',', header=TRUE, row.names = 1, check.names = F)

meta.data <- read.table("../data/metadata_batch.csv",
                        sep='\t', header=TRUE, check.names = F)

total.samples <- intersect(colnames(expr.matrix), row.names(meta.data))

meta.data["TCGA-AU-6004",] <-	c("TCGA-COAD",	"female",	1941,	"white	not hispanic or latino",	
                                "Alive",	"Not reported",	824,	"Cecum",	"MSI",	"Stage I",	"Stage I",	
                                "T2",	"T2",	"N0",	"N0",	"M0",	"M0",	"Radiation Therapy, NOS",
                                -11.498749233762132,	-26.350431957105005,	8.892427444458008,	-2.485907554626465,	"AU",	"AU",	"St. Joseph's Medical Center-(MD)",	"Colon adenocarcinoma",	"IGC")

meta.data <- meta.data[total.samples,]

mod <- model.matrix(as.formula(~`Source Site`), data=meta.data)

mod0 <- model.matrix(~1, data = meta.data)

expr.matrix <- expr.matrix[,total.samples]

combat_edata <- ComBat(as.matrix(expr.matrix), batch = meta.data$`Source Site`, mod=mod0)

write.csv("../data/filtered_genes_FPKM_batchfree.csv", x = combat_edata)
