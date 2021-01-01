my_libs = .libPaths()
my_libs
my_libs = c("/home/sashkoah/soft/rpackages", my_libs)
.libPaths(my_libs)
# BiocManager::install("Rtsne")

library(Rtsne)
library(affy)
library(ArrayExpress)
# library(arrayQualityMetrics)
library(sva)
library(limma)
library(stringr)
library(RColorBrewer)
library(factoextra)
library(stats)
library(usedist)
library(easypackages)
# libraries(getDependencies("STRINGdb"))



mergedpath = "/home/sashkoah/a/r/igea-r/merged/article_3/merged_dec_2020/23_no_55439"
# mergedpath = "/home/sashkoah/a/r/igea-r/merged/article_3/merged_earlier"
setwd('/home/sashkoah/a/r/igea-r')
source(paste(getwd(),'scripts/plots.R',sep='/'))
plotsqcpath = paste(getwd(), 'plots/qc/illumina', sep='/')

# preprocess affymetrix
# cbind
# 


# read merged expression data and metadata created with unite.R
mrgd = NULL
pdata = NULL
mrgd = read.table(file.path("temp", "mrgd.tsv"), sep="\t", header=TRUE)
pdata = read.table(file.path("temp","pdata.tsv"), sep="\t", header=TRUE)
ncol(mrgd)
nrow(pdata)



# # align metadata in pdata and expression data in mrgd so that rows in pdata match cols in mrgd
pdata$Array.Data.File = pdata$arraydatafile_exprscolumnnames

propercolnames = as.character(make.names(pdata$Array.Data.File))
propercolnames
rownames(pdata) = propercolnames
propercolnames
colnames(mrgd) = as.character(colnames(mrgd))
colnames(mrgd)
propercolnames
mrgd = mrgd[,propercolnames]
# check correctness
colnames(mrgd) == propercolnames
colnames(mrgd) == rownames(pdata)


pdata$Gestational.Age.Category
plotspath = "/home/sashkoah/a/r/igea-r/methods_article/plots"
datapath = "/home/sashkoah/a/r/igea-r/methods_article/data"
write.csv(mrgd, file.path(datapath, "mrgd_1_2.csv"), quote = FALSE)
write.csv(pdata, file.path(datapath, "pdata_1_2.csv"))

mrgd = read.csv(file.path(datapath, "mrgd_1_2.csv"), row.names = 1)
pdata = read.csv(file.path(datapath, "pdata_1_2.csv"),row.names = 1)


# remove technical batch effect caused by different datasets 
batch = as.factor(pdata$secondaryaccession)

mod = model.matrix(~ as.factor(trim_term)+as.factor(Combined.Fetus.Sex), data=pdata)
colnames(mrgd) == rownames(pdata)

exprs = mrgd
exprs = ComBat(dat=as.matrix(mrgd), batch=batch, mod=mod, par.prior=TRUE, prior.plots=FALSE)

dim(exprs)
sub_pdata = pdata
nrow(pdata)
ncol(exprs)
sub_exprs = exprs

table(sub_pdata[,c("secondaryaccession","Combined.Fetus.Sex")])

pca = prcomp(t(na.omit(sub_exprs)))
pca = prcomp(t(na.omit(mrgd)))

svg_width = 5
svg_height = 5
svg(file.path(plotspath, "pca_accession_exprs.svg"), width=svg_width, height=svg_height)
# png(file.path(plotspath, "pca_accession_mrgd1.png"), width=1000, height=1000)
pl = fviz_pca_ind(pca, axes = c(3,4),
                  label=c("var"),
                  habillage=sub_pdata$Combined.Fetus.Sex,
                  repel = TRUE,
                  palette = "Set4",
                  # col.ind = new_ga,
                  # gradient.cols = rainbow(4),
                  # select.var = list(name=as.character(merged_sel$SYMBOL)),
                  # select.var = list(contrib=1),
                  addEllipses = TRUE
)
pl
dev.off()

write.table(sub_exprs,file.path(datapath, "sub_exprs.tsv"), sep="\t", quote=FALSE)
write.table(sub_pdata,file.path(datapath, "sub_pdata.tsv"), sep="\t", quote=FALSE)

# Differential Expression Analysis
###################################################################################################

library(data.table)
library(openxlsx)
library(limma)

sub_exprs = read.table(file.path(datapath, "sub_exprs.tsv"), sep="\t")
sub_pdata = read.table(file.path(datapath, "sub_pdata.tsv"), sep="\t")

dim(sub_exprs)
dim(sub_pdata)
table(sub_pdata[,c( "trim")])
as.factor(sub_pdata$trim_term)

fit_mod = model.matrix(~ 0 + as.factor(sub_pdata$trim_term), data=sub_pdata)
colnames(fit_mod)
colnames(fit_mod) = c('i','ii')
fit <- lmFit(sub_exprs, fit_mod)  # fit each probeset to model
contrast.matrix <- makeContrasts(ii-i, levels=fit_mod)
fit2 <- contrasts.fit(fit, contrast.matrix)
efit <- eBayes(fit2)        # empirical Bayes adjustment

t = topTable(efit, number = nrow(sub_exprs))
nrow(t)
d = t
d = d[d$adj.P.Val<0.05,]
d = d[which(abs(d$logFC)>1),]
nrow(d)
diff = d
nrow(diff)
rownames(diff)
diff

library(org.Hs.eg.db)
keys(org.Hs.eg.db)

sel = AnnotationDbi::select(org.Hs.eg.db, rownames(diff), c("SYMBOL","GENENAME"))
sel = sel[!is.na(sel$SYMBOL),]

sel2 = AnnotationDbi::select(org.Hs.eg.db, rownames(sub_exprs), c("SYMBOL","GENENAME"))
sel2 = sel2[!is.na(sel2$SYMBOL),]

length(unique(rownames(sub_exprs)))==length(rownames(sub_exprs))
length(sel2$SYMBOL) == length(unique(sel2$SYMBOL))

NA %in% sel2$SYMBOL
NA %in% rownames(sub_exprs)
sub_exprs = sub_exprs[sel2$ENTREZID,]

nrow(sub_exprs)== nrow(sel2)
!(FALSE %in% (rownames(sub_exprs) == sel2$ENTREZID))
nrow(diff)
rownames(sub_exprs) = sel2$SYMBOL

write.csv(sub_exprs, file.path(datapath, "sub_exprs_tmp.csv"))



merged_sel = merge(sel,diff, by.x = "ENTREZID", by.y = "row.names")
nrow(merged_sel)

length(levels(as.factor(sub_pdata$trim_term)))

# add Average column to diffexp table
for (i in 1:length(levels(as.factor(sub_pdata$trim_term))))
{
  trim = levels(as.factor(sub_pdata$trim_term))[i]
  trim
  sub_exprs[which(rownames(sub_exprs) %in% merged_sel$SYMBOL),]
  difexp_exprs = sub_exprs[which(rownames(sub_exprs) %in% merged_sel$SYMBOL),]
  
  setdiff(rownames(difexp_exprs), merged_sel$SYMBOL)
  merged_sel[is.na(merged_sel$SYMBOL),]
  nrow(difexp_exprs)
  nrow(merged_sel)
  !(FALSE %in% (merged_sel$SYMBOL == rownames(difexp_exprs)))
  trim_sample_names = sub_pdata[sub_pdata$trim_term==trim,]$arraydatafile_exprscolumnnames
  
  difexp_exprs_trim = difexp_exprs[,which(colnames(difexp_exprs) %in% trim_sample_names)]
  merged_sel[,paste(trim,"Average",sep = " ")] = rowMeans(difexp_exprs_trim)
  merged_sel
}

nrow(merged_sel)
length(intersect(as.character(merged_sel$SYMBOL), as.character(row.names(sub_exprs))))
length((merged_sel$SYMBOL))

write.table(merged_sel,file.path(datapath,"difexp.tsv"), sep= '\t', row.names = FALSE)
write.table(sub_exprs, file.path(datapath,"sub_exprs_mapped.tsv"), sep= '\t')
