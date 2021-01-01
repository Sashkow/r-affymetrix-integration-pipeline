library(affy) 
library(ArrayExpress)
library(arrayQualityMetrics)
library(sva)
# library(stringr)
# library(ggplot2)
# library(ggfortify)
library(limma)
library(stats)
library(mvnfast)
# install.packages('mvnfast')


mergedpath =""
source(paste(getwd(),'scripts/plots.R',sep='/'))

setwd('/home/sashkoah/a/r/igea-r')
getwd()

# source('scripts/addon.R')


microarray_platform = 'illumina'

# rawspath = 'raws/affymetrix'
# prepath = 'preprocessed/affymetrix'
# pdatapath = 'pdata/'
# plotsqcpath = paste(getwd(), 'plots/qc', sep='/')
# mappedpath = 'mapped/affymetrix'
# mergedpath = 'merged'


# rawspath = 'raws/illumina'
# prepath = '/home/sashkoah/a/r/igea-r/preprocessed/article_3/placenta'
# pdatapath = 'pdata/illumina'
plotsqcpath = paste(getwd(), 'plots/illumina/qc', sep='/')
# mappedpath = '/home/sashkoah/a/r/igea-r/mapped/article_3/third article placenta 16676  92'
# mappedpath = "/home/sashkoah/a/r/igea-r/mapped/article_3/third article placenta december 2020"
# mappedpath = "/home/sashkoah/a/r/igea-r/mapped/article_3/third article placenta december 2020/placenta 1 2 3 trim no 55439"

mappedpath = '/home/sashkoah/a/r/igea-r/mapped/article_3/third article placenta december 2020/placenta 1 2 trim no 55439/'
# mappedpath = '/home/sashkoah/a/r/igea-r/mapped/article_3/placenta 2 3 trim'

# mappedpath = '/home/sashkoah/a/r/igea-r/mapped/lumiaffy'
# mappedpath = '/home/sashkoah/a/r/igea-r/mapped/article_4/placenta'

# mappedpath = '/home/sashkoah/a/r/igea-r/mapped/lumiaffy'

# mergedpath = 'merged'


# Load studies description
# studies <- read.table("general/affymetrix_placenta_studies.tsv", header = TRUE, sep = "\t")
studies <- read.table("general/illumina_placenta_studies.tsv", header = TRUE, sep = ",")

all_studies = read.table("general/accession_platform", header = TRUE, sep = '\t')
all_studies

# load sample metadata from csv.file from IGEA database which is based on ArrayExpress data mostly
# igea = read.table('igea_tsv/samples.csv',header = TRUE, sep = ',', fill = TRUE)
igea = read.table('igea_tsv/smoking_preeclampsia.csv',header = TRUE, sep = ',', fill = TRUE)

unique(igea[igea$Gestational.Age.Category=="First Trimester",]$accession)

decidua_exps = unique(igea[igea$Biological.Specimen=="Decidua",]$accession)

# get list of gene expression matrix files
exprs_files = list.files(mappedpath)
exprs_files
# f  = function(str){
#   return(strsplit(str,"_")[[1]][1])
# }
# exprs_files_short = lapply(exprs_files,f)
# exprs_files_short

mrgd = read.table(paste(mappedpath, exprs_files[1], sep = '/'), header = TRUE, sep = '\t')


# read first expression matrix
current_exprs = read.table(paste(mappedpath, exprs_files[1], sep = '/'), header = TRUE, sep = '\t')
current_exprs2 = read.table(paste(mappedpath, exprs_files[2], sep = '/'), header = TRUE, sep = '\t')

ncol(current_exprs)

# rownames(current_exprs)
# rownames(current_exprs2)
# gene_lists = list()
# 
# for (i in 1:length(exprs_files)){
#   current_exprs = read.table(paste(mappedpath, exprs_files[i], sep = '/'), header = TRUE, sep = '\t')
#   gene_lists[[i]] = rownames(current_exprs)
#   
# }
# names(gene_lists) = exprs_files_short
# 
# library(RAM)
# 
# # gene_lists = list(a=c("a","b"), b=c("b","c"))
# 
# 
# group.venn(gene_lists[1:5], label=FALSE,
#            # fill = c("orange", "blue"),
#            cat.pos = c(0, 0,1000,50,0),
#            # lab.cex=0,
#            cex=2)



# read each expression matrix and concatenate its columns, leaving only common rows
excluded_genes = c()

mrgd = current_exprs

for (exprs_file in exprs_files[2:length(exprs_files)]){
  current_exprs = read.table(paste(mappedpath, exprs_file, sep = '/'), header = TRUE, sep = '\t')
  # print(c(exprs_file,nrow(current_exprs)))
  mrgd = merge(mrgd, current_exprs, by = "row.names")
  rownames(mrgd) = mrgd$Row.names
  mrgd = mrgd[,!(colnames(mrgd) == "Row.names")]
  print(paste(exprs_file, nrow(mrgd), sep = ' '))
}

# excluded_genes = c()
# for (exprs_file in exprs_files[1:length(exprs_files)]){
#   current_exprs = read.table(paste(mappedpath, exprs_file, sep = '/'), header = TRUE, sep = '\t')
#   current_excluded_genes = setdiff(rownames(current_exprs),rownames(mrgd))
#   excluded_genes = union(as.character(excluded_genes),as.character(current_excluded_genes))
#   print(paste(exprs_file, length(excluded_genes), sep = ' '))
# }
# 
# all_possible_genes = c()
# for (exprs_file in exprs_files[1:length(exprs_files)]){
#   current_exprs = read.table(paste(mappedpath, exprs_file, sep = '/'), header = TRUE, sep = '\t')
#   all_possible_genes = union(as.character(rownames(current_exprs)),as.character(all_possible_genes))
#   print(paste(exprs_file, length(rownames(current_exprs)), length(all_possible_genes), sep = ' '))
# }
# 36-19
# 36481 - 19805
# length(excluded_genes)==length(unique(excluded_genes))
# 
# 
# library(org.Hs.eg.db)
# 
# keys(org.Hs.eg.db)
# sel = AnnotationDbi::select(org.Hs.eg.db, all_possible_genes, c("SYMBOL","GENENAME"))
# sel = AnnotationDbi::select(org.Hs.eg.db, excluded_genes, c("SYMBOL","GENENAME"))
# sel = AnnotationDbi::select(org.Hs.eg.db, rownames(mrgd), c("SYMBOL","GENENAME"))
# 
# 
# sel = sel[!is.na(sel$SYMBOL),]
# nrow(sel)

# map_exprs_entrez_to_symbol = function(exprs_file){
#   current_exprs = read.table(paste(mappedpath, exprs_file, sep = '/'), header = TRUE, sep = '\t')
#   sel = AnnotationDbi::select(org.Hs.eg.db, rownames(current_exprs), c("SYMBOL"))
#   sel = sel[!is.na(sel$SYMBOL),]
#   # print(c(exprs_file,length(rownames(current_exprs)),nrow(sel)))
#   current_exprs_trimmed = current_exprs[sel$ENTREZID,]
#   rownames(current_exprs_trimmed) = sel$SYMBOL
#   return(current_exprs_trimmed)
# }
# 
# current_exprs = read.table(paste(mappedpath, exprs_files[1], sep = '/'), header = TRUE, sep = '\t')
# current_exprs = map_exprs_entrez_to_symbol(current_exprs)
# dim(current_exprs)
# mrgd = current_exprs
# for (exprs_file in exprs_files[2:length(exprs_files)]){
#   # exprs_file = exprs_files[2]
#   current_exprs = map_exprs_entrez_to_symbol(exprs_file)
#   # print(c(exprs_file,nrow(current_exprs)))
#   mrgd = merge(mrgd, current_exprs, by = "row.names")
#   rownames(mrgd) = mrgd$Row.names
#   mrgd = mrgd[,!(colnames(mrgd) == "Row.names")]
#   print(paste(exprs_file, nrow(mrgd), sep = ' '))
#   
# }




ncol(mrgd)
nrow(mrgd)
# for affymetrix only
# pdata = igea[make.names(igea$Array.Data.File) %in% colnames(mrgd),]

# arrayexpress does not store processed exprs sample name column in a single place
# so this column was created manually for illumina
make.names(igea$exprs_column_names)

# arrayexpress stores affymetrix exprs raw file name in columns in Array.Data.File
make.names(igea$Array.Data.File)



# merge exprs samples file/colum names for illumina and affymetrix into a single column
new_column = character()
for (i in 1:nrow(igea)) {
  if (igea$Array.Data.File[i] != "_" & igea$Array.Data.File[i] != ""){
    value = as.character(igea$Array.Data.File[i])
    print(1)
  } else if (igea$exprs_column_names[i] != "") {
    value = as.character(igea$exprs_column_names[i])
    print(2)
  } else {
    print(paste(as.character(igea$secondaryaccession[i]), igea$exprs_column_names[i]))
    print(3)
    value = ""
  }
  new_column = c(new_column, value)
}

new_column
igea$arraydatafile_exprscolumnnames = new_column


getwd()
write(excluded_genes, "tmp.txt")




length(intersect(make.names(igea$arraydatafile_exprscolumnnames), colnames(mrgd)))
# length(intersect(make.names(igea$arraydatafile_exprscolumnnames),colnames(mrgd))) == length(colnames(mrgd))
# colnames(mrgd)[which(!(colnames(mrgd) %in% make.names(igea$arraydatafile_exprscolumnnames)))]


# get sample metadata only for expression data in mrgd

pdata = igea[make.names(igea$arraydatafile_exprscolumnnames) %in% colnames(mrgd),]
nrow(pdata)


# pdata = read.table("/home/sashkoah/a/r/igea-r/mapped/smoking/pdata.csv", sep="\t", header=TRUE)
# pdata$arraydatafile_exprscolumnnames = pdata$Array.Data.File
# pdata$secondaryaccession

# read-write metadata to elliminate unneeded factor levels
write.table(pdata,file.path("temp","pdata.tsv"), sep="\t", quote=FALSE)
pdata = read.table(file.path("temp","pdata.tsv"), sep="\t", header=TRUE)

# check how many samples we have in expression data and metadata
# must be equal
setdiff(colnames(mrgd),make.names(pdata$arraydatafile_exprscolumnnames))

ncol(mrgd)
nrow(pdata)

# pdata[pdata$trim=="Third Trimester",]$Gestational.Age.Appr



pdata$Diagnosis

# pdata = pdata[which(pdata$Diagnosis == "Healthy" &
#                     pdata$Biological.Specimen != "Maternal Blood" &
#                     pdata$Biological.Specimen != "Umbilical Cord Blood"),]
# 
# pdata$Smoking.Status
pdata$Biological.Specimen
# pdata = pdata[which(pdata$Biological.Specimen=="Placenta" & pdata$Diagnosis=="Healthy" & (pdata$Gestational.Age.Category=="Second Trimester" | pdata$Gestational.Age.Category=="Early Preterm" | pdata$Gestational.Age.Category=="Late Preterm" | pdata$Gestational.Age.Category=="Term")),]
pdata$Gestational.Age.Category
# pdata = pdata[which(pdata$Biological.Specimen=="Placenta" & pdata$Diagnosis=="Healthy" ),]
# pdata = pdata[which( pdata$Biological.Specimen=="Placenta" & pdata$Diagnosis=="Healthy" & (pdata$Gestational.Age.Category=="Second Trimester" | pdata$Gestational.Age.Category=="First Trimester")),]
pdata$Gestational.Age.Category

# pdata = pdata[which( pdata$Biological.Specimen=="Placenta" &
#                        pdata$Diagnosis=="Healthy" &
#                        (pdata$Gestational.Age.Category=="Term" | pdata$Gestational.Age.Category=="Second Trimester" )),]


# pdata = pdata[which( pdata$Biological.Specimen=="Placenta" &
#                      pdata$Diagnosis=="Healthy" &
#                      (pdata$Gestational.Age.Category=="Second Trimester" | pdata$Gestational.Age.Category=="First Trimester" | pdata$Gestational.Age.Category=="Term" )),]

table(pdata$Gestational.Age.Category)
      
# pdata = pdata[which(
#   pdata$Biological.Specimen=="Placenta" &
#   (pdata$Smoking.Status=="smoker" | pdata$Smoking.Status=="non-smoker" | pdata$Diagnosis=="Pre-Eclampsia" | pdata$Diagnosis=="Healthy") &
#   (pdata$Gestational.Age.Category %in% c("Early Preterm", "Late Preterm", "Term"))
# ),]

# for (accession in as.character(unique(pdata$secondaryaccession))) {
#   print(table(pdata[pdata$secondaryaccession==accession,]$Biological.Specimen))
# }

# pdata$Biological.Specimen

# pdata = pdata[which(
#   (pdata$Biological.Specimen=="Decidua" | pdata$Biological.Specimen=="Chorion" |  pdata$Biological.Specimen=="Placenta") & 
#     (pdata$Diagnosis=="Healthy" )
#     # (pdata$Gestational.Age.Category %in% c("Early Preterm", "Late Preterm", "Term"))
# ),]

# third article
pdata = pdata[which(
  (pdata$Biological.Specimen=="Placenta") & (pdata$Gestational.Age.Category=="First Trimester" |pdata$Gestational.Age.Category=="Second Trimester") & (pdata$Diagnosis=="Healthy")
),]


unique(pdata$secondaryaccession)
table(pdata$Diagnosis)

table(pdata[which(pdata$Diagnosis=="Pre-Eclampsia" | pdata$Diagnosis=="Severe Pre-Eclampsia"),]$Gestational.Age.Category)




# pdata = pdata[which(pdata$tissue.ch1=="placenta" & pdata$secondaryaccession=="GSE18044"),]

# outliers = c("X11761", "X11420", "X11670", "X13497", "X13521", "X14258")
# pdata = pdata[which((pdata$arraydatafile_exprscolumnnames %in% outliers)| pdata$accession=="E-GEOD-73685"),]

# pdata = pdata[which(pdata$accession=="E-GEOD-73685" &
#                     pdata$Diagnosis=="Healthy" &
#                       (pdata$Biological.Specimen == "Decidua" |
#                        pdata$Biological.Specimen == "Maternal Blood" |
#                        pdata$Biological.Specimen == "Umbilical Cord Blood")
# ),]
# pdata = pdata[which(pdata$Diagnosis=="Healthy" & 
#                     pdata$Gestational.Age.Category!="First Trimester" & pdata$Gestational.Age.Category!="Second Trimester" &
#                     pdata$Biological.Specimen!="Umbilical Cord Blood" &
#                     pdata$Biological.Specimen!="Maternal Blood" &
#                     pdata$Biological.Specimen!="Amnion" &
#                     pdata$Biological.Specimen!="Lower Segment" &
#                     pdata$Biological.Specimen!="Uterus Fundus")
# ,]

# allign expression data in mrgd with metadata in pdata
mrgd = mrgd[,make.names(pdata$arraydatafile_exprscolumnnames)]

unique(pdata$secondaryaccession)
# check colnames in mrgd match col arraydatafile_exprscolumnnames in pdata
setdiff(colnames(mrgd), make.names(pdata$arraydatafile_exprscolumnnames))

write.table(pdata,file.path("temp","pdata.tsv"), sep="\t", quote=FALSE)
# write.table(here,file.path("temp","pdata.tsv"), sep="\t", quote=FALSE)
write.table(mrgd, file.path("temp", "mrgd.tsv"), sep="\t", quote=FALSE)
nrow(pdata)
ncol(mrgd)
unique(pdata$secondaryaccession)
nrow(mrgd)
dim(mrgd)


table(pdata[,c("secondaryaccession","Fetus.Sex")])

for (acc in unique(pdata$secondaryaccession)){
  print(acc)
  print(table(pdata[pdata$secondaryaccession==acc,]$Smoking.Status))
  
}

mrgd = read.table(file.path("temp", "mrgd.tsv"), sep="\t", header=TRUE)
pdata = read.table(file.path("temp","pdata.tsv"), sep="\t", header=TRUE)
table(pdata[,c("secondaryaccession","Diagnosis")])

table(pdata$microarrays)
# 
# 
# # create Label Classes for CIBERSORT algorigthm for tissue mixture deconvolution
# pdata = read.table(file.path("temp","pdata.tsv"), sep="\t", header=TRUE)
# pdata$arraydatafile_exprscolumnnames
# 
# 
# 
# 
# # decidua 1 1 1 1 0 0 0 0 0 
# # blood1  0 0 0 0 1 1 0 0 0
# # blood2  0 0 0 0 0 0 1 1 1
# labels_df = data.frame(row.names = levels(pdata$Biological.Specimen))
# 
# for (tissue in rownames(labels_df)){
#   for (sample in pdata$arraydatafile_exprscolumnnames) {
#     if (pdata[pdata$arraydatafile_exprscolumnnames==sample,]$Biological.Specimen == tissue){
#       labels_df[tissue,sample] = 1
#     } else {
#       labels_df[tissue,sample] = 2
#     }
#   }
# }
# 
# labels_df
# colnames(labels_df) == colnames(mrgd)
# 
# decidua_samples_pdata = pdata[pdata$Biological.Specimen=="Decidua",]$arraydatafile_exprscolumnnames 
# decidua_samples_labels_df = colnames(labels_df["Decidua",which(labels_df["Decidua",] == 1)])
# decidua_samples_labels_df
# decidua_samples_pdata == decidua_samples_labels_df
#
# write.table(labels_df,file.path("temp","labels.tsv"), sep="\t", quote=FALSE)


 
# write.table(mrgd,file.path("temp","mrgd.tsv"), sep="\t", quote=FALSE)
# 
# setdiff(colnames(mrgd), make.names(igea$arraydatafile_exprscolumnnames))
# 
# 
# # end, next ask.R

