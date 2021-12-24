#GBM PSets:
# work_dir = "C:/BHK/PSet/GBM/Pachyderm/"
# setwd(work_dir)

input_dir <- "/pfs/downloadGBMData/"
input_annotation <- "/pfs/downAnnotations/"
out_dir <- "/pfs/out/" 

# ==== Data ====
load(paste0(input_dir, "Ensembl.v99.annotation.RData"))
# GSE152160_RAW.tar -> "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE152160&format=file" 
# mmc2.xlsx -> "https://www.cell.com/cms/10.1016/j.celrep.2020.107897/attachment/eb4d5373-4278-40fd-8377-8bd36387bacf/mmc2.xlsx",rowNames = TRUE)
# HK_genes.txt -> "https://www.tau.ac.il/~elieis/HKG/HK_genes.txt"
# HGCC_DNA_copy_number_gene_level.txt -> "http://portal.hgcc.se/data/HGCC_DNA_copy_number_gene_level.txt"
# HGCC_WES_mutations_variants.txt -> "http://portal.hgcc.se/data/HGCC_WES_mutations_variants.txt"
# HGCC_DNA_methylation.txt -> "http://portal.hgcc.se/data/HGCC_DNA_methylation.txt"
# 1-s2.0-S221359601630071X-mmc1.txt -> "https://ars.els-cdn.com/content/image/1-s2.0-S221359601630071X-mmc1.txt"
# 1-s2.0-S221359601630071X-mmc2.txt -> "https://ars.els-cdn.com/content/image/1-s2.0-S221359601630071X-mmc2.txt"
# drug_with_ids.csv -> From pachyderm
# mmc3.xlsx -> "https://www.cell.com/cms/10.1016/j.celrep.2020.107897/attachment/9241fbaa-f27a-430f-82c9-3f38c4b7e062/mmc3.xlsx"
# HGCC_drug_response_AUC.txt -> "http://portal.hgcc.se/data/HGCC_drug_response_AUC.txt"

# ======================== Packages ========================
#Install brainarray package:
#download the library from here: http://mbni.org/customcdf/24.0.0/ensg.download/hta20hsensgcdf_24.0.0.tar.gz
#If required run ----> BiocManager::install("AnnotationDbi")

# Troubleshooting -- Error; return code from  pthread_create() is 22 when using affy::justRMA()
# https://support.bioconductor.org/p/122925/#124701 suggests that this issue comes from an affy dependency called preprocessCore library.
# The following removes affy and preprocessCore libraries, reinstalls preprocessCore with threading disabled.
remove.packages("preprocessCore")
BiocManager::install("preprocessCore", configure.args="--disable-threading")

# Loading packages
library(data.table)
library(openxlsx)
library(stringr)
library(Biobase)
library(GEOquery)
library(affy)
library(oligo)
library(affycoretools)
library(pd.hta.2.0)
library(RUVnormalize)
library(RUVnormalizeData)
library(PharmacoGx)
library(abind)
library(reshape2)
library(CoreGx)
library(AnnotationDbi)
install.packages(paste0(input_dir, "hta20hsensgcdf_24.0.0.tar.gz"), sep="", repos = NULL, type = "source")
#Functions
# ======================== drug_name_correction ========================
#Function to find the synonym drug names by not considering the salts and unwanted punctuation
drug_name_correction <-function(table, column){
  
  table$clean.ids <- tolower(str_replace_all(table[, column] , "[^[:alnum:]]", ""))
  
  salts <- c("malate","sulfate","dihydrochloride","hydrochloride","citrate","ethanolamine",
             "2hcl", "oxalate" , "sodiumsalt" ,"bromide" , "monohydrate" ,"isethionate","sodium",
             "furoate","hcl","ca","dimaleate", "oxalate","dihydrate","maleate","fumarate","lactate","di")
  for (salt in salts){ 
    table[, "clean.ids"] <- gsub(salt,"", table[, "clean.ids"], ignore.case=T)
  }
  return(table)
}

# ======================== fdata_builder ========================
#Function to annotate the genes when creating feature data for a PSet. If annotation for a gene is not available 
#it will not be removed but its annotation will be presented by "NA".

fdata_builder<-function(annotation, assay,ID_column="gene_name"){
  feat_empty<-data.frame(matrix(NA, nrow = nrow(assay), ncol = ncol(annotation)))
  colnames(feat_empty)<-colnames(annotation)
  feat_empty[,ID_column]<-rownames(assay)
  
  annotated<-annotation[annotation[,ID_column] %in% rownames(assay),] #Subseting the genes from features_gene that belong to cnv
  feat<-merge(x=feat_empty, y=annotated, all=TRUE) #We keep all the genes even the ones that don't have annotations.
  feat<-feat[!duplicated(feat[,ID_column]),]
  
  #Reformatting feature data to match the PSet requirements:
  feat<-feat[match(rownames(assay), feat[,ID_column]),]
  rownames(feat)<-feat[,ID_column]
  feat$Symbol<-feat$gene_name
  return(feat)}

# ======================== ph_data_builder ========================
#Function to create the pheno data from cell-line object (given it includes the meta data)

ph_data_builder<- function(annotation,assay){
  
  phen<-as.data.frame(colnames(assay))
  phen$temp<-sub("U","",phen[,1])
  phen$temp<-substring(phen$temp, 1, 4)
  annotation$temp<-substring(rownames(annotation), 2, 5)
  phen<-merge(phen, annotation, by="temp" , all.x=TRUE)
  
  phen$batchid <- NA
  phen$cellid<-phen$`colnames(assay)`
  phen$Replicate <- substring(phen$cellid,6)
  phen$Replicate[phen$Replicate ==""]<- NA
  rownames(phen)<-phen$`colnames(assay)`
  
  phen<-phen[,c(15,18,16,17, 9:13)]
  return(phen)}

# ======================== eSetToSE ========================
# A function converting ExpressionSet to SummarizedExperiment

eSetToSE <- function(eSet , annot_name) {
  
  BiocGenerics::annotation(eSet) <- annot_name
  stopifnot(all(rownames(fData(eSet)) == rownames(exprs(eSet))))
  stopifnot(all(rownames(pData(eSet)) == colnames(exprs(eSet))))
  
  # Build summarized experiment from eSet
  SE <- SummarizedExperiment::SummarizedExperiment(
    assays=SimpleList(as.list(Biobase::assayData(eSet))),
    # Switch rearrange columns so that IDs are first, probes second
    rowData=S4Vectors::DataFrame(Biobase::fData(eSet)),
    colData=S4Vectors::DataFrame(Biobase::pData(eSet)),
    metadata=list("experimentData" = eSet@experimentData, 
                  "annotation" = Biobase::annotation(eSet), 
                  "protocolData" = Biobase::protocolData(eSet)))
  # Extract names from expression set                  
  SummarizedExperiment::assayNames(SE) <- Biobase::assayDataElementNames(eSet)
  
  stopifnot(all(rownames(colData(SE)) == rownames(pData(eSet))))
  stopifnot(all(rownames(rowData(SE)) == rownames(fData(eSet))))
  return(SE)
}

# ======================== GSM_maping ========================
# Data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152160
# which includes GSM ids of the samples mapped to their corresponding cell ids

GSM_html<-readLines("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152160")

GSM_map<- data.frame( accession_id= GSM_html[grep(">GSM", GSM_html)] , cellid= GSM_html[grep(">GSM", GSM_html) + 1])
GSM_map $accession_id <- str_sub( GSM_map $accession_id, -19,-10)
GSM_map $cellid <- str_sub(sub('<td valign=\"top\">' ,'', GSM_map $cellid) , end = -6)

# Creating replicate column including replicate markers (e.g. "A") 
GSM_map$Replicate <- gsub(".*cells", "", GSM_map$cellid)
GSM_map$Replicate [GSM_map$Replicate =="" | GSM_map$Replicate =="human astrocytes"]<-NA
GSM_map$Replicate <- gsub(" ", "", GSM_map$Replicate)


# Creating Patient-id column consistant with cell names in the main paper
GSM_map$Patient_id <-gsub("cells.*", "", GSM_map$cellid)
GSM_map$Patient_id <- gsub(" ", "", GSM_map$Patient_id)
GSM_map$Patient_id [GSM_map$Patient_id =="humanastrocytes"]<-"human_astrocytes"

# Creating unique cell-ids based on patient ids and replicates
GSM_map$cellid<-sub("MG cells", "", GSM_map$cellid)
GSM_map$cellid<-sub(" ", "_", GSM_map$cellid)

# GSM ids are only used in the expression data
GSM_map<-GSM_map[, c("Patient_id","Replicate","cellid", "accession_id")]

# ======================== Molecular Profiles ========================
# ============= expression data =============
#Creating eset from raw expression data 
# GSE152160_RAW.tar file to be downloaded from here "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE152160&format=file" 
untar(paste(input_dir, "GSE152160_RAW.tar", sep=""), exdir="GSE152160_RAW")#Unpack the CEL files
cels<-list.files("GSE152160_RAW/", pattern = "CEL")
for(cel in cels){
  if(!file.exists(paste0("GSE152160_RAW/", gsub("[.]gz$", "", cel)))){
    GEOquery::gunzip(filename=paste0("GSE152160_RAW/", cel), overwrite = TRUE, remove = FALSE)
  }
}
# sapply(paste("GSE152160_RAW", cels, sep="/"), GEOquery::gunzip)

####Gene-level expression
cdf <- "hta20hsensgcdf"
cels <- list.celfiles("GSE152160_RAW/", full.names = TRUE)#Raw_expression Folder includes 145 CEL files
expr_cel <- justRMA(filenames = cels, verbose = TRUE, cdfname = cdf)

# Assay data expression 
assay_exp<-as.data.frame(exprs(expr_cel))
rownames(assay_exp)<-sub("_at","",rownames(assay_exp))
colnames(assay_exp)<-sub("_PA.*","",colnames(assay_exp))

# Removing control genes (which start with "AFFX")
assay_exp<-assay_exp[-c(grep("AFFX", rownames(assay_exp))),] 


# Feature data expression
feat_exp<-fdata_builder(annotation=features_gene, assay=assay_exp,ID_column="gene_id")#features_gene from Ensembel.v99.annotation.RData
feat_exp$BEST<-NA

# Pheno data expression
cell<-read.xlsx(paste(input_dir, "mmc2.xlsx", sep=""),rowNames = TRUE)#Cell_line names
cell$Patient_id<-rownames(cell)
phen_exp<-merge(GSM_map , cell, by="Patient_id" , all.x=TRUE)
phen_exp<-phen_exp[, c(1:4 , 11:15)]
phen_exp$batchid <- NA
rownames(phen_exp)<-phen_exp$accession_id

# Protocol data expression normalized based on RMA only
protocol_exp<-as.data.frame(rep("Affymetrix HTA 2.0 array",ncol(assay_exp)),row.names=colnames(assay_exp))
colnames(protocol_exp)<-"Array"
protocol_exp$URL_raw_data<-"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152160" #Link to access raw expression data
protocol_exp$Folder_raw_data<-"GSE152160_RAW.tar" #Folder including the raw data in the above link

#Creating ExpressionSet normalized based on RMA only
assay_exp<-assay_exp[,rownames(phen_exp)] #rearranging the colnames so it is similar to pheno data
protocol_exp<-protocol_exp[rownames(phen_exp),]#rearranging the rownames so it is similar to pheno data
exp_eSet<- ExpressionSet(assayData = as.matrix(assay_exp), phenoData = AnnotatedDataFrame(phen_exp), 
                         featureData = AnnotatedDataFrame(feat_exp),
                         protocolData=AnnotatedDataFrame(protocol_exp)) 

# ============= Normalizing based on negative control genes =============
ctrl<-read.delim(paste(input_dir, "HK_genes.txt", sep=""),stringsAsFactors = FALSE, header = FALSE) # Negative control genes 
ctrl$gene_name<-gsub(" ","",ctrl$V1) #Removing spaces from the gene names
ctrl<-merge(ctrl, feat_exp[,c("gene_name","gene_id")], by="gene_name")
ctrl_ind<-which(rownames(assay_exp) %in% ctrl$gene_id) 

## Prepare control samples
# A table that has as many columns as the largest set of replicates for one sample. Each
# row corresponds to a set of replicates of the same sample and gives the row indices of the
# replicates in the gene expression matrix, padded with -1 entries.
# See https://www.bioconductor.org/packages/release/bioc/vignettes/RUVnormalize/inst/doc/RUVnormalize.pdf

YY <- t(assay_exp) 
rep_cIdx <- matrix(-1,145,2) # 145 is number of cell lines and 2 is max number of replicates a cell has
added_pat<-c()
count <-1#row number of rep_cIdx

for (i in 1:nrow(rep_cIdx)){
  pat = phen_exp$Patient_id[i]
  if (pat %in% added_pat)
    next
  GSM=phen_exp$accession_id[phen_exp$Patient_id==pat]
  if(length(GSM)==2){
    rep_cIdx[count,1]=which(rownames(YY)==GSM[1])
    rep_cIdx[count,2]=which(rownames(YY)==GSM[2])
  }
  else(rep_cIdx[count,1]=which(rownames(YY)==GSM))
  
  added_pat=c(added_pat,pat)  
  count=count+1
}

## Replicate-based
# Set k to the number of samples / 4 or to the number of replicates, if the latter is smaller than the former. 
# See https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4896498/

Res <- naiveReplicateRUV(YY, ctrl_ind, rep_cIdx, k=17)

# Assay expression RUV
assay_exp_ruv<-data.frame(t(Res$cY))

# Protocol data expression RUV
protocol_exp_ruv<-as.data.frame(rep("Affymetrix HTA 2.0 array",ncol(assay_exp_ruv)),row.names=colnames(assay_exp_ruv))
colnames(protocol_exp_ruv)<-"Array"
protocol_exp_ruv$URL_raw_data<-"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152160" #Link to access raw expression data
protocol_exp_ruv$Folder_raw_data<-"GSE152160_RAW.tar" #Folder including the raw data in the above link
protocol_exp_ruv$Negative_control_genes<-"https://www.tau.ac.il/~elieis/HKG/"

# Creating ExpressionSet RUV
assay_exp_ruv<-assay_exp_ruv[,rownames(phen_exp)] #rearranging the colnames so it is similar to pheno data
protocol_exp_ruv<-protocol_exp_ruv[rownames(phen_exp),]#rearranging the rownames so it is similar to pheno data
exp_eSet_ruv<- ExpressionSet(assayData = as.matrix(assay_exp_ruv), phenoData = AnnotatedDataFrame(phen_exp), 
                             featureData = AnnotatedDataFrame(feat_exp),
                             protocolData=AnnotatedDataFrame(protocol_exp_ruv)) 

# ============= Probe-level expression =============
raw_data<- oligo::read.celfiles(cels) # If required run this: memory.limit(size=76000) , memory.limit()
norm_data<-rma(raw_data,target="core") # Perform RMA normalization
norm_main<- getMainProbes(norm_data, level = "core")#Remove the control transcripts (Controls do not match to any genes)
norm_main_annot<- annotateEset(norm_main, pd.hta.2.0)#Annotating the eset on transcript level : "https://support.bioconductor.org/p/89308/"

# Assay data expression-probe level
assay_exp_probe<-as.data.frame(exprs(norm_main_annot))
colnames(assay_exp_probe)<-sub("_PA.*","",colnames(assay_exp_probe))

# Feature data expression-probe level
feat_exp_probe<-fData(norm_main_annot)
feat_exp_probe$BEST<-NA

# Pheno data expression-probe level is same as Pheno data expression

# Protocol data expression-probe level
protocol_exp_probe<-as.data.frame(rep("Affymetrix HTA 2.0 array",ncol(assay_exp_probe)),row.names=colnames(assay_exp_probe))
colnames(protocol_exp_probe)<-"Array"
protocol_exp_probe$URL_raw_data<-"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152160" #Link to access raw expression data
protocol_exp_probe$annotation<-"pd.hta.2.0"

#Creating ExpressionSet 
assay_exp_probe<-assay_exp_probe[,rownames(phen_exp)]#rearranging the colnames so it is similar to pheno data
protocol_exp_probe<-protocol_exp_probe[rownames(phen_exp),]#rearranging the rownames so it is similar to pheno data
exp_eSet_probe<- ExpressionSet(assayData = as.matrix(assay_exp_probe), phenoData = AnnotatedDataFrame(phen_exp), 
                               featureData = AnnotatedDataFrame(feat_exp_probe),
                               protocolData=AnnotatedDataFrame(protocol_exp_probe)) 

# ============= CNV data =============

assay_cnv<- read.delim(paste(input_dir,"HGCC_DNA_copy_number_gene_level.txt", sep=""), header=T, sep="\t", stringsAsFactors = FALSE)
rownames(assay_cnv)<-assay_cnv$Row
assay_cnv<-assay_cnv[,-1]

feat_cnv<-fdata_builder(annotation=features_gene, assay=assay_cnv,ID_column="gene_name")
phen_cnv<-ph_data_builder(annotation=cell,assay=assay_cnv)
phen_cnv$Replicate[phen_cnv$cellid=="U10000" |phen_cnv$cellid=="U10001" ]<-NA

#Creating ExpressionSet 
assay_cnv<-assay_cnv[,rownames(phen_cnv)]#rearranging the colnames so it is similar to pheno data
cnv_eSet<- ExpressionSet(assayData = as.matrix(assay_cnv), phenoData = AnnotatedDataFrame(phen_cnv), featureData = AnnotatedDataFrame(feat_cnv)) 

# ============= Mutation data =============
#Assay data
mutation<- read.delim(paste(input_dir, "HGCC_WES_mutations_variants.txt", sep=""), header=T, sep="\t", stringsAsFactors = FALSE,na.strings=c("", " ", "NA","NaN"))

mutation$CELL_LINE<- paste("U", mutation $CELL_LINE, sep="")
assay_mut<-data.frame(matrix("wt", nrow = length(unique(mutation$Gene.refGene[!is.na(mutation$Gene.refGene)])), ncol = length(unique(mutation$CELL_LINE))))
rownames(assay_mut)<-unique(mutation$Gene.refGene[!is.na(mutation$Gene.refGene)])
colnames(assay_mut)<-unique(mutation$CELL_LINE)

#Labeling wild type cells as "wt"
for (i in 1:nrow(assay_mut)){
  if(i%%100 ==0 ){print(i)}
  gene<-rownames(assay_mut)[i]
  
  indS<-which(mutation$Gene.refGene == gene)
  for (ind in indS){
    cell <- mutation$CELL_LINE [ind]
    func <- mutation$ExonicFunc.refGene [ind]
    
    if (!is.na(func)){
      if(assay_mut[gene, cell]=="wt"){assay_mut[gene, cell]=func}
      else if(grepl(func, assay_mut[gene, cell], perl=TRUE)) {next} #If the mutation effect is already reported for a gene then we skip
      else{assay_mut[gene, cell] = paste(assay_mut[gene,cell],func, sep="///")
      }} }}

assay_mut<-assay_mut[,-2]

feat_mutation<-fdata_builder(annotation=features_gene, assay=assay_mut)
phen_mutation<-ph_data_builder(annotation=cell,assay=assay_mut)
phen_mutation$Patient_id[phen_mutation$cellid=="U3067"]<-"U3067MG"# Based on "GSM_map" data frame

#Creating ExpressionSet 
assay_mut<-assay_mut[,rownames(phen_mutation)]#rearranging the colnames so it is similar to pheno data
mutation_eSet<- ExpressionSet(assayData = as.matrix(assay_mut), phenoData = AnnotatedDataFrame(phen_mutation), featureData = AnnotatedDataFrame(feat_mutation)) 

# ============= Methylation data =============
#Assay data
assay_methyl<- read.delim(paste(input_dir, "HGCC_DNA_methylation.txt", sep=""), header=T, sep="\t",stringsAsFactors = FALSE)
illum<-read.csv(paste(input_dir,"MethylationEPIC_v-1-0_B2.csv", sep=""),skip=7, stringsAsFactors = FALSE) #The first 7 rows are not informative 

#Removing SNP probes and sex chromosome probes
assay_methyl<-assay_methyl[-c(grep("rs",rownames(assay_methyl))),]
assay_methyl<-assay_methyl[-c(grep("ch.X",rownames(assay_methyl))),]

#Removing cross-reactive and polymorphic probes
polymorph_probes<-read.delim(paste(input_dir, "1-s2.0-S221359601630071X-mmc1.txt", sep=""), stringsAsFactors = FALSE)
cross_probes<-read.delim(paste(input_dir,"1-s2.0-S221359601630071X-mmc2.txt", sep=""), header=F , stringsAsFactors = FALSE)
assay_methyl<-assay_methyl[rownames(assay_methyl) %in% polymorph_probes$IlmnID == FALSE,]
assay_methyl<-assay_methyl[rownames(assay_methyl) %in% cross_probes$V1 == FALSE,]
colnames(assay_methyl)[colnames(assay_methyl) == "hAstro"]<-"human_astrocytes"

#Assay data gene-level
assay_methyl_gene<-read.delim(paste(input_dir ,"methMat.txt", sep=""), sep="\t", stringsAsFactors = FALSE)
colnames(assay_methyl_gene)[colnames(assay_methyl_gene) == "hAstro"]<-"human_astrocytes"

#Feature data
feat_methyl<-fdata_builder(annotation=illum, assay=assay_methyl,ID_column="IlmnID")
feat_methyl_gene<-fdata_builder(annotation=features_gene, assay=assay_methyl_gene)

#Pheno data
phen_methyl<-ph_data_builder(annotation=cell, assay=assay_methyl)
phen_methyl$Patient_id[phen_methyl$cellid=="human_astrocytes"]<-"human_astrocytes"
phen_methyl$Replicate[phen_methyl$Replicate=="_astrocytes"]<-NA

#Pheno data methylation gene level is same as Pheno data probe level

#Protocol data
protocol_methyl<-as.data.frame(rep("Infinium® MethylationEPIC BeadChip Kit",ncol(assay_methyl)),row.names=colnames(assay_methyl))
colnames(protocol_methyl)<-"Array"
protocol_methyl$Provider<-"Illumina"
protocol_methyl$URL_GpG<-"http://portal.hgcc.se/" #Link to access methylation data provided by the authors
protocol_methyl$File_GpG<-"HGCC_DNA_methylation.txt" #File including the data in the above link
protocol_methyl$URL_annotation<-"https://support.illumina.com/array/array_kits/infinium-methylationepic-beadchip-kit/downloads.html"#Link to access annotations by illumina
protocol_methyl$File_annotation<-"Infinium MethylationEPIC v1.0 B5 Manifest File (CSV Format)"#File including annotations in the above link

protocol_methyl_gene<-as.data.frame(rep("Infinium® MethylationEPIC BeadChip Kit",ncol(assay_methyl_gene)),row.names=colnames(assay_methyl_gene))
colnames(protocol_methyl_gene)<-"Array"
protocol_methyl_gene$Provider<-"Illumina"

#Creating ExpressionSet 
assay_methyl<-assay_methyl[,rownames(phen_methyl)]#rearranging the colnames so it is similar to pheno data
protocol_methyl<-protocol_methyl[rownames(phen_methyl),]#rearranging the rownames so it is similar to pheno data
methyl_eSet<- ExpressionSet(assayData = as.matrix(assay_methyl), phenoData = AnnotatedDataFrame(phen_methyl), 
                            featureData = AnnotatedDataFrame(feat_methyl), protocolData=AnnotatedDataFrame(protocol_methyl)) 

assay_methyl_gene<-assay_methyl_gene[,rownames(phen_methyl)]#rearranging the colnames so it is similar to pheno data
protocol_methyl_gene<-protocol_methyl_gene[rownames(phen_methyl),]#rearranging the rownames so it is similar to pheno data
methyl_gene_eSet<- ExpressionSet(assayData = as.matrix(assay_methyl_gene), 
                                 phenoData = AnnotatedDataFrame(phen_methyl), 
                                 featureData = AnnotatedDataFrame(feat_methyl_gene),
                                 protocolData=AnnotatedDataFrame(protocol_methyl_gene)) 

# ============= SE objects =============
#Checks are included in the eSetToSE function
expression_SE<- eSetToSE(exp_eSet,annot_name="rna")
expression_ruv_SE<- eSetToSE(exp_eSet_ruv,annot_name="rna_ruv")
expression_probe_SE<- eSetToSE(exp_eSet_probe,annot_name="rna_probe")
cnv_SE <- eSetToSE(cnv_eSet,annot_name="cnv")
mutation_SE <- eSetToSE(mutation_eSet,annot_name="mut")
methyl_gene_SE <- eSetToSE(methyl_gene_eSet,annot_name="methyl_gene")
methyl_SE <- eSetToSE(methyl_eSet,annot_name="methyl_probe")

# ======================== Cell object ========================
#Gathering cell lines from all the experiments
phen_exp<-phen_exp[,colnames(phen_cnv)] # Reordeing the columns of phen_exp for rbind
all_cell_obj<-as.data.frame(unique(rbindlist( list(phen_exp,phen_cnv,phen_mutation,phen_methyl))))
rownames(all_cell_obj)<-all_cell_obj$cellid

# ======================== Drug annotation data from Pachy annotation ========================
drug_with_ids_gbm <- read.csv(paste(input_annotation, "drugs_with_ids.csv", sep = "") , stringsAsFactors = FALSE, na.strings = "") #this is the drug with ids file downloaded form pachy annotation repo

#To be corrected in "drug_with_ids.csv"
conc.name <- "Doxorubicin1///Doxorubicin2///Doxorubicin3///Doxorubicin4///Doxorubicin5///Doxorubicin6///Doxorubicin7///Doxorubicin8///Doxorubicin hydrochloride"
ind= which(drug_with_ids_gbm$unique.drugid == "Doxorubicin")
drug_with_ids_gbm$GBM.drugid[ind] <- conc.name # TO BE DELETED ONCE CORRECTED
drug_with_ids_gbm$GBM.drugid[which(drug_with_ids_gbm$unique.drugid == "Carfilzomib (PR-171) (combination with valproic acid)")] <- NA # TO BE DELETED ONCE CORRECTED

names <- unlist(strsplit(conc.name, "///"))

for(j in 1:length(names)){
  drug_with_ids_gbm[1+nrow(drug_with_ids_gbm) , ] = drug_with_ids_gbm[ind, ]
  drug_with_ids_gbm$GBM.drugid[nrow(drug_with_ids_gbm)] = names[j]
}

drug_with_ids_gbm <- drug_with_ids_gbm[-ind,]
drug_with_ids_gbm$unique.drugid <- ifelse(grepl("Doxorubicin", drug_with_ids_gbm$GBM.drugid), drug_with_ids_gbm$GBM.drugid, drug_with_ids_gbm$unique.drugid) 


# ======================== Drug object ========================
drugs<- read.xlsx(paste(input_dir, "mmc3.xlsx", sep="") ,rowNames = TRUE , startRow = 2)
drugs[duplicated(drugs[ , "Compound.name"]),] #Checking duplications in drug names: 3 duplications found

#Replacing duplicated drug names in "Compound.name" column with the correct drug names:
drugs$Compound.name[drugs$Molecular.Formula=="C26H27ClN2O"]<-"LOFEPRAMINE"
drugs$Compound.name[drugs$`Chemical.Abstracts.Service.(CAS).code`=="101477-54-7"]<-"Lomerizine 2hcl"
drugs$Compound.name [drugs$ Compound.name == "lomerizine"] <-"Mitoxantrone dihydrochloride" 

#To make sure that the drug-object contains ALL the drug names used in screen2 and screen3 
#the drug-names in drug-object are mapped to drug-names from screen2 and screen3.

drugs$Compound.name[drugs$Compound.name == "bortezomib [velcade ]"]<-"Bortezomib  "
drugs$Compound.name[drugs$Compound.name == "5-azacytidine"]<-"Azacytidine-5"
drugs$Compound.name[drugs$Compound.name == "TV-001"]<-"prm-116"
drugs$Compound.name[drugs$Compound.name == "TV-002"]<-"prm-122"
drugs$Compound.name[drugs$Compound.name == "TV-003"]<-"PRM-123"
drugs$Compound.name[drugs$Compound.name == "tenovin-6 derivative 39"] <- "Tenovin-6 derivat A"#To be double checked
drugs$Compound.name[drugs$Compound.name == "tenovin-6 derivative 50"] <- "Tenovin-6 derivat B"#To be double checked 
drugs$Compound.name[drugs$Compound.name == "GLN-1001"]<-"Vakuinol-1" # Based on molecular formula (ClC1=CC=C(C2=NC(C=CC=C3)=C3C(C(O)C4NCCCC4)=C2)C=C1), The correct dictation is "Vacquinol-1"
drugs$Compound.name[drugs$Compound.name == "AM404"]<-"AMA404"
drugs$Compound.name[drugs$Compound.name == "dichlorbenzamide"]<-"Dichlorbenzamil"#To be double checked 
drugs$Compound.name[drugs$Compound.name == "8-azaguanine"]<-"azaguanine-8"
drugs$Compound.name[drugs$Compound.name == "duloxetine"]<-"(R)-Duloxetine" # this is mapped because we know drugs from "drug"file and "dose-resp" file must be the same
drugs$Compound.name[drugs$Compound.name == "propranolol-(S)"]<-"(S)-propranolol" # this is mapped because we know drugs from "drug"file and "dose-resp" file must be the same

# Cleaning the drug names for mapping
drugs<-drug_name_correction(table=drugs , column = "Compound.name") 
drug_with_ids_gbm<-drug_name_correction(table = drug_with_ids_gbm, column = "GBM.drugid")
setdiff(drug_with_ids_gbm$clean.ids , drugs$clean.ids)
setdiff(drugs$clean.ids , drug_with_ids_gbm$clean.ids)

# Drugs_with_ids$GBM.drugid includes drugs from both screen2 and screen3, there is no info available for some of these drugs in the drugs df
# These drugs are added to the drug object. 
drugs <- merge(drug_with_ids_gbm[!is.na(drug_with_ids_gbm$GBM.drugid), c("GBM.drugid", "clean.ids", "unique.drugid")],drugs, by="clean.ids",all.x = TRUE)
rownames(drugs) <- drugs$unique.drugid
drugs <- drugs[, 3:ncol(drugs)]

# ======================== Sensitivity data ======================== 
#Published AUC info
drug_cell<-read.delim(paste(input_dir , "HGCC_drug_response_AUC.txt", sep=""), header=T, sep="\t", stringsAsFactors = FALSE)
colnames(drug_cell)[colnames(drug_cell)=="X"]<-"GBM.drugid"
drug_cell<-drug_name_correction(table= drug_cell, column="GBM.drugid")
temp<-drug_name_correction(table= drug_with_ids_gbm, column="GBM.drugid")
temp<-merge(temp[,c("unique.drugid","clean.ids")], drug_cell, by="clean.ids", all.y=T)
drug_cell<-data.frame(temp[,-c(1,3)])
rm(temp)


#Transposing from wide to long 
drug_cell_long<- melt(setDT(drug_cell), id.vars = "unique.drugid", variable.name = "Pat")
drug_cell_long$EXP_details<-paste(drug_cell_long$unique.drugid, drug_cell_long$Pat , sep ="_")

# ======================== Screens ========================

screen_objects<-function(screen, cell_obj, drug_obj,drug_with_ids,d_c_l){
  
  #Keeping the rows that are not all NaNs (2 represents the two first columns which are drugs and pat names).
  screen<-screen[rowSums(is.na(screen[,3:ncol(screen)])) <ncol(screen)-2, ] 
  
  #Multiplying the viabilites by 100
  screen[,3:ncol(screen)]<-apply(screen[,3:ncol(screen)],2,function(x){x*100})
  
  screen<-drug_name_correction(table = screen , column = "Drug")
  temp<-drug_name_correction(table= drug_with_ids_gbm, column="GBM.drugid")
  temp<-merge(temp[,c("unique.drugid","clean.ids")], screen, by="clean.ids", all.y=T)
  #temp$unique.drugid <- ifelse(is.na(temp$unique.drugid) & grepl("Doxorubicin", temp$Drug), temp$Drug , temp$unique.drugid)
  screen<-data.frame(temp[-c(1,3)])
  rm(temp)
  
  screen$EXP_details<-paste(screen$unique.drugid, screen$Pat , sep ="_")
  colnames(screen) <- c("drugid","cellid", colnames(screen)[3:length(colnames(screen))])# With this, drugids in the screen are already unique.ids
  rownames(screen)<-screen$EXP_details
  
  
  sen_info_scr<-screen[,c(1,2,14)]
  sen_info_scr$Duration<-"3 Days"
  
  ###sensitivity_raw
  dose_scr<-as.data.frame(matrix(nrow=nrow(screen)))
  i=1
  for (col in colnames(screen)[3:13]){
    micro_mol<-as.numeric(substr(col, start = 6, stop = nchar(col)))
    name<-paste("doses",i,sep="") #doses need to be named as (doses1, doses2,..)
    dose_scr[,name]<-micro_mol
    i=i+1}
  
  dose_scr<-dose_scr[-1]
  rownames(dose_scr)<-rownames(screen)
  viability_scr<-screen[,3:13] #Only keeping the doses 
  colnames(viability_scr)<-colnames(dose_scr) 
  
  sen_raw_scr<-abind(dose_scr, viability_scr, along = 3) #Creating a 3d array
  dimnames(sen_raw_scr)[[3]] <- c("Dose","Viability")
  
  ###Sensitivity_profile
  print("Running calculateFromRaw function... This takes a while")
  AAC_IC50_scr<-PharmacoGx:::.calculateFromRaw(sen_raw_scr, cap = NA, nthread = 1, family = c("normal","Cauchy"), scale = 0.07, n=1 )
  sen_profile_scr<-cbind(as.data.frame(AAC_IC50_scr["AUC"]), as.data.frame(AAC_IC50_scr["IC50"]))
  colnames(sen_profile_scr)<-c("aac_recomputed","ic50_recomputed")
  
  sen_profile_scr$auc_published<-d_c_l$value[match(rownames(sen_profile_scr),d_c_l$EXP_details)]
  sen_profile_scr$auc_published<-sen_profile_scr$auc_published
  
  ########### Adapting the cell object for each screen
  scr_cell_obj<-cell_obj
  diff_scr<-unique(screen$cellid[screen$cellid %in% scr_cell_obj$cellid ==FALSE])
  
  #For i in diff_scr, adds a new row 
  for(i in 1:length(diff_scr)){
    cid=substring(diff_scr[i],1,5) 
    if (cid %in% scr_cell_obj$cellid){
      scr_cell_obj[nrow(scr_cell_obj)+1,]<- scr_cell_obj[scr_cell_obj$cellid==cid,][1,]
      scr_cell_obj[nrow(scr_cell_obj),2]<- substring(diff_scr[i],6) #Replicate name
      scr_cell_obj[nrow(scr_cell_obj),4]<- diff_scr[i] #cell id
    }
    else{scr_cell_obj[nrow(scr_cell_obj)+1,]<- c(NA, NA,NA, diff_scr[i], rep(NA,length(colnames(scr_cell_obj))-4))}
  }
  
  rownames(scr_cell_obj)<-scr_cell_obj$cellid
  scr_cell_obj$tissueid <- "GBM"  # Adding mandatory tissue Id column 
  scr_cell_obj[,2]<-gsub("_","",scr_cell_obj[,2],fixed = T)# Removing "_"
  
  ########## Adapting the drug object for each screen
  scr_drugs<-drugs[unique(screen$drugid),]
  
  ########## Curation Cell dataframe
  scr_cur_cell <- data.frame(unique.cellid=rownames(scr_cell_obj),
                             GBM.cellid=rownames(scr_cell_obj),
                             row.names=rownames(scr_cell_obj))
  
  ########## Curation drug dataframe
  scr_cur_drug <- data.frame(drugid = scr_drugs$drugid,
                             GBM.drugid = scr_drugs$GBM.drugid,
                             row.names = scr_drugs$drugid)
  
  ########## Curation tissue dataframe
  scr_cur_tissue <- data.frame(data.frame(tissueid = rep("GBM", nrow(scr_cell_obj)),
                                          GBM.tissueid = rep("GBM", nrow(scr_cell_obj)),
                                          row.names = scr_cell_obj$cellid))
  
  
  return(list("sen_info"= sen_info_scr, "sen_raw"= sen_raw_scr, "sen_profile"= sen_profile_scr,
              "cell_obj"= scr_cell_obj, "drug_obj"= scr_drugs,
              "cur_cell"= scr_cur_cell, "cur_drug"= scr_cur_drug,"cur_tissue"= scr_cur_tissue))
}

# ======================== Creating PSet ========================
# =============Screen2 =============

screen2<-read.delim(paste(input_dir, "Screen2-drugData.txt", sep=""), stringsAsFactors = FALSE)#Drug_dose_scr2 response data

#Doxirubicin is a typo of Doxorubicin
screen2$Drug[screen2$Drug == "Doxirubicin1"]<-"Doxorubicin1"
screen2$Drug[screen2$Drug == "Doxirubicin2"]<-"Doxorubicin2"
screen2$Drug[screen2$Drug == "Doxirubicin3"]<-"Doxorubicin3"
screen2$Drug[screen2$Drug == "Doxirubicin4"]<-"Doxorubicin4"
screen2$Drug[screen2$Drug == "Doxirubicin5"]<-"Doxorubicin5"
screen2$Drug[screen2$Drug == "Doxirubicin6"]<-"Doxorubicin6"
screen2$Drug[screen2$Drug == "Doxirubicin7"]<-"Doxorubicin7"
screen2$Drug[screen2$Drug == "Doxirubicin8"]<-"Doxorubicin8"

scr2_objects<-screen_objects(screen=screen2, cell_obj=all_cell_obj, drug_obj=drugs, drug_with_ids=drug_with_ids_gbm, d_c_l=drug_cell_long)

GBM_scr2_PSet<- PharmacoGx::PharmacoSet("GBM_scr2_PSet",
                                        molecularProfiles = list( "rna" = expression_SE, "rna_probe" = expression_probe_SE, 
                                                                  "rna_ruv"=expression_ruv_SE, "cnv"= cnv_SE, "mut" = mutation_SE, 
                                                                  "methyl_probe" = methyl_SE, "methyl_gene" = methyl_gene_SE),
                                        
                                        cell = scr2_objects[["cell_obj"]],
                                        drug = scr2_objects[["drug_obj"]],
                                        sensitivityInfo = scr2_objects[["sen_info"]],
                                        sensitivityRaw =  scr2_objects[["sen_raw"]],
                                        sensitivityProfiles = scr2_objects[["sen_profile"]],
                                        curationDrug = scr2_objects[["cur_drug"]],
                                        curationCell = scr2_objects[["cur_cell"]],
                                        curationTissue = scr2_objects[["cur_tissue"]],
                                        datasetType = "sensitivity",
                                        verify = TRUE)

GBM_scr2_PSet@annotation$notes <- "This PSet includes drug-dose information from phase II screening of the paper. 1. All cellids in the PSet have prefix of 'U' and suffix of 'MG' (expect for 'human_astrocytes'). 2. Types of mutations affecting mutant cells are concatenated by '///' in the assay data of mutation ESet from 'molecular-profiles' object. 3. All cell and drug metadata can be found in 'cell' and 'drug' objects, respectively. 4. Dose values are based on micromolar. 5. Throughout the 'sensitivity' object, a unique identifier has been created by concatenating drugid-cellid. 6. All raw dose and viability values are in the 'sensitivity-raw' object. 7. 'sensitivity-profiles' includes published-AUC, recomputed_AAC, and recomputed_IC50."

# ============= Screen3 =============
screen3<-read.delim(paste(input_dir,"Screen3-drugData.txt", sep = ""), stringsAsFactors = FALSE)#Drug_dose_scr2 response data

#"Vinorelbinetartrate" is differently spelled in GBM_scr2 (with space) and GBM_scr3 (without space)
#To avoid defining a unique id for a same drug it is manually added to scr3_cur_drug through the below line.
screen3$Drug[screen3$Drug == "Vinorelbinetartrate"] <- "Vinorelbine tartrate"

scr3_objects<-screen_objects(screen=screen3, cell_obj=all_cell_obj, drug_obj=drugs, drug_with_ids=drug_with_ids_gbm, d_c_l=drug_cell_long)

#Removing the published AUC values for the cell-drug pairs mutual between screen2 and screen3
#The published AUCs have a higher correlation with calculated AUCs from screen2
#Since the original paper has not specified what screen do the published AUCs belong to
#It is more conservative to report the published AUCs for screen2 

Mutual_pairs<-merge(scr2_objects[["sen_profile"]],scr3_objects[["sen_profile"]], by="row.names")
cor(Mutual_pairs$aac_recomputed.x , Mutual_pairs$auc_published.x , use= "complete.obs")##Correlation between auc_published and aac_recomputed from screen2 (-0.9132723)
cor(Mutual_pairs$aac_recomputed.y , Mutual_pairs$auc_published.y , use= "complete.obs")##Correlation between auc_published and aac_recomputed from screen3 (-0.6995755)

scr3_objects[["sen_profile"]]$auc_published[rownames(scr3_objects[["sen_profile"]]) %in% Mutual_pairs$Row.names]<-NA


GBM_scr3_PSet<- PharmacoGx::PharmacoSet("GBM_scr3_PSet",
                                        molecularProfiles = list( "rna" = expression_SE, "rna_probe" = expression_probe_SE,
                                                                  "rna_ruv" = expression_ruv_SE,"cnv"= cnv_SE, "mut" = mutation_SE,
                                                                  "methyl_probe" = methyl_SE, "methyl_gene" = methyl_gene_SE),
                                        
                                        
                                        cell = scr3_objects[["cell_obj"]],
                                        drug = scr3_objects[["drug_obj"]],
                                        sensitivityInfo = scr3_objects[["sen_info"]],
                                        sensitivityRaw = scr3_objects[["sen_raw"]],
                                        sensitivityProfiles <- scr3_objects[["sen_profile"]],
                                        curationDrug = scr3_objects[["cur_drug"]],
                                        curationCell = scr3_objects[["cur_cell"]],
                                        curationTissue = scr3_objects[["cur_tissue"]],
                                        datasetType = "sensitivity",
                                        verify = TRUE)

GBM_scr3_PSet@annotation$notes <- "This PSet includes drug-dose information from phase III screening of the paper. 1. All cellids in the PSet have prefix of 'U' and suffix of 'MG' (expect for 'human_astrocytes'). 2. Types of mutations affecting mutant cells are concatenated by '///' in the assay data of mutation ESet from 'molecular-profiles' object. 3. All cell and drug metadata can be found in 'cell' and 'drug' objects, respectively. 4. Dose values are based on micromolar. 5. Throughout the 'sensitivity' object, a unique identifier has been created by concatenating drugid-cellid. 6. All raw dose and viability values are in the 'sensitivity-raw' object. 7. 'sensitivity-profiles' includes published-AUC, recomputed_AAC, and recomputed_IC50. 8. Numbers in 'replicate' column from the 'cell' object are not interpretable as there are merely dummy numbers emphasizing that the cell line is a replicate."  

saveRDS(GBM_scr3_PSet, paste0(out_dir, "GBM.rds"))