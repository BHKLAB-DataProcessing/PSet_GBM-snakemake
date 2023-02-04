args <- commandArgs(trailingOnly = TRUE)
download_dir <- paste0(args[1], "download")
processed_dir <- paste0(args[1], "processed")

functions <- paste0(args[1], "scripts/functions.R")
source(functions)

cell <- readRDS(file.path(processed_dir, "cell.rds"))
load(file.path(download_dir, "Ensembl.v99.annotation.RData"))

# ============= Methylation data =============
# Assay data
assay_methyl <- read.delim(file.path(download_dir, "HGCC_DNA_methylation.txt"), header = T, sep = "\t", stringsAsFactors = FALSE)
illum <- read.csv(file.path(download_dir, "MethylationEPIC_v-1-0_B2.csv"), skip = 7, stringsAsFactors = FALSE) # The first 7 rows are not informative

# Removing SNP probes and sex chromosome probes
assay_methyl <- assay_methyl[-c(grep("rs", rownames(assay_methyl))), ]
assay_methyl <- assay_methyl[-c(grep("ch.X", rownames(assay_methyl))), ]

# Removing cross-reactive and polymorphic probes
polymorph_probes <- read.delim(file.path(download_dir, "1-s2.0-S221359601630071X-mmc1.txt"), stringsAsFactors = FALSE)
cross_probes <- read.delim(file.path(download_dir, "1-s2.0-S221359601630071X-mmc2.txt"), header = F, stringsAsFactors = FALSE)
assay_methyl <- assay_methyl[rownames(assay_methyl) %in% polymorph_probes$IlmnID == FALSE, ]
assay_methyl <- assay_methyl[rownames(assay_methyl) %in% cross_probes$V1 == FALSE, ]
colnames(assay_methyl)[colnames(assay_methyl) == "hAstro"] <- "human_astrocytes"

# Assay data gene-level
assay_methyl_gene <- read.delim(file.path(download_dir, "methMat.txt"), sep = "\t", stringsAsFactors = FALSE)
colnames(assay_methyl_gene)[colnames(assay_methyl_gene) == "hAstro"] <- "human_astrocytes"

# Feature data
feat_methyl <- fdata_builder(annotation = illum, assay = assay_methyl, ID_column = "IlmnID")
feat_methyl_gene <- fdata_builder(annotation = features_gene, assay = assay_methyl_gene)

# Pheno data
phen_methyl <- ph_data_builder(annotation = cell, assay = assay_methyl)
phen_methyl$Patient_id[phen_methyl$cellid == "human_astrocytes"] <- "human_astrocytes"
phen_methyl$Replicate[phen_methyl$Replicate == "_astrocytes"] <- NA

# Pheno data methylation gene level is same as Pheno data probe level

# Protocol data
protocol_methyl <- as.data.frame(rep("Infinium® MethylationEPIC BeadChip Kit", ncol(assay_methyl)), row.names = colnames(assay_methyl))
colnames(protocol_methyl) <- "Array"
protocol_methyl$Provider <- "Illumina"
protocol_methyl$URL_GpG <- "http://portal.hgcc.se/" # Link to access methylation data provided by the authors
protocol_methyl$File_GpG <- "HGCC_DNA_methylation.txt" # File including the data in the above link
protocol_methyl$URL_annotation <- "https://support.illumina.com/array/array_kits/infinium-methylationepic-beadchip-kit/downloads.html" # Link to access annotations by illumina
protocol_methyl$File_annotation <- "Infinium MethylationEPIC v1.0 B5 Manifest File (CSV Format)" # File including annotations in the above link

protocol_methyl_gene <- as.data.frame(rep("Infinium® MethylationEPIC BeadChip Kit", ncol(assay_methyl_gene)), row.names = colnames(assay_methyl_gene))
colnames(protocol_methyl_gene) <- "Array"
protocol_methyl_gene$Provider <- "Illumina"

# Creating ExpressionSet
assay_methyl <- assay_methyl[, rownames(phen_methyl)] # rearranging the colnames so it is similar to pheno data
protocol_methyl <- protocol_methyl[rownames(phen_methyl), ] # rearranging the rownames so it is similar to pheno data
methyl_eSet <- ExpressionSet(
  assayData = as.matrix(assay_methyl), phenoData = AnnotatedDataFrame(phen_methyl),
  featureData = AnnotatedDataFrame(feat_methyl), protocolData = AnnotatedDataFrame(protocol_methyl)
)

assay_methyl_gene <- assay_methyl_gene[, rownames(phen_methyl)] # rearranging the colnames so it is similar to pheno data
protocol_methyl_gene <- protocol_methyl_gene[rownames(phen_methyl), ] # rearranging the rownames so it is similar to pheno data
methyl_gene_eSet <- ExpressionSet(
  assayData = as.matrix(assay_methyl_gene),
  phenoData = AnnotatedDataFrame(phen_methyl),
  featureData = AnnotatedDataFrame(feat_methyl_gene),
  protocolData = AnnotatedDataFrame(protocol_methyl_gene)
)
print("Methylation: done")

methyl_gene_SE <- eSetToSE(methyl_gene_eSet, annot_name = "methyl_gene")
methyl_SE <- eSetToSE(methyl_eSet, annot_name = "methyl_probe")

saveRDS(methyl_gene_SE, file.path(processed_dir, "methyl_gene_SE.rds"))
saveRDS(methyl_SE, file.path(processed_dir, "methyl_SE.rds"))
saveRDS(phen_methyl, file.path(processed_dir, "phen_methyl.rds"))
