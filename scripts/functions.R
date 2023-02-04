library(SummarizedExperiment)
library(stringr)
library(data.table)
library(BiocGenerics)

# ======================== drug_name_correction ========================
# Function to find the synonym drug names by not considering the salts and unwanted punctuation
drug_name_correction <- function(table, column) {
  table$clean.ids <- tolower(str_replace_all(table[, column], "[^[:alnum:]]", ""))

  salts <- c(
    "malate", "sulfate", "dihydrochloride", "hydrochloride", "citrate", "ethanolamine",
    "2hcl", "oxalate", "sodiumsalt", "bromide", "monohydrate", "isethionate", "sodium",
    "furoate", "hcl", "ca", "dimaleate", "oxalate", "dihydrate", "maleate", "fumarate", "lactate", "di"
  )
  for (salt in salts) {
    table[, "clean.ids"] <- gsub(salt, "", table[, "clean.ids"], ignore.case = T)
  }
  return(table)
}

# ======================== fdata_builder ========================
# Function to annotate the genes when creating feature data for a PSet. If annotation for a gene is not available
# it will not be removed but its annotation will be presented by "NA".

fdata_builder <- function(annotation, assay, ID_column = "gene_name") {
  feat_empty <- data.frame(matrix(NA, nrow = nrow(assay), ncol = ncol(annotation)))
  colnames(feat_empty) <- colnames(annotation)
  feat_empty[, ID_column] <- rownames(assay)

  annotated <- annotation[annotation[, ID_column] %in% rownames(assay), ] # Subseting the genes from features_gene that belong to cnv
  feat <- merge(x = feat_empty, y = annotated, all = TRUE) # We keep all the genes even the ones that don't have annotations.
  feat <- feat[!duplicated(feat[, ID_column]), ]

  # Reformatting feature data to match the PSet requirements:
  feat <- feat[match(rownames(assay), feat[, ID_column]), ]
  rownames(feat) <- feat[, ID_column]
  feat$Symbol <- feat$gene_name
  return(feat)
}

# ======================== ph_data_builder ========================
# Function to create the pheno data from cell-line object (given it includes the meta data)

ph_data_builder <- function(annotation, assay) {
  phen <- as.data.frame(colnames(assay))
  phen$temp <- sub("U", "", phen[, 1])
  phen$temp <- substring(phen$temp, 1, 4)
  annotation$temp <- substring(rownames(annotation), 2, 5)
  phen <- merge(phen, annotation, by = "temp", all.x = TRUE)

  phen$batchid <- NA
  phen$cellid <- phen$`colnames(assay)`
  phen$Replicate <- substring(phen$cellid, 6)
  phen$Replicate[phen$Replicate == ""] <- NA
  rownames(phen) <- phen$`colnames(assay)`

  phen <- phen[, c(15, 18, 16, 17, 9:13)]
  return(phen)
}


# ======================== eSetToSE ========================
# A function converting ExpressionSet to SummarizedExperiment

eSetToSE <- function(eSet, annot_name) {
  BiocGenerics::annotation(eSet) <- annot_name
  stopifnot(all(rownames(fData(eSet)) == rownames(exprs(eSet))))
  stopifnot(all(rownames(pData(eSet)) == colnames(exprs(eSet))))

  # Build summarized experiment from eSet
  SE <- SummarizedExperiment::SummarizedExperiment(
    assays = SimpleList(as.list(Biobase::assayData(eSet))),
    # Switch rearrange columns so that IDs are first, probes second
    rowData = S4Vectors::DataFrame(Biobase::fData(eSet)),
    colData = S4Vectors::DataFrame(Biobase::pData(eSet)),
    metadata = list(
      "experimentData" = eSet@experimentData,
      "annotation" = Biobase::annotation(eSet),
      "protocolData" = Biobase::protocolData(eSet)
    )
  )
  # Extract names from expression set
  SummarizedExperiment::assayNames(SE) <- Biobase::assayDataElementNames(eSet)

  stopifnot(all(rownames(colData(SE)) == rownames(pData(eSet))))
  stopifnot(all(rownames(rowData(SE)) == rownames(fData(eSet))))
  return(SE)
}
