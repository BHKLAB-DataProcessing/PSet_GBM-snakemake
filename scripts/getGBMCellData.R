library(stringr)
library(openxlsx)

args <- commandArgs(trailingOnly = TRUE)
download_dir <- paste0(args[1], "download")
processed_dir <- paste0(args[1], "processed")

cell <- read.xlsx(file.path(download_dir, "mmc2.xlsx"), rowNames = TRUE) # Cell_line names
cell$Patient_id <- rownames(cell)

saveRDS(cell, file.path(processed_dir, "cell.rds"))
