library(PharmacoGx)

out_dir <- "/pfs/out/" 
# out_dir <- "~/Documents/pfs/getGBMPSet/" 

expression_SE <- readRDS("/pfs/getGBMGeneExpression/expression_SE.rds")
expression_probe_SE <- readRDS("/pfs/getGBMGeneExpression/expression_probe_SE.rds")
expression_ruv_SE <- readRDS("/pfs/getGBMGeneExpression/expression_ruv_SE.rds")
cnv_SE <- readRDS("/pfs/getGBMCNV/cnv_SE.rds")
mutation_SE <- readRDS("/pfs/getGBMMutation/mutation_SE.rds")
methyl_SE <- readRDS("/pfs/getGBMMethylation/methyl_SE.rds")
methyl_gene_SE <- readRDS("/pfs/getGBMMethylation/methyl_gene_SE.rds")
scr2_objects <- readRDS("/pfs/getGBMDrugScreen/scr2_objects.rds")
scr3_objects <- readRDS("/pfs/getGBMDrugScreen/scr3_objects.rds")

# expression_SE <- readRDS("~/Documents/pfs/getGBMGeneExpression/expression_SE.rds")
# expression_probe_SE <- readRDS("~/Documents/pfs/getGBMGeneExpression/expression_probe_SE.rds")
# expression_ruv_SE <- readRDS("~/Documents/pfs/getGBMGeneExpression/expression_ruv_SE.rds")
# cnv_SE <- readRDS("~/Documents/pfs/getGBMCNV/cnv_SE.rds")
# mutation_SE <- readRDS("~/Documents/pfs/getGBMMutation/mutation_SE.rds")
# methyl_SE <- readRDS("~/Documents/pfs/getGBMMethylation/methyl_SE.rds")
# methyl_gene_SE <- readRDS("~/Documents/pfs/getGBMMethylation/methyl_gene_SE.rds")
# scr2_objects <- readRDS("~/Documents/pfs/getGBMDrugScreenObjects/scr2_objects.rds")
# scr3_objects <- readRDS("~/Documents/pfs/getGBMDrugScreenObjects/scr3_objects.rds")

# =============Screen2 =============
print("Creating GBM_scr2_PSet")
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
saveRDS(GBM_scr2_PSet, paste0(out_dir, "GBM_scr2.rds"))

# ============= Screen3 =============
print("Creating GBM_scr3_PSet")
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
saveRDS(GBM_scr3_PSet, paste0(out_dir, "GBM_scr3.rds"))

