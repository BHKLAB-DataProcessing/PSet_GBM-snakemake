from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
# S3 = S3RemoteProvider(
#     access_key_id=config["key"],
#     secret_access_key=config["secret"],
#     host=config["host"],
#     stay_on_remote=False
# )

prefix = config["prefix"]

rule get_gbm_pset:
    input:
        prefix + "processed/expression_SE.rds",
        prefix + "processed/expression_probe_SE.rds",
        prefix + "processed/expression_ruv_SE.rds",
        prefix + "processed/cnv_SE.rds",
        prefix + "processed/mutation_SE.rds",
        prefix + "processed/methyl_SE.rds",
        prefix + "processed/methyl_gene_SE.rds",
        prefix + "processed/scr2_objects.rds",
        prefix + "processed/scr3_objects.rds"
    output:
        prefix + "GBM_scr2.rds",
        prefix + "GBM_scr3.rds"
    shell:
        """
        Rscript {prefix}scripts/getGBMPSet.R {prefix}
        """

rule process_drug_screen:
    input:
        prefix + "download/drugs_with_ids.csv",
        prefix + "download/mmc3.xlsx",
        prefix + "download/HGCC_drug_response_AUC.txt",
        prefix + "download/Screen2-drugData.txt",
        prefix + "download/Screen3-drugData.txt",
        prefix + "processed/phen_exp.rds",
        prefix + "processed/phen_cnv.rds",
        prefix + "processed/phen_mutation.rds",
        prefix + "processed/phen_methyl.rds"
    output:
        prefix + "processed/scr2_objects.rds",
        prefix + "processed/scr3_objects.rds"
    shell:
        """
        Rscript {prefix}scripts/getGBMDrugScreen.R {prefix}
        """

rule process_expression:
    input:
        prefix + "download/hta20hsensgcdf_24.0.0.tar.gz",
        prefix + "download/Ensembl.v99.annotation.RData",
        prefix + "download/GSE152160_RAW.tar",
        prefix + "download/HK_genes.txt",
        prefix + "processed/cell.rds"
    output:
        prefix + "processed/expression_SE.rds",
        prefix + "processed/expression_ruv_SE.rds",
        prefix + "processed/expression_probe_SE.rds",
        prefix + "processed/phen_exp.rds"
    shell:
        """
        Rscript {prefix}scripts/getGBMGeneExpression.R {prefix}
        """

rule process_cnv:
    input:
        prefix + "download/Ensembl.v99.annotation.RData",
        prefix + "download/HGCC_DNA_copy_number_gene_level.txt",
        prefix + "processed/cell.rds"
    output:
        prefix + "processed/cnv_SE.rds",
        prefix + "processed/phen_cnv.rds"
    shell:
        """
        Rscript {prefix}scripts/getGBMCNV.R {prefix}
        """

rule process_methylation:
    input:
        prefix + "download/Ensembl.v99.annotation.RData",
        prefix + "download/HGCC_DNA_methylation.txt",
        prefix + "download/MethylationEPIC_v-1-0_B2.csv",
        prefix + "download/1-s2.0-S221359601630071X-mmc1.txt",
        prefix + "download/1-s2.0-S221359601630071X-mmc2.txt",
        prefix + "download/methMat.txt",
        prefix + "processed/cell.rds"
    output:
        prefix + "processed/methyl_gene_SE.rds",
        prefix + "processed/methyl_SE.rds",
        prefix + "processed/phen_methyl.rds"
    shell:
        """
        Rscript {prefix}scripts/getGBMMethylation.R {prefix}
        """

rule process_mutation:
    input:
        prefix + "download/Ensembl.v99.annotation.RData",
        prefix + "download/HGCC_WES_mutations_variants.txt",
        prefix + "processed/cell.rds"
    output:
        prefix + "processed/mutation_SE.rds",
        prefix + "processed/phen_mutation.rds"
    shell:
        """
        Rscript {prefix}scripts/getGBMMutation.R {prefix}
        """

rule process_cell_data:
    input:
        prefix + "download/mmc2.xlsx"
    output:
        prefix + "processed/cell.rds"
    shell:
        """
        Rscript {prefix}scripts/getGBMCellData.R {prefix}
        """

rule download_data:
    output:
        prefix + "download/mmc2.xlsx",
        prefix + "download/Ensembl.v99.annotation.RData",
        prefix + "download/HGCC_WES_mutations_variants.txt",
        prefix + "download/HGCC_DNA_methylation.txt",
        prefix + "download/MethylationEPIC_v-1-0_B2.csv",
        prefix + "download/1-s2.0-S221359601630071X-mmc1.txt",
        prefix + "download/1-s2.0-S221359601630071X-mmc2.txt",
        prefix + "download/methMat.txt",
        prefix + "download/hta20hsensgcdf_24.0.0.tar.gz",
        prefix + "download/GSE152160_RAW.tar",
        prefix + "download/HK_genes.txt",
        prefix + "download/HGCC_DNA_copy_number_gene_level.txt",
        prefix + "download/drugs_with_ids.csv",
        prefix + "download/mmc3.xlsx",
        prefix + "download/HGCC_drug_response_AUC.txt",
        prefix + "download/Screen2-drugData.txt",
        prefix + "download/Screen3-drugData.txt"
    shell:
        """
        Rscript {prefix}scripts/downloadGBMData.R {prefix}
        """
