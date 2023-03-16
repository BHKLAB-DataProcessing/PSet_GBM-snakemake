from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
S3 = S3RemoteProvider(
    access_key_id=config["key"],
    secret_access_key=config["secret"],
    host=config["host"],
    stay_on_remote=False
)

prefix = config["prefix"]

rule get_gbm_pset:
    input:
        S3.remote(prefix + "processed/expression_SE.rds"),
        S3.remote(prefix + "processed/expression_probe_SE.rds"),
        S3.remote(prefix + "processed/expression_ruv_SE.rds"),
        S3.remote(prefix + "processed/cnv_SE.rds"),
        S3.remote(prefix + "processed/mutation_SE.rds"),
        S3.remote(prefix + "processed/methyl_SE.rds"),
        S3.remote(prefix + "processed/methyl_gene_SE.rds"),
        S3.remote(prefix + "processed/scr2_objects.rds"),
        S3.remote(prefix + "processed/scr3_objects.rds")
    output:
        prefix + "PSet_GBM_scr2.rds",
        prefix + "PSet_GBM_scr3.rds"
    shell:
        """
        Rscript scripts/getGBMPSet.R {prefix}
        """

rule process_drug_screen:
    input:
        S3.remote(prefix + "download/drugs_with_ids.csv"),
        S3.remote(prefix + "download/mmc3.xlsx"),
        S3.remote(prefix + "download/HGCC_drug_response_AUC.txt"),
        S3.remote(prefix + "download/Screen2-drugData.txt"),
        S3.remote(prefix + "download/Screen3-drugData.txt"),
        S3.remote(prefix + "processed/phen_exp.rds"),
        S3.remote(prefix + "processed/phen_cnv.rds"),
        S3.remote(prefix + "processed/phen_mutation.rds"),
        S3.remote(prefix + "processed/phen_methyl.rds"),
        S3.remote(prefix + "scripts/functions.R")
    output:
        S3.remote(prefix + "processed/scr2_objects.rds"),
        S3.remote(prefix + "processed/scr3_objects.rds")
    shell:
        """
        Rscript scripts/getGBMDrugScreen.R {prefix}
        """

rule process_expression:
    input:
        S3.remote(prefix + "download/hta20hsensgcdf_24.0.0.tar.gz"),
        S3.remote(prefix + "download/Ensembl.v99.annotation.RData"),
        S3.remote(prefix + "download/GSE152160_RAW.tar"),
        S3.remote(prefix + "download/HK_genes.txt"),
        S3.remote(prefix + "processed/cell.rds"),
        S3.remote(prefix + "scripts/functions.R")
    output:
        S3.remote(prefix + "processed/expression_SE.rds"),
        S3.remote(prefix + "processed/expression_ruv_SE.rds"),
        S3.remote(prefix + "processed/expression_probe_SE.rds"),
        S3.remote(prefix + "processed/phen_exp.rds")
    shell:
        """
        Rscript scripts/getGBMGeneExpression.R {prefix}
        """

rule process_cnv:
    input:
        S3.remote(prefix + "download/Ensembl.v99.annotation.RData"),
        S3.remote(prefix + "download/HGCC_DNA_copy_number_gene_level.txt"),
        S3.remote(prefix + "processed/cell.rds"),
        S3.remote(prefix + "scripts/functions.R")
    output:
        S3.remote(prefix + "processed/cnv_SE.rds"),
        S3.remote(prefix + "processed/phen_cnv.rds")
    shell:
        """
        Rscript scripts/getGBMCNV.R {prefix}
        """

rule process_methylation:
    input:
        S3.remote(prefix + "download/Ensembl.v99.annotation.RData"),
        S3.remote(prefix + "download/HGCC_DNA_methylation.txt"),
        S3.remote(prefix + "download/MethylationEPIC_v-1-0_B2.csv"),
        S3.remote(prefix + "download/1-s2.0-S221359601630071X-mmc1.txt"),
        S3.remote(prefix + "download/1-s2.0-S221359601630071X-mmc2.txt"),
        S3.remote(prefix + "download/methMat.txt"),
        S3.remote(prefix + "processed/cell.rds"),
        S3.remote(prefix + "scripts/functions.R")
    output:
        S3.remote(prefix + "processed/methyl_gene_SE.rds"),
        S3.remote(prefix + "processed/methyl_SE.rds"),
        S3.remote(prefix + "processed/phen_methyl.rds")
    shell:
        """
        Rscript scripts/getGBMMethylation.R {prefix}
        """

rule process_mutation:
    input:
        S3.remote(prefix + "download/Ensembl.v99.annotation.RData"),
        S3.remote(prefix + "download/HGCC_WES_mutations_variants.txt"),
        S3.remote(prefix + "processed/cell.rds"),
        S3.remote(prefix + "scripts/functions.R")
    output:
        S3.remote(prefix + "processed/mutation_SE.rds"),
        S3.remote(prefix + "processed/phen_mutation.rds")
    shell:
        """
        Rscript scripts/getGBMMutation.R {prefix}
        """

rule download_script:
    output:
        S3.remote(prefix + "scripts/functions.R")
    shell:
        """
        wget 'https://raw.githubusercontent.com/BHKLAB-DataProcessing/PSet_GBM-snakemake/main/scripts/functions.R' \
            -O {prefix}scripts/functions.R
        """

rule process_cell_data:
    input:
        S3.remote(prefix + "download/mmc2.xlsx")
    output:
        S3.remote(prefix + "processed/cell.rds")
    shell:
        """
        Rscript scripts/getGBMCellData.R {prefix}
        """

rule download_data:
    output:
        S3.remote(prefix + "download/mmc2.xlsx"),
        S3.remote(prefix + "download/Ensembl.v99.annotation.RData"),
        S3.remote(prefix + "download/HGCC_WES_mutations_variants.txt"),
        S3.remote(prefix + "download/HGCC_DNA_methylation.txt"),
        S3.remote(prefix + "download/MethylationEPIC_v-1-0_B2.csv"),
        S3.remote(prefix + "download/1-s2.0-S221359601630071X-mmc1.txt"),
        S3.remote(prefix + "download/1-s2.0-S221359601630071X-mmc2.txt"),
        S3.remote(prefix + "download/methMat.txt"),
        S3.remote(prefix + "download/hta20hsensgcdf_24.0.0.tar.gz"),
        S3.remote(prefix + "download/GSE152160_RAW.tar"),
        S3.remote(prefix + "download/HK_genes.txt"),
        S3.remote(prefix + "download/HGCC_DNA_copy_number_gene_level.txt"),
        S3.remote(prefix + "download/drugs_with_ids.csv"),
        S3.remote(prefix + "download/mmc3.xlsx"),
        S3.remote(prefix + "download/HGCC_drug_response_AUC.txt"),
        S3.remote(prefix + "download/Screen2-drugData.txt"),
        S3.remote(prefix + "download/Screen3-drugData.txt")
    shell:
        """
        Rscript scripts/downloadGBMData.R {prefix}
        """
