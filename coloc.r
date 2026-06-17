##sc-eqtl and ESCC GWAS colocalisation analysis 
library(coloc)
 coloc_result <- coloc.abf(
      dataset1 = list(
        snp = data$ID,
        beta = data$beta_gwas,
        varbeta = data$varbeta_gwas,
        type = "cc",
        N = 8073,
        s=0.487
      ),
      dataset2 = list(
        snp = data$ID,
        pvalues = data$pvalue_eqtl,
        beta = data$beta_eqtl,
        varbeta = data$varbeta_eqtl,
        type = "quant",
        N = 46,
        MAF = data$MAF_eqtl
      )
    )
