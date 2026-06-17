##twas
./SPrediXcan.py \
--model_db_path data/sc_models.db \
--covariance data/covariances.txt \
--gwas_folder data/GWAS \
--gwas_file_pattern ".*txt" \
--snp_column ID \
--effect_allele_column A1 \
--non_effect_allele_column A2 \
--beta_column beta \
--pvalue_column P \
--output_file results/sc_output.csv
