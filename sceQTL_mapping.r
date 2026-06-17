##matrixeqtl
rm(list=ls())
library(MatrixEQTL)
useModel = modelLINEAR

# Genotype file name  
SNP_file_name = "snp.txt"
snps_location_file_name ="snploc.txt"
snps_location_file = data.table::fread(snps_location_file_name, header=T) 
head(snps_location_file)

# Gene expression file name
expression_file_name = "expr.txt"
gene_location_file_name = "geneloc.txt"
gene_location_file = data.table::fread(gene_location_file_name, header=T) 
head(gene_location_file)

# Covariates file name
# Set to character() for no covariates
covariates_file_name = "covar.txt"

# Output file name

output_file_name_cis = tempfile()
output_file_name_tra = tempfile()

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1
pvOutputThreshold_tra = 0

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

## Load genotype data

snps = SlicedData$new()
snps$fileDelimiter = "\t"      # the TAB character
snps$fileOmitCharacters = "NA" # denote missing values;
snps$fileSkipRows = 1          # one row of column labels
snps$fileSkipColumns = 1       # one column of row labels
snps$fileSliceSize = 2000      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name)

## Load gene expression data

gene = SlicedData$new()
gene$fileDelimiter = "\t"      # the TAB character
gene$fileOmitCharacters = "NA" # denote missing values;
gene$fileSkipRows = 1          # one row of column labels
gene$fileSkipColumns = 1       # one column of row labels
gene$fileSliceSize = 2000      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name)

## Load covariates

cvrt = SlicedData$new()
cvrt$fileDelimiter = "\t"      # the TAB character
cvrt$fileOmitCharacters = "NA" # denote missing values;
cvrt$fileSkipRows = 1          # one row of column labels
cvrt$fileSkipColumns = 1       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name)
}

## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE)
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE)

me = Matrix_eQTL_main(
  snps = snps, 
  gene = gene, 
  cvrt = cvrt,
  output_file_name     = output_file_name_tra,
  pvOutputThreshold     = pvOutputThreshold_tra,
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE, 
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos, 
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = TRUE,
  min.pv.by.genesnp =TRUE,
  noFDRsaveMemory = FALSE)

plot(me)

all_cis_eqtl<-me[["cis"]][["eqtls"]]
head(all_cis_eqtl)

write.table(all_cis_eqtl,"all_cis_eqtl.txt",row.names=F,col.names = T,quote = FALSE,sep = "\t") 
saveRDS(me, file = "me.rds")
