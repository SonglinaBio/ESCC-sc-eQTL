setwd("/data/ESCC_ATAC")

library(Signac)
library(Seurat)
library(GenomicRanges)
library(ggplot2)
library(patchwork)
library(hdf5r)
library(future)
library(AnnotationHub)
library(biovizBase)
library(cowplot)
library(cols4all)
library(harmony)
library(TFBSTools)
library(JASPAR2020)
library(BSgenome.Hsapiens.UCSC.hg38)

#### Reading files and merging datasets ####
dir.ls <- c(
  P1_T = './P1_T/outs',
  P2_T = './P2_T/outs',
  P3_T = './P3_T/outs',
  SRR27532331 = './SRR27532331/outs',
  SRR27532332 = './SRR27532332/outs',
  SRR27532333 = './SRR27532333/outs',
  SRR27532334 = './SRR27532334/outs'
)

# Read peak sets and convert to genomic ranges.
for (sname in names(dir.ls)) {
  peak_file <- file.path(dir.ls[sname], "peaks.bed")
  peaks_obj <- read.table(peak_file, col.names = c("chr", "start", "end"))
  assign(paste0("peaks.", sname), peaks_obj)
  assign(paste0("gr.", sname), makeGRangesFromDataFrame(peaks_obj))
}

# Create a unified set of peaks to quantify in each dataset
combined.peaks <- disjoin(x = c(gr.P1_T, gr.P2_T, gr.P3_T,
                                gr.SRR27532331, gr.SRR27532332, gr.SRR27532333, gr.SRR27532334))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths < 10000 & peakwidths > 20]

# Create objects
obj.ls <- list()
for (i in names(dir.ls)) {
  metadata <- read.table(
    file = paste(dir.ls[[i]],'singlecell.csv',sep = '/'),
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )

  metadata <- metadata[metadata$passed_filters > 1000, ]

  frags <- CreateFragmentObject(
    path = paste(dir.ls[[i]],'fragments.tsv.gz',sep = '/'),
    cells = rownames(metadata)
  )

  counts <- FeatureMatrix(
    fragments = frags,
    features = combined.peaks,
    cells = rownames(metadata)
  )

  assay <- CreateChromatinAssay(counts, fragments = frags)
  escc <- CreateSeuratObject(assay, assay = "ATAC", meta.data = metadata)

  escc$dataset <- i
  obj.ls[[i]] <- escc
}
obj.ls

# Merge objects
combined <- merge(x = obj.ls[[1]], y = obj.ls[2:7], add.cell.ids = names(obj.ls))

# remove chromosome scaffolds and add gene annotations
peaks.keep <- seqnames(granges(combined)) %in% standardChromosomes(granges(combined))
combined <- combined[as.vector(peaks.keep), ]

# add gene annotations
ah <- AnnotationHub()
query(ah, "EnsDb.Hsapiens.v98") # retrieve record with 'object[["AH75011"]]'
ensdb_v98 <- ah[["AH75011"]]
annotations <- GetGRangesFromEnsDb(ensdb = ensdb_v98)

# change to UCSC style since the data was mapped to hg38
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"

# add the gene information to the object
Annotation(combined) <- annotations

#### QC ####
# compute nucleosome signal score per cell
combined <- NucleosomeSignal(object = combined)
combined$nucleosome_group <- ifelse(combined$nucleosome_signal > 3, 'NS > 3', 'NS < 3')

# compute TSS enrichment score per cell
combined <- TSSEnrichment(object = combined)

# add fraction of reads in peaks
combined$pct_reads_in_peaks <- combined$peak_region_fragments / combined$passed_filters * 100

# add blacklist ratio
combined$blacklist_ratio <- FractionCountsInRegion(
  object = combined, 
  assay = 'ATAC',
  regions = blacklist_hg38_unified
)

combined <- subset(
  x = combined,
  subset = nCount_ATAC > 1000 &
    nCount_ATAC < 50000 &
    pct_reads_in_peaks > 10 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 3 &
    TSS.enrichment > 2
)

#### LSI ####
# Normalization and linear dimensional reduction
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 'q5')
combined <- RunSVD(combined)

# Non-linear dimension reduction and clustering
combined <- RunUMAP(object = combined, reduction = 'lsi', dims = 2:30)
combined <- FindNeighbors(object = combined, reduction = 'lsi', dims = 2:30)
combined <- FindClusters(object = combined, 
                         resolution = 0.5,
                         verbose = FALSE, 
                         algorithm = 3 # SLM algorithm
)
colnames(combined@meta.data) <- gsub("_snn_","_combined_",colnames(combined@meta.data))

# Harmony Integration
combined <- RunHarmony(
  object = combined,
  group.by.vars = 'dataset',
  reduction.use = 'lsi',
  assay.use = 'ATAC'
)

combined <- RunUMAP(combined, reduction = "harmony", reduction.name = "umap_harmony", dims = 2:30)

combined <- FindNeighbors(object = combined, reduction = 'harmony', dims = 2:30)
combined <- FindClusters(object = combined,
                         resolution = c(0.3,0.5,1),
                         verbose = FALSE, 
                         algorithm = 3 # SLM algorithm
)
colnames(combined@meta.data) <- gsub("_snn_","_harmony_",colnames(combined@meta.data))

#### Gene Activity Matrix ####
# Create a gene activity matrix and check markers
gene.activities <- GeneActivity(combined)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
combined[['ACTIVITY']] <- CreateAssayObject(counts = gene.activities)

# normalize gene activities
DefaultAssay(combined) <- "ACTIVITY"
combined <- NormalizeData(
  object = combined,
  assay = 'ACTIVITY',
  normalization.method = 'LogNormalize',
  scale.factor = median(combined$nCount_ACTIVITY)
)

combined <- ScaleData(combined, features = rownames(combined))

#### Cell type annotation and peak calling ####
# Identifying anchors between scRNA-seq and scATAC-seq
DefaultAssay(combined) <- "ACTIVITY"

# Identify anchors
transfer.anchors <- FindTransferAnchors(reference = merge_sc, 
                                        query = combined, 
                                        features = VariableFeatures(object = merge_sc),
                                        reference.assay = "RNA", 
                                        # reference.reduction = 'harmony',
                                        query.assay = "ACTIVITY", 
                                        reduction = "cca")
                                        
# Annotate scATAC-seq cells via label transfer
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = merge_sc$celltype_SLN,
                                     weight.reduction = combined[["harmony"]], dims = 2:30)
                                     
combined <- AddMetaData(combined, metadata = celltype.predictions)

# Call peaks
DefaultAssay(combined) <- "ATAC"

peaks <- CallPeaks(
  object = combined,
  group.by = "predicted.id",
  macs2.path = "~/miniconda3/envs/R4/bin/macs3"
)

peak_data <- as.data.frame(peaks)

#### TF motif enrichment analyses ####
# Differentially accessible peaks
Idents(combined) <- "predicted.id"

da_peaks <- FindAllMarkers(
  object = combined,
  only.pos = T,
  test.use = 'wilcox',
  min.pct = 0.20
)
da_peaks_sig <- da_peaks[da_peaks$p_val_adj < 0.05 & da_peaks$avg_log2FC > 1,]

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
combined <- AddMotifs(
  object = combined,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

# Motif enrichmen
cell <- unique(da_peaks_sig$cluster)
enriched_motifs_list <- list()

for (i in cell) {
  data <- da_peaks_sig[da_peaks_sig$cluster == i,]
  peak_id <- data$gene
  
  # test enrichment
  enriched.motifs <- FindMotifs(
    object = combined,
    features = peak_id
  )
  
  enriched.motifs$celltype <- i
  enriched_motifs_list[[i]] <- enriched.motifs
}

all_enriched <- matrix(nrow = 0,ncol = 10)
colnames(all_enriched) <- colnames(enriched.motifs)
for (i in cell) {
  data <- enriched_motifs_list[[i]]
  all_enriched <- rbind(all_enriched,data)
}

sig_enriched <- all_enriched[all_enriched$p.adjust<0.05 & all_enriched$fold.enrichment>1,]
