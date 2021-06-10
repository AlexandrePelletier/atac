library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
source("../methyl/scripts/utils/new_utils.R")
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 50000 * 1024^2) # for 50 Gb RAM

out<-"analyses/cbps_merged"

cbps<-readRDS(fp(out,"cbps_atac1-4_merged_qc_rna.rds"))

DefaultAssay(cbps) <- 'ATAC'
Idents(cbps) <- 'inferred.lineage'

da_peaks <- FindAllMarkers(
  object = cbps,
  min.pct = 0.2,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)

fwrite(data.table(da_peaks),fp(out,"lineage_specific_peaks.csv.gz"),sep=";")