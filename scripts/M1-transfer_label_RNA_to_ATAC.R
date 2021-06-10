
library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
source("../methyl/scripts/utils/new_utils.R")
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 50000 * 1024^2) # for 50 Gb RAM


out<-"analyses/cbps_merged"
cbps<-readRDS(fp(out,"cbps_atac1-4_merged_qc.rds"))

DefaultAssay(cbps) <- 'RNA'

# Load cbps scRNA-seq data 
cbps_rna <- readRDS("../singlecell/analyses/02-hematopo_datasets_integration/cbp0-4/all_cbps.rds")
message("find transfer Anchors..")
transfer.anchors <- FindTransferAnchors(
  reference = cbps_rna,
  query = cbps,
  reduction = 'cca'
)

message("Transfer Data..")
predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = cbps_rna$cell_type,
  weight.reduction = cbps[['harmony']],
  dims = 2:30
)

predicted.labels$predicted.cell_type<-predicted.labels$predicted.id
predicted.labels$predicted.id<-NULL
cbps <- AddMetaData(object = cbps, metadata = predicted.labels)


predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = cbps_rna$new.lineage2,
  weight.reduction = cbps[['harmony']],
  dims = 2:30
)

predicted.labels$predicted.lineage<-predicted.labels$predicted.id
predicted.labels$predicted.id<-NULL

cbps <- AddMetaData(object = cbps, metadata = predicted.labels)
message("saving object with cell type prediction")
saveRDS(cbps,"analyses/cbps_merged/cbps_atac1-4_merged_qc.rds")


p1 <- DimPlot(
  object = cbps_rna,
  group.by = 'cell_type',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')

p2 <- DimPlot(
  object = cbps,
  group.by = 'predicted.cell_type',
  reduction = "humap",
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')


p_all<-p1 |p2

ggsave(fp(out,"label_transfer_RNA_to_ATAC.png"),plot = p_all)