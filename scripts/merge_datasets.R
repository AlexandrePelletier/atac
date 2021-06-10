
#merge datasets and peaks (ref = https://satijalab.org/signac/articles/merging.html)
library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
out<-"analyses/cbps_merged"
source("../methyl/scripts/utils/new_utils.R")
dir.create(out)
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 50000 * 1024^2) # for 50 Gb RAM
# read in peak sets
peaks.atac1 <- read.table(
  file = "datasets/run_542_atac1/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.atac2 <- read.table(
  file = "datasets/run_542_atac2/peaks.bed",
  col.names = c("chr", "start", "end")
  )
peaks.atac3 <- read.table(
  file = "datasets/run_542_atac3/peaks.bed",
  col.names = c("chr", "start", "end")
  )

peaks.atac4 <- read.table(
  file = "datasets/run_542_atac4/peaks.bed",
  col.names = c("chr", "start", "end")
  )

# convert to genomic ranges
gr.atac1 <- makeGRangesFromDataFrame(peaks.atac1)
gr.atac2 <- makeGRangesFromDataFrame(peaks.atac2)
gr.atac3 <- makeGRangesFromDataFrame(peaks.atac3)
gr.atac4 <- makeGRangesFromDataFrame(peaks.atac4)

# Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = c(gr.atac1, gr.atac2, gr.atac3, gr.atac4))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

# load metadata
md.atac1 <- read.table(
  file = "datasets/run_542_atac1/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row

md.atac2 <- read.table(
  file = "datasets/run_542_atac2/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.atac3 <- read.table(
  file = "datasets/run_542_atac3/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.atac4 <- read.table(
  file = "datasets/run_542_atac4/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

# perform an initial filtering of low count cells
md.atac1 <- md.atac1[md.atac1$passed_filters > 500, ]
md.atac2 <- md.atac2[md.atac2$passed_filters > 500, ]
md.atac3 <- md.atac3[md.atac3$passed_filters > 500, ]
md.atac4 <- md.atac4[md.atac4$passed_filters > 500, ] 

# create fragment objects
frags.atac1 <- CreateFragmentObject(
  path = "datasets/run_542_atac1/fragments.tsv.gz",
  cells = rownames(md.atac1)
)

frags.atac2 <- CreateFragmentObject(
  path = "datasets/run_542_atac2/fragments.tsv.gz",
  cells = rownames(md.atac2)
)

frags.atac3 <- CreateFragmentObject(
  path = "datasets/run_542_atac3/fragments.tsv.gz",
  cells = rownames(md.atac3)
)

frags.atac4 <- CreateFragmentObject(
  path = "datasets/run_542_atac4/fragments.tsv.gz",
  cells = rownames(md.atac4)
)

cbp.atac1.counts <- FeatureMatrix(
  fragments = frags.atac1,
  features = combined.peaks,
  cells = rownames(md.atac1)
)

cbp.atac2.counts <- FeatureMatrix(
  fragments = frags.atac2,
  features = combined.peaks,
  cells = rownames(md.atac2)
)

cbp.atac3.counts <- FeatureMatrix(
  fragments = frags.atac3,
  features = combined.peaks,
  cells = rownames(md.atac3)
)

cbp.atac4.counts <- FeatureMatrix(
  fragments = frags.atac4,
  features = combined.peaks,
  cells = rownames(md.atac4)
)

cbp.atac1_assay <- CreateChromatinAssay(cbp.atac1.counts, fragments = frags.atac1)
cbp.atac1 <- CreateSeuratObject(cbp.atac1_assay, assay = "ATAC")

cbp.atac2_assay <- CreateChromatinAssay(cbp.atac2.counts, fragments = frags.atac2)
cbp.atac2 <- CreateSeuratObject(cbp.atac2_assay, assay = "ATAC")

cbp.atac3_assay <- CreateChromatinAssay(cbp.atac3.counts, fragments = frags.atac3)
cbp.atac3 <- CreateSeuratObject(cbp.atac3_assay, assay = "ATAC")

cbp.atac4_assay <- CreateChromatinAssay(cbp.atac4.counts, fragments = frags.atac4)
cbp.atac4 <- CreateSeuratObject(cbp.atac4_assay, assay = "ATAC")


# add information to identify dataset of origin
cbp.atac1$dataset <- 'cbp.atac1'
cbp.atac2$dataset <- 'cbp.atac2'
cbp.atac3$dataset <- 'cbp.atac3'
cbp.atac4$dataset <- 'cbp.atac4'

# merge all datasets, adding a cell ID to make sure cell names are unique
combined <- merge(
  x = cbp.atac1,
  y = list(cbp.atac2, cbp.atac3, cbp.atac4),
  add.cell.ids = c("atac1", "atac2", "atac3", "atac4")
)
combined[["ATAC"]]

combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 20)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi')

saveRDS(combined,fp(out,"cbps_atac1-4_merged.rds"))

DimPlot(combined, group.by = 'dataset', pt.size = 0.1)
ggsave(fp(out,"umap_combined.png"))

CoveragePlot(
  object = combined,
  group.by = 'dataset',
  region = "chr14-99700000-99760000"
)

ggsave(fp(out,"region_check_coverage.png"))


#QC filtering
source("../methyl/scripts/utils/new_utils.R")
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)

cbps<-readRDS("analyses/cbps_merged/cbps_atac1-4_merged.rds")
cbps
# 126925 features across 232072 samples within 1 assay 
# Active assay: ATAC (126925 features, 126920 variable features)
#  2 dimensional reductions calculated: lsi, umap
head(cbps@meta.data)

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)


# change to UCSC style because the data was mapped to hg38/UCSC
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

saveRDS(annotations,"ref/gene_annotations_hg38_GRanges.rds")

# add the gene information to the object
Annotation(cbps) <- annotations

#add metadata to the object

md.atac1 <- read.table(
  file = "datasets/run_542_atac1/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row

md.atac2 <- read.table(
  file = "datasets/run_542_atac2/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.atac3 <- read.table(
  file = "datasets/run_542_atac3/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.atac4 <- read.table(
  file = "datasets/run_542_atac4/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

head(md.atac1)
head(cbps@meta.data)
rownames(md.atac1)<-paste0("atac1_",rownames(md.atac1))
rownames(md.atac2)<-paste0("atac2_",rownames(md.atac2))
rownames(md.atac3)<-paste0("atac3_",rownames(md.atac3))
rownames(md.atac4)<-paste0("atac4_",rownames(md.atac4))

head(md.atac1)
md.cbps<-Reduce(rbind,list(md.atac1,md.atac2,md.atac3,md.atac4))
head(md.cbps)

cbps <- AddMetaData(object = cbps, metadata = md.cbps)

#QC Metrics

# compute nucleosome signal score per cell
#approximate ratio of mononucleosomal to nucleosome-free fragments (stored as nucleosome_signal)
cbps <- NucleosomeSignal(object = cbps)

# compute TSS enrichment score per cell : ratio of fragments centered at the TSS to fragments in TSS-flanking regions
cbps <- TSSEnrichment(object = cbps, fast = FALSE)



#add fraction of reads in peaks (cellular sequencing depth / complexity)
cbps$pct_reads_in_peaks <- cbps$peak_region_fragments / cbps$passed_filters * 100
#and blacklist ratio (reads which are often associated with artefactual signal
cbps$blacklist_ratio <- cbps$blacklist_region_fragments / cbps$peak_region_fragments

#validate that TSS enrichment scores compute really represent enrichment for cells around TSS :
cbps$high.tss <- ifelse(cbps$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(cbps, group.by = 'high.tss') + NoLegend()
TSSPlot(cbps, group.by = 'dataset') + NoLegend()

#fragment length periodicity for all the cells :
#need mononucleosomal / nucleosome-free ratio < 4 to have a good Fragment length profile for ATACseq exp :
sum(cbps$nucleosome_signal > 4) #143
VlnPlot(cbps,"nucleosome_signal",group.by="dataset")
cbps$nucleosome_group <- ifelse(cbps$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = cbps, group.by = 'nucleosome_group')


VlnPlot(cbps,"peak_region_fragments",pt.size = 0.1,log=T,group.by="dataset")+geom_hline(yintercept = 5000)

cbps<-subset(cbps,peak_region_fragments>5000)
cbps #126925 features across 9699 samples
head(cbps@meta.data)

VlnPlot(
  object = cbps,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)

VlnPlot(cbps,"pct_reads_in_peaks",group.by="dataset",pt.size = 0.1)+geom_hline(yintercept = 15)


VlnPlot(cbps,"TSS.enrichment",pt.size = 0.1,group.by="dataset")+geom_hline(yintercept = 2)

VlnPlot(cbps,"nucleosome_signal",pt.size = 0.1,group.by="dataset")+geom_hline(yintercept = 1)

VlnPlot(cbps,"blacklist_ratio",pt.size = 0.1,group.by="dataset")+geom_hline(yintercept = 0.0015)


table(subset(
  x = cbps,
  subset = peak_region_fragments > 5000 & 
    peak_region_fragments < 60000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.0015 &
    nucleosome_signal < 1&
    TSS.enrichment > 2& TSS.enrichment <10 
)$dataset)

cbps<-subset(
  x = cbps,
  subset = peak_region_fragments > 5000 & 
    peak_region_fragments < 60000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.0015 &
    nucleosome_signal < 1&
    TSS.enrichment > 2& TSS.enrichment <10 
)

cbps #126925 features across 8733 samples 

#Normalization and linear dimensional reduction
#1) term frequency-inverse document frequency (TF-IDF) normalization.
#     across cells to correct for differences in cellular sequencing depth, 
#     and across peaks to give higher values to more rare peaks
        
cbps <- RunTFIDF(cbps)
cbps <- FindTopFeatures(cbps, min.cutoff = 'q0') #use only the top n% of features (peaks) for dimensional reduction, or remove features present in less than n cells
cbps <- RunSVD(cbps) #singular value decomposition (SVD) on the TD-IDF matrix
#TF-IDF + SVD = latent semantic indexing (LSI),

# first LSI component often captures sequencing depth (technical variation) rather than biological variation
DepthCor(cbps) #++ correl comp.1 and depth => perform downstream steps without this component

#Non-linear dimension reduction and clustering
#same than for scRNA-seq
cbps <- RunUMAP(object = cbps, reduction = 'lsi', dims = 2:30)

DimPlot(object = cbps, group.by = 'dataset') #batch effect
p1<-DimPlot(object = cbps, group.by = 'dataset') +ggtitle("unintegrated")

#Integration with Harmony
library(harmony)
#error https://github.com/satijalab/seurat/issues/1849
#need adding assay=assay.use
# and @scale.data slot

head(cbps@assays$ATAC@scale.data)

cbps <- RunHarmony(
  object = cbps,
  group.by.vars = 'dataset',
  reduction = 'lsi',
  assay.use = 'ATAC',
  project.dim = FALSE
)

# re-compute the UMAP using corrected LSI embeddings
cbps <- RunUMAP(cbps, dims = 2:30, reduction = 'harmony',reduction.name = "humap")
p2 <- DimPlot(cbps, group.by = 'dataset',reduction = "humap", pt.size = 0.1) + ggplot2::ggtitle("Harmony integration")
p1 + p2

cbps <- FindNeighbors(object = cbps, reduction = 'harmony', dims = 2:30)
cbps <- FindClusters(object = cbps, verbose = FALSE, algorithm = 3)
DimPlot(object = cbps, reduction = "humap",label = TRUE) + NoLegend()


saveRDS(cbps,fp(out,"cbps_atac1-4_merged_qc.rds"))

# labeltransfer and annotate clusters
#need first calculate geneactivity (count of fragments in gene region)
#need update signac first for geneactivity(cbps) works on merged object
# renv::install("bioc::Signac")
# #restart
# packageVersion("Signac") #1.1.1 !!

# renv::install("timoast/signac@develop")

library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
source("../methyl/scripts/utils/new_utils.R")
options(future.globals.maxSize = 50000 * 1024^2) # for 50 Gb RAM
plan("multiprocess", workers = 4)

out<-"analyses/cbps_merged"
cbps<-readRDS(fp(out,"cbps_atac1-4_merged_qc.rds"))
DefaultAssay(cbps) <- 'ATAC'
gene.activities <- GeneActivity(cbps)
sum(gene.activities>0)
saveRDS(gene.activities,fp(out,"gene_activities.rds"))
# add the gene activity matrix to the Seurat object as a new assay and normalize it
gene.activities<-readRDS(fp(out,"gene_activities.rds"))
ncol(gene.activities)
gene.activities<-gene.activities[,colnames(cbps)]
cbps[["RNA"]] <- CreateAssayObject(counts = gene.activities)

cbps <- NormalizeData(
  object = cbps,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(cbps$nCount_RNA)
)
DefaultAssay(cbps) <- 'RNA'
FeaturePlot(cbps,"GATA1",reduction = 'humap')
ggsave("analyses/cbps_merged/check_RNA_ok.png")
saveRDS(cbps,fp(out,"cbps_atac1-4_merged_qc_rna.rds"))

# Load cbps scRNA-seq data 
cbps_rna <- readRDS("../singlecell/analyses/02-hematopo_datasets_integration/cbp0-4/all_cbps.rds")
message("find transfer Anchors..")
DefaultAssay(cbps_rna)<-"RNA"
features<-SelectIntegrationFeatures(list(cbps_rna,cbps),nfeatures = 3000)
transfer.anchors <- FindTransferAnchors(
  reference = cbps_rna,
  query = cbps,
  features = features,
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
  dims = 1:30
)

predicted.labels$predicted.lineage<-predicted.labels$predicted.id
predicted.labels$predicted.id<-NULL

cbps <- AddMetaData(object = cbps, metadata = predicted.labels)
message("saving object with cell type prediction")

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


saveRDS(cbps,"analyses/cbps_merged/cbps_atac1-4_merged_qc_rna.rds")

#and annotate : 
library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
source("../methyl/scripts/utils/new_utils.R")
options(future.globals.maxSize = 50000 * 1024^2) # for 50 Gb RAM

out<-"analyses/cbps_merged"

cbps<-readRDS(fp(out,"cbps_atac1-4_merged_qc_rna.rds"))

head(cbps@meta.data)
DefaultAssay(cbps)<-"ATAC"

cbps <- FindClusters(object = cbps, algorithm = 3,resolution = 0.3)

p1<-DimPlot(object = cbps, reduction = "humap",label = TRUE,group.by = "seurat_clusters") + NoLegend()
p2<-DimPlot(object = cbps, reduction = "humap",group.by=c("predicted.lineage"))

p1+p2

mtd<-data.table(cbps@meta.data,keep.rownames = "bc")
mtd[,n.cluster:=.N,by="seurat_clusters"]
mtd[,pct.lin.cluster:=.N/n.cluster,by=.(seurat_clusters,predicted.lineage)]
p3<-ggplot(mtd)+
  geom_bar(aes(x=seurat_clusters,fill=predicted.lineage),position = 'fill')+
  scale_y_continuous(labels = scales::percent)
p1+p2+p3

p4<-DimPlot(object = cbps, reduction = "humap",group.by=c("predicted.cell_type"),label = TRUE)
mtd[,pct.ct.cluster:=.N/n.cluster,by=.(seurat_clusters,predicted.cell_type)]
p5<-ggplot(unique(mtd[seurat_clusters==5],by=c("seurat_clusters","predicted.cell_type")))+
  geom_col(aes(x=predicted.cell_type,y=pct.ct.cluster,fill=predicted.cell_type))

p1+p4+p5
Idents(cbps)<-"seurat_clusters"
cbps <- RenameIdents(
  object = cbps,
  '0' = 'HSC',
  '1'="MPP",
  '2'='Erythroid',
  '3'='LMPP',
  '4' = 'Myeloid',
  '5' = 'HSC/LMPP',
  '6'="Lymphoid",
  '8' = 'proB',
  '11' = 'DC',
  '7' = 'Mas',
  '10' = 'T cell',
  '9' = 'B cell'
)

cbps[["inferred.lineage"]]<-Idents(cbps)
cbps@meta.data$cell_type<-NULL
cbps@meta.data$infered.lineage<-NULL

DimPlot(cbps,reduction = 'humap',label = T)

mtd1<-data.table(cbps@meta.data,keep.rownames = "bc")
mtd1[,n.sample:=.N,by="dataset"]
mtd1[,pct.lin.sample:=.N/n.sample,by=.(dataset,inferred.lineage)]

ggplot(unique(mtd1[,.(inferred.lineage,pct.lin.sample,dataset)]))+
  geom_col(aes(x=dataset,y=pct.lin.sample,fill=dataset))+
  facet_wrap("inferred.lineage",scales = "free")

fwrite(mtd1,fp(out,"metadata.csv"))

saveRDS(cbps,fp(out,"cbps_atac1-4_merged_qc_rna.rds"))


#plot transfer label 
library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
source("../methyl/scripts/utils/new_utils.R")
options(future.globals.maxSize = 50000 * 1024^2) # for 50 Gb RAM

out<-"analyses/cbps_merged"

cbps<-readRDS(fp(out,"cbps_atac1-4_merged_qc_rna.rds"))

cbps_rna <- readRDS("../singlecell/analyses/02-hematopo_datasets_integration/cbp0-4/all_cbps.rds")

p1<-DimPlot(cbps_rna,group.by = "new.lineage2",label = T)+ggtitle("scRNA-seq")
p2<-DimPlot(cbps,group.by = "predicted.lineage",reduction = 'humap')+ggtitle("scATAC-seq")
p1+p2+plot_layout(guides = "collect")


#find lineage spe peaks
# change back to working with peaks instead of gene activities

DefaultAssay(cbps) <- 'ATAC'
DimPlot(cbps,reduction = 'humap',label = T)


da_peaks <- FindMarkers(
  object = cbps,
  ident.1 = "Erythroid",
  ident.2 = "Myeloid",
  min.pct = 0.2,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)
da_peaks<-data.table(da_peaks,keep.rownames = "peak")
da_peaks[p_val_adj<10^-50]
fwrite(da_peaks,fp(out,"erythroid_vs_myeloid_peaks.csv"),sep=";")
head(da_peaks)
nrow(da_peaks)
p1 <- VlnPlot(
  object = cbps,
  features = da_peaks$peak[1],
  pt.size = 0.1,
  idents = c("Erythroid","Myeloid")
)
p2 <- FeaturePlot(
  object = cbps,
  features = da_peaks$peak[1],
  reduction = 'humap',
  pt.size = 0.2,
  label=T,
  label.size = 3
)

p1 | p2

#watch closest genes from this peaks 
open_erythroid <- da_peaks[avg_log2FC > 0.5&p_val_adj<10^-50, ]$peak
open_myeloid <- da_peaks[avg_log2FC < -0.5&p_val_adj<10^-50, ]$peak

closest_genes_erythroid <- ClosestFeature(cbps, regions = open_erythroid)
closest_genes_myeloid <- ClosestFeature(cbps, regions = open_myeloid)
head(closest_genes_erythroid)
closest_genes_erythroid<-data.table(closest_genes_erythroid)
closest_genes_erythroid[gene_name%in%c("GATA1","HBD","HBB","")]
closest_genes_erythroid$gene_name

levels(pbmc) <- c("CD4 Naive","CD4 Memory","CD8 Naive","CD8 Effector","DN T","NK CD56bright","NK CD56Dim","pre-B",'pro-B',"pDC","DC","CD14 Mono",'CD16 Mono')

CoveragePlot(
  object = cbps,
  region = da_peaks$peak[1],
  extend.upstream = 40000,
  extend.downstream = 20000
)
  
#find TF motif
# renv::install("bioc::JASPAR2020")
library(JASPAR2020)
library(TFBSTools)

# renv::install("bioc::BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome.Hsapiens.UCSC.hg38)
# Get a list of motif position frequency matrices from the JASPAR database
?getMatrixSet
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = "Homo sapiens", all_versions = FALSE)
)

# add motif information
?AddMotifs
cbps <- AddMotifs(
  object = cbps,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

top_erythroid <- da_peaks[avg_log2FC > 0.5&p_val_adj<0.05 ]$peak
nrow(top_erythroid) #
enriched.motifs <- FindMotifs(
  object = cbps,
  features = top_erythroid
)

enriched.motifs<-data.table(enriched.motifs)
enriched.motifs[,motif.name:=factor(motif.name,levels=enriched.motifs[order(pvalue)]$motif.name)]
ggplot(enriched.motifs[pvalue<10^-10][order(pvalue)])+
  geom_col(aes(x=motif.name,y=fold.enrichment,fill=pvalue))

MotifPlot(
  object = cbps,
  motifs = head(rownames(enriched.motifs))
)

saveRDS(cbps,"analyses/cbps_merged/cbps_atac1-4_merged_qc_rna.rds")

#run find_lineage_specific_peaks.R

lin_peaks<-fread("analyses/cbps_merged/lineage_specific_peaks.csv.gz")

lin_peaks[cluster=="LT-HSC",cluster:="Lymphoid"]
fwrite(lin_peaks,"analyses/cbps_merged/lineage_specific_peaks.csv.gz",sep=";")

unique(lin_peaks$cluster)
lin_peaks
lin_peaks[,lineage:=factor(cluster,levels=c("HSC","MPP",'LMPP',"Lymphoid","Myeloid","Erythroid","proB",'B cell','T cell','Mas','Mo'))]
lin_peaks[,accessibility:=ifelse(avg_log2FC>0,"increased","decreased")]

ggplot(lin_peaks[p_val_adj<0.05&lineage%in%c("Erythroid","Lymphoid","Myeloid","HSC","MPP","LMPP","proB")])+
       geom_bar(aes(x=lineage,fill=accessibility),position = "dodge")

ggplot(lin_peaks[p_val_adj<0.05&lineage%in%c("HSC","MPP")])+
       geom_bar(aes(x=lineage,fill=accessibility),position = "dodge")

#correlate with scRNA-seq

cbps_rna <- readRDS("../singlecell/analyses/02-hematopo_datasets_integration/cbp0-4/all_cbps.rds")
DimPlot(cbps_rna,label=T,group.by="new.lineage2")
DefaultAssay(cbps_rna)<-"SCT"
FeaturePlot(cbps_rna,"PBX1",max.cutoff = 2,label=T,label.size = 3,pt.size = 0.1)
VlnPlot(subset(cbps_rna,new.lineage2!="undefined"),"PBX1",group.by = "new.lineage2" ,pt.size = 0.1,log=T,sort = T)

FeaturePlot(cbps_rna,c("GATA1"),max.cutoff = 2,label=T,label.size = 3,pt.size = 0.1)

#peaks match with markers 
DimPlot(cbps_rna,group.by = 'cell_type')

cbps_rna$new.lineage3<-cbps_rna$new.lineage2
cbps_rna@meta.data[cbps_rna@meta.data$cell_type=="11-proB","new.lineage3"]<-"proB"

DimPlot(cbps_rna,group.by = 'new.lineage3')
Idents(cbps_rna)<-"new.lineage3"
markers<-FindAllMarkers(cbps_rna, min.pct = 0.3, logfc.threshold = 0.4)
source("../singlecell/scripts/utils/seurat_utils.R")
source("../singlecell/scripts/utils/scoreCluster.R")
markers<-scoreMarquageCluster(markers,cbps,seuil = "intraClusterFixe",filtreMin = 2)
markers<-annotMarkers(markers)
fwrite(markers,"../singlecell/analyses/02-hematopo_datasets_integration/cbp0-4/markers_new.lineage3.csv",sep=";")

#find closest gene from peaks
lin_peaks<-lin_peaks[,peak:=gene][,-"gene"]
closest_genes<-ClosestFeature(object = cbps,
                              regions = unique(lin_peaks$peak))
closest_genes<-data.table(closest_genes)
lin_peaks<-merge(lin_peaks,closest_genes[,gene:=gene_name][,peak:=query_region][,.(peak,gene,gene_biotype,type,distance)])
fwrite(lin_peaks,"analyses/cbps_merged/lineage_specific_peaks_and_closest_gene.csv.gz",sep=";")

unique(lin_peaks$cluster)
unique(markers$cluster)
lin_peaks[,lineage:=cluster]
lin_peaks[cluster=="DC",lineage:="Mo"]

markers[,lineage:=cluster]
markers[cluster%in%c("HSC","LT-HSC"),lineage:="HSC"]


lin_peaks_markers<-merge(lin_peaks,markers,all.x=T,by=c("gene","lineage"))
fwrite(lin_peaks_markers,"analyses/cbps_merged/lineage_specific_peaks_and_closest_gene_integr_with_scRNAseq_markers.csv.gz",sep=";")

lin_peaks_markers[,n.peak.markers:=sum(!is.na(cluster.y)),by='lineage']
lin_peaks_markers[,lineage:=factor(lineage,levels = unique(lin_peaks_markers[order(n.peak.markers)]$lineage))]
ggplot(lin_peaks_markers[!is.na(cluster.y)])+geom_bar(aes(x=lineage,fill=lineage))+scale_y_log10()
lin_peaks_markers[,correlation:=ifelse(sign(avg_log2FC.x)==sign(avg_log2FC.y),"positive","negative")]
lin_peaks_markers[correlation=="positive",candidate.regulation:=ifelse(avg_log2FC.x>0,"active upregulation","active downregulation")]
lin_peaks_markers[correlation=="negative",candidate.regulation:=ifelse(avg_log2FC.x>0,"passive upregulation","passive downregulation")]
fwrite(lin_peaks_markers,"analyses/cbps_merged/lineage_specific_peaks_and_closest_gene_integr_with_scRNAseq_markers.csv.gz",sep=";")


lin_peaks_markers[,markers.matching.peak:=!is.na(cluster.y)]
lin_peaks_markers[,lineage:=factor(lineage,levels=c("HSC","MPP",'LMPP',"Lymphoid","Myeloid","Erythroid","proB",'B cell','T cell','Mas','Mo'))]

ggplot(lin_peaks_markers[!is.na(cluster.y)][lineage%in%c('LMPP',"MPP","HSC","Erythroid","Myeloid","Lymphoid",'proB')])+
  geom_bar(aes(x=lineage,fill=candidate.regulation),position = "dodge")


lin_peaks_markers[,pct.peaks.link.to.degs:=n.peak.markers/.N,by="lineage"]
ggplot(unique(lin_peaks_markers[lineage%in%c('LMPP',"MPP","HSC","Erythroid","Myeloid","Lymphoid",'proB')],by="lineage"))+
  geom_col(aes(x=lineage,y=pct.peaks.link.to.degs))+
  scale_y_continuous(labels = scales::percent)

ggplot(lin_peaks_markers[!is.na(cluster.y)][lineage%in%c('LMPP',"MPP","HSC","Erythroid","Myeloid","Lymphoid")])+geom_bar(aes(x=lineage,fill=correlation))+scale_fill_manual(values = c("red3","green3"))


ggplot(lin_peaks_markers[!is.na(cluster.y)][lineage%in%c('LMPP',"MPP","HSC","Erythroid","Myeloid","Lymphoid")])+
       geom_bar(aes(x=lineage,fill=candidate.regulation),position = "fill")+scale_y_continuous(labels = scales::percent)


ggplot(lin_peaks_markers[!is.na(cluster.y)])+geom_bar(aes(x=lineage,fill=avg_log2FC.x>0))+scale_y_log10()


#where are this peaks ?

ggplot(lin_peaks_markers[lineage%in%c('LMPP',"MPP","HSC","Erythroid","Myeloid","Lymphoid")])+
       geom_bar(aes(x=type,fill=markers.matching.peak))


                                                      
