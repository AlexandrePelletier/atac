#explorate Signac workflow 
#1) adapted from cbp1_hi vignette ( https://satijalab.org/signac/articles/cbp1_hi_vignette.html )
source("../methyl/scripts/utils/new_utils.R")
library(Signac)
library(Seurat)
library(GenomeInfoDb)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# BiocManager::install("EnsDb.Hsapiens.v86")
library(EnsDb.Hsapiens.v86)
set.seed(1234)

sample<-"atac1_low_depth"
out<-fp("analyses",sample)
dir.create(out)
#Peak/Cell matrix.
counts <- Read10X_h5(filename = "datasets/run_516_atac1/filtered_peak_bc_matrix.h5")
dim(counts) #88536  2778
metadata <- read.csv(
  file = "datasets/run_516_atac1/singlecell.csv",
  header = TRUE,
  row.names = 1
)
head(metadata)

#                       total duplicate chimeric unmapped lowmapq mitochondrial passed_filters
# NO_BARCODE         13357381   4713636   191229  2987172 1369785        296828        3798731
# AAACGAAAGAAACGCC-1        4         1        0        3       0             0              0
# AAACGAAAGAAAGCAG-1     3766      1222       56      102     442            28           1916
# AAACGAAAGAAAGGGT-1     4190      1361       81       67     485            23           2173
# AAACGAAAGAAATACC-1       29         1        0       15       2             1             10
# AAACGAAAGAAATCTG-1        5         0        0        3       0             0              2


#Fragments file (all fragments )

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = 'datasets/run_516_atac1/fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)

#create ChromatinAssay "peaks" in the seurat object
cbp1_hi <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

#additional information that can be contained in the ChromatinAssay:
cbp1_hi[["peaks"]]

#can call granges on a Seurat object with a ChromatinAssay set as the active assay
granges(cbp1_hi)

#See the object interaction vignette (https://satijalab.org/signac/articles/data_structures.html) 
  #for more information about the ChromatinAssay class.

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style if the data was mapped to hg19/UCSC
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

# add the gene information to the object
Annotation(cbp1_hi) <- annotations

#QC Metrics

# compute nucleosome signal score per cell
#approximate ratio of mononucleosomal to nucleosome-free fragments (stored as nucleosome_signal)
cbp1_hi <- NucleosomeSignal(object = cbp1_hi)

head(cbp1_hi@meta.data)
# compute TSS enrichment score per cell : ratio of fragments centered at the TSS to fragments in TSS-flanking regions
cbp1_hi <- TSSEnrichment(object = cbp1_hi, fast = FALSE)

#add fraction of reads in peaks (cellular sequencing depth / complexity)
cbp1_hi$pct_reads_in_peaks <- cbp1_hi$peak_region_fragments / cbp1_hi$passed_filters * 100
#and blacklist ratio (reads which are often associated with artefactual signal
cbp1_hi$blacklist_ratio <- cbp1_hi$blacklist_region_fragments / cbp1_hi$peak_region_fragments

#validate that TSS enrichment scores compute really represent enrichment for cells around TSS :
cbp1_hi$high.tss <- ifelse(cbp1_hi$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(cbp1_hi, group.by = 'high.tss') + NoLegend()

#fragment length periodicity for all the cells :
#need mononucleosomal / nucleosome-free ratio > 4 to have a good Fragment length profile for ATACseq exp :

sum(cbp1_hi$nucleosome_signal > 4) #0
#cbp1_hi$nucleosome_group <- ifelse(cbp1_hi$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = cbp1_hi, group.by = 'high.tss')

VlnPlot(
  object = cbp1_hi,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)

cbp1_hi <- subset(
  x = cbp1_hi,
  subset = peak_region_fragments > 5000 &
    peak_region_fragments < 40000 &
    pct_reads_in_peaks > 20 &
    blacklist_ratio < 0.0015 &
    nucleosome_signal < 0.75&
    TSS.enrichment > 2& TSS.enrichment <10
)
cbp1_hi
# An object of class Seurat 
# 85377 features across 2196 samples within 1 assay 
# Active assay: peaks (85377 features, 0 variable features)

saveRDS(cbp1_hi,fp(out,ps(sample,"_QC_filtered.rds")))

#Normalization and linear dimensional reduction
#1) term frequency-inverse document frequency (TF-IDF) normalization.
#     across cells to correct for differences in cellular sequencing depth, 
#     and across peaks to give higher values to more rare peaks
        
cbp1_hi <- RunTFIDF(cbp1_hi)
cbp1_hi <- FindTopFeatures(cbp1_hi, min.cutoff = 'q0') #use only the top n% of features (peaks) for dimensional reduction, or remove features present in less than n cells
cbp1_hi <- RunSVD(cbp1_hi) #singular value decomposition (SVD) on the TD-IDF matrix
#TF-IDF + SVD = latent semantic indexing (LSI),
# first LSI component often captures sequencing depth (technical variation) rather than biological variation

DepthCor(cbp1_hi) #++ correl comp.1 and depth => perform downstream steps without this component

#Non-linear dimension reduction and clustering
#same than for scRNA-seq
cbp1_hi <- RunUMAP(object = cbp1_hi, reduction = 'lsi', dims = 2:30)
cbp1_hi <- FindNeighbors(object = cbp1_hi, reduction = 'lsi', dims = 2:30)
cbp1_hi <- FindClusters(object = cbp1_hi, verbose = FALSE, algorithm = 3)
DimPlot(object = cbp1_hi, label = TRUE) + NoLegend()

#Create a gene activity matrix
#= summing the fragments intersecting the gene body and promoter region (
# gene coordinates + 2 kb upstream region

gene.activities <- GeneActivity(cbp1_hi)
# add the gene activity matrix to the Seurat object as a new assay and normalize it
cbp1_hi[['RNA']] <- CreateAssayObject(counts = gene.activities)
cbp1_hi <- NormalizeData(
  object = cbp1_hi,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(cbp1_hi$nCount_RNA)
)

DefaultAssay(cbp1_hi) <- 'RNA'
lymph_markers<-c('SELL',"LTB","VPREB1","IGHM")
myelo_markers<-c("MPO","LYZ","PRTN3")
erythro_markers<-c("PF4","HBB","GATA2","GATA1","CPA3","HDC","FCER1A")
hsc_markers<-c('SOCS3', 'ID1','ID2',"FOS","JUN",'ZFP36',"HES1")


FeaturePlot(
  object = cbp1_hi,
  features =lymph_markers,
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)

FeaturePlot(
  object = cbp1_hi,
  features =myelo_markers,
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)


FeaturePlot(
  object = cbp1_hi,
  features =erythro_markers,
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 4
)

FeaturePlot(
  object = cbp1_hi,
  features =hsc_markers,
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 4
)

#no markers so stop here.
#save
saveRDS(cbp1_hi,fp(out,ps(sample,"_QC_filtered.rds")))


#basic Signac workflow for cbp1_hi 
#explorate Signac workflow 
#1) adapted from pbmc vignette ( https://satijalab.org/signac/articles/cbp1_hi_vignette.html )
source("../methyl/scripts/utils/new_utils.R")
library(Signac)
library(Seurat)
library(GenomeInfoDb)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# BiocManager::install("EnsDb.Hsapiens.v86")
library(EnsDb.Hsapiens.v86)
set.seed(1234)

sample<-"atac1_hi_depth"
out<-fp("analyses",sample)
dir.create(out)
#Peak/Cell matrix.
counts <- Read10X_h5(filename = "datasets/run_542_atac1/filtered_peak_bc_matrix.h5")
dim(counts) #88536  2778 > 89606  2982
metadata <- read.csv(
  file = "datasets/run_542_atac1/singlecell.csv",
  header = TRUE,
  row.names = 1
)
head(metadata)

#Fragments file (all fragments )

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = 'datasets/run_542_atac1/fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)

#create ChromatinAssay "peaks" in the seurat object
cbp1_hi <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

#additional information that can be contained in the ChromatinAssay:
cbp1_hi[["peaks"]]

#can call granges on a Seurat object with a ChromatinAssay set as the active assay
granges(cbp1_hi)

#See the object interaction vignette (https://satijalab.org/signac/articles/data_structures.html) 
  #for more information about the ChromatinAssay class.

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style if the data was mapped to hg19/UCSC
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

# add the gene information to the object
Annotation(cbp1_hi) <- annotations

#QC Metrics

# compute nucleosome signal score per cell
#approximate ratio of mononucleosomal to nucleosome-free fragments (stored as nucleosome_signal)
cbp1_hi <- NucleosomeSignal(object = cbp1_hi)

# compute TSS enrichment score per cell : ratio of fragments centered at the TSS to fragments in TSS-flanking regions
cbp1_hi <- TSSEnrichment(object = cbp1_hi, fast = FALSE)

#add fraction of reads in peaks (cellular sequencing depth / complexity)
cbp1_hi$pct_reads_in_peaks <- cbp1_hi$peak_region_fragments / cbp1_hi$passed_filters * 100
#and blacklist ratio (reads which are often associated with artefactual signal
cbp1_hi$blacklist_ratio <- cbp1_hi$blacklist_region_fragments / cbp1_hi$peak_region_fragments

#validate that TSS enrichment scores compute really represent enrichment for cells around TSS :
cbp1_hi$high.tss <- ifelse(cbp1_hi$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(cbp1_hi, group.by = 'high.tss') + NoLegend()

#fragment length periodicity for all the cells :
#need mononucleosomal / nucleosome-free ratio > 4 to have a good Fragment length profile for ATACseq exp :
sum(cbp1_hi$nucleosome_signal > 4) #0
#cbp1_hi$nucleosome_group <- ifelse(cbp1_hi$nucleosome_signal > 4, 'NS > 4', 'NS < 4')

FragmentHistogram(object = cbp1_hi, group.by = 'high.tss')

VlnPlot(
  object = cbp1_hi,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)

VlnPlot(cbp1_hi,"peak_region_fragments",log=T)+geom_hline(yintercept = 5750)

cbp1_hi <- subset(
  x = cbp1_hi,
  subset = peak_region_fragments > 5750 & #instead of >5000
    peak_region_fragments < 60000 &
    pct_reads_in_peaks > 20 &
    blacklist_ratio < 0.0015 &
    nucleosome_signal < 0.75&
    TSS.enrichment > 2.5& TSS.enrichment <10 #instead of >2
)
cbp1_hi
# An object of class Seurat 
#87467 features across 2248 samples 

#instead of  85377 features across 2196 samples 


#Normalization and linear dimensional reduction
#1) term frequency-inverse document frequency (TF-IDF) normalization.
#     across cells to correct for differences in cellular sequencing depth, 
#     and across peaks to give higher values to more rare peaks
        
cbp1_hi <- RunTFIDF(cbp1_hi)
cbp1_hi <- FindTopFeatures(cbp1_hi, min.cutoff = 'q0') #use only the top n% of features (peaks) for dimensional reduction, or remove features present in less than n cells
cbp1_hi <- RunSVD(cbp1_hi) #singular value decomposition (SVD) on the TD-IDF matrix
#TF-IDF + SVD = latent semantic indexing (LSI),

# first LSI component often captures sequencing depth (technical variation) rather than biological variation
DepthCor(cbp1_hi) #++ correl comp.1 and depth => perform downstream steps without this component

#Non-linear dimension reduction and clustering
#same than for scRNA-seq
cbp1_hi <- RunUMAP(object = cbp1_hi, reduction = 'lsi', dims = 2:30)
cbp1_hi <- FindNeighbors(object = cbp1_hi, reduction = 'lsi', dims = 2:30)
cbp1_hi <- FindClusters(object = cbp1_hi, verbose = FALSE, algorithm = 3)
DimPlot(object = cbp1_hi, label = TRUE) + NoLegend()

#Create a gene activity matrix
#= summing the fragments intersecting the gene body and promoter region (
# gene coordinates + 2 kb upstream region

gene.activities <- GeneActivity(cbp1_hi)
# add the gene activity matrix to the Seurat object as a new assay and normalize it
cbp1_hi[['RNA']] <- CreateAssayObject(counts = gene.activities)
cbp1_hi <- NormalizeData(
  object = cbp1_hi,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(cbp1_hi$nCount_RNA)
)

DefaultAssay(cbp1_hi) <- 'RNA'
lymph_markers<-c('SELL',"LTB","VPREB1","IGHM")
myelo_markers<-c("MPO","LYZ","PRTN3")
erythro_markers<-c("PF4","HBB","GATA2","GATA1","CPA3","HDC","FCER1A")
hsc_markers<-c('SOCS3', 'ID1','ID2',"FOS","JUN",'ZFP36',"HES1")


FeaturePlot(
  object = cbp1_hi,
  features =lymph_markers,
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)

FeaturePlot(
  object = cbp1_hi,
  features =myelo_markers,
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)


FeaturePlot(
  object = cbp1_hi,
  features =erythro_markers,
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 4
)

FeaturePlot(
  object = cbp1_hi,
  features =hsc_markers,
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 4
)

#no markers so stop here.
#save

saveRDS(cbp1_hi,fp(out,ps(sample,"_QC_filtered.rds")))


#scATACseq compare lodepth and hidepth
source("../methyl/scripts/utils/new_utils.R")
library(Seurat)
library(Signac)

cbp1_hi<-readRDS("analyses/atac1_hi_depth/atac1_hi_depth_QC_filtered.rds")
cbp1_lo<-readRDS("analyses/atac1_low_depth/atac1_low_depth_QC_filtered.rds")


p1<-FeaturePlot(
  object = cbp1_lo,
  features =c("GATA1"),
  pt.size = 0.2,
  max.cutoff = 'q95'
)

p2<-FeaturePlot(
  object = cbp1_hi,
  features =c("GATA1"),
  pt.size = 0.2,
  max.cutoff = 'q95'
)
p1|p2


p1<-DimPlot(cbp1_lo,label = T)

p2<-DimPlot(cbp1_hi,label = T)

p1|p2

m6_lo<-FindMarkers(cbp1_lo,ident.1 = "6",logfc.threshold = 0.25,only.pos = T,min.diff.pct = 0.2)
m6_hi<-FindMarkers(cbp1_hi,ident.1 = "6",logfc.threshold = 0.25,only.pos = T,min.diff.pct = 0.2)

nrow(m6_lo) #23
nrow(m6_hi) #117

#sc ATACseq annotate and compare ctrl LGA in clusters


