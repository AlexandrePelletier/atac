#Can we pool atac samples based on sex differences ?
source("../methyl/scripts/utils/new_utils.R")
library(Seurat)
library(Signac)

cbps<-readRDS("outputs/cbps_merged/cbps_atac1-4_merged_qc.rds")
cbps

head(rownames(cbps))

sum(str_detect(rownames(cbps),'chrY'))#152 ChrY site
chry_peaks<-rownames(cbps)[str_detect(rownames(cbps),'chrY')]
q95<-quantile(rowSums(cbps@assays$ATAC@counts[chry_peaks,]>0),0.95)
top_y_peaks<-chry_peaks[rowSums(cbps@assays$ATAC@counts[chry_peaks,]>0)>=q95]
VlnPlot(cbps,features =top_y_peaks,group.by="dataset" )

cbps$chry_counts<-colSums(cbps@assays$ATAC@counts[chry_peaks,])
VlnPlot(cbps,features ="chry_counts",group.by="dataset" )


cbps$top_chry_counts<-colSums(cbps@assays$ATAC@counts[top_y_peaks,])
VlnPlot(cbps,features ="top_chry_counts",group.by="dataset" )

q90<-quantile(rowSums(cbps@assays$ATAC@counts[chry_peaks,]>0),0.90)
top90_y_peaks<-chry_peaks[rowSums(cbps@assays$ATAC@counts[chry_peaks,]>0)>=q90]
cbps$top90_chry_counts<-colSums(cbps@assays$ATAC@counts[top90_y_peaks,])
VlnPlot(cbps,features ="top90_chry_counts",group.by="dataset" )


q90_var_y_peaks<-quantile(apply(cbps@assays$ATAC@counts[chry_peaks,],1,sd),0.90)
tops_var_y_peaks<-chry_peaks[apply(cbps@assays$ATAC@counts[chry_peaks,],1,sd)>q90_var_y_peaks]
cbps$top_var_chry_counts<-colSums(cbps@assays$ATAC@counts[tops_var_y_peaks,])
VlnPlot(cbps,features ="top_var_chry_counts",group.by="dataset" )

cbps$chry_features<-colSums(cbps@assays$ATAC@counts[chry_peaks,]>0)
VlnPlot(cbps,features ="chry_features",group.by="dataset" )+geom_hline(yintercept = 3,color="red")


#%efficiency
cells_male<-str_detect(cbps$dataset,"1|4")
efficiency<-sum(cbps@meta.data[cells_male,]$chry_features>3)/sum(cells_male)
efficiency #91.3%

#%specificity
specificity<-sum(cbps@meta.data[!cells_male,]$chry_features>3)/sum(!cells_male)
specificity #3.3%

head(cbps@meta.data)

#femelles specific peaks ? 
female_peaks<-FindMarkers(cbps,group.by = "dataset",ident.1 = c("cbp.atac2","cbp.atac3"),only.pos = T)
female_peaks<-data.table(female_peaks,keep.rownames = "gene")
head(female_peaks,100)
female_peaks<-female_peaks[pct.2<0.02&p_val_adj]

top_femal_peaks<-female_peaks[p_val_adj==0&avg_log2FC>0]$gene


cbps$femal_spe_chrX_peaks_counts<-colSums(cbps@assays$ATAC@counts[top_femal_peaks,])
VlnPlot(cbps,features ="femal_spe_chrX_peaks_counts",group.by="dataset" )+geom_hline(yintercept = 1,color="red")

#%efficiency
cells_female<-str_detect(cbps$dataset,"2|3")
efficiency<-sum(cbps@meta.data[cells_female,]$femal_spe_chrX_peaks_counts>0)/sum(cells_female)
efficiency #89%

#%specificity
specificity<-sum(cbps@meta.data[!cells_female,]$femal_spe_chrX_peaks_counts>0)/sum(!cells_female)
specificity #5%


# test#1 OR test#2 ?
efficiency<-sum(cbps@meta.data[cells_male,]$chry_features>3|cbps@meta.data[cells_male,]$femal_spe_chrX_peaks_counts==0)/sum(cells_male)
efficiency #98.6%


specificity<-sum(cbps@meta.data[!cells_male,]$chry_features>3|cbps@meta.data[!cells_male,]$femal_spe_chrX_peaks_counts==0)/sum(!cells_male)
specificity #13.6% 



# test#1 AND test#2 ?
efficiency<-sum(cbps@meta.data[cells_male,]$chry_features>3&cbps@meta.data[cells_male,]$femal_spe_chrX_peaks_counts==0)/sum(cells_male)
efficiency #87%

specificity<-sum(cbps@meta.data[!cells_male,]$chry_features>3&cbps@meta.data[!cells_male,]$femal_spe_chrX_peaks_counts==0)/sum(!cells_male)
specificity #0.2% 

#cells designation 1 : male chry>3, female X >0
cbps$male<-cbps$chry_features>3
cbps$female<-cbps$femal_spe_chrX_peaks_counts>0
cbps$ambigous<-cbps$male&cbps$female
sum(cbps$ambigous/ncol(cbps)) #3.5%
cbps$undetermined<-!(cbps$male|cbps$female)
sum(cbps$undetermined/ncol(cbps)) #8.8%


#cells designation : male if chry>2, female if Xpeaks_count>0 | chry==0 sinon unassigned 
cbps$male<-cbps$chry_features>2
cbps$female<-cbps$femal_spe_chrX_peaks_counts>0|cbps$chry_features==0
cbps$ambigous<-cbps$male&cbps$female
sum(cbps$ambigous/ncol(cbps)) #5.6%
cbps$undetermined<-!(cbps$male|cbps$female)
sum(cbps$undetermined/ncol(cbps)) #3.4%
