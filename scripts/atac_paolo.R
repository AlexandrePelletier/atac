
#Are the singlenuclei based bulk ATAC-seq from hypocampus of mouse  brain exploitable or not ?
#Can be integrated with the bulk RNA-seq from same brain regions and condition ?

out<-"../atac/outputs/atac_paolo"
dir.create(out)

#Design of the experiments
# condition,sample
# TREATMENT,pamh1
# TREATMENT,pamh2
# TREATMENT,pamh3
# TREATMENT,pamh4
# CONTROL,pbs1
# CONTROL,pbs2
# CONTROL,pbs3
# CONTROL,pbs4


source("../methyl/scripts/utils/new_utils.R")


#Are the singlenuclei based bulk ATAC-seq from hypocampus of mouse  brain exploitable ?####
#run ATAC-seq pipeline (run by Medhi). pipeline: https://nf-co.re/atacseq/2.0/output

#PEAKS calling
peaks<-fread("~/RUN/Run_669_ATAC/Output/Bulk/bwa/mergedLibrary/macs/broadPeak/consensus/consensus_peaks.mLb.clN.bed",col.names = c("chr","start","end","peak_id","V5",'strand'))
peaks #207423 peaks
peaks[,peak.length:=end-start]
summary(peaks$peak.length)
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 131.0   303.5   564.0   731.2   995.0  6667.0 
unique(peaks$strand) #+

#peaks annotation : are the peaks enriched for TSS, known regulatory region ? based on http://homer.ucsd.edu/homer/ngs/annotation.html 
peaks_anno<-fread("~/RUN/Run_669_ATAC/Output/Bulk/bwa/mergedLibrary/macs/broadPeak/consensus/consensus_peaks.mLb.clN.annotatePeaks.txt")

#  PeakID (cmd=annotatePeaks.pl consensus_peaks.mLb.clN.bed genome.fa -gid -gtf genes.gtf -cpu 6)   Chr     Start       End Strand Peak Score Focus Ratio/Region Size
#      1:                                                                                Interval_197154  chr9  60584449  60584614      +          0                      NA
#      2:                                                                                Interval_202036  chr9 110596990 110597389      +          0                      NA
#      3:                                                                                 Interval_66699 chr15  10463075  10464068      +          0                      NA
#      4:                                                                                Interval_162400  chr6  77485439  77486035      +          0                      NA
#      5:                                                                                Interval_153615  chr5 122473472 122473872      +          0                      NA
#     ---                                                                                                                                                                   
# 207419:                                                                                 Interval_63087 chr14  70309876  70310155      +          0                      NA
# 207420:                                                                                 Interval_67941 chr15  32422471  32423407      +          0                      NA
# 207421:                                                                                 Interval_51166 chr13  43587854  43589435      +          0                      NA
# 207422:                                                                                 Interval_91029 chr18   3382598   3383858      +          0                      NA
# 207423:                                                                                 Interval_15377 chr10   6705450   6705810      +          0                      NA
#                               Annotation               Detailed Annotation Distance to TSS Nearest PromoterID     Entrez ID Nearest Unigene Nearest Refseq Nearest Ensembl
#      1:                       Intergenic                        Intergenic          -38649      9230112J17Rik 9230112J17Rik       NR_040464             NA              NA
#      2:  intron (Setd2, intron 15 of 20)  -intron (Setd2, intron 15 of 20)           27204              Nradd         Nradd       NM_026012             NA              NA
#      3: intron (Dnajc21, intron 2 of 11) -intron (Dnajc21, intron 2 of 11)            6945            Dnajc21       Dnajc21       NM_030046             NA              NA
#      4:  intron (Ctnna2, intron 8 of 18)  -intron (Ctnna2, intron 8 of 19)         -103262      6330415B21Rik 6330415B21Rik       NR_045141             NA              NA
#      5:  intron (Atp2a2, intron 5 of 20)  -intron (Atp2a2, intron 5 of 19)           28553             Atp2a2        Atp2a2    NM_001110140             NA              NA
#     ---                                                                                                                                                                   
# 207419: intron (Slc39a14, intron 7 of 9) -intron (Slc39a14, intron 7 of 9)          -20518             Ppp3cc        Ppp3cc       NM_008915             NA              NA
# 207420:  intron (Sema5a, intron 3 of 22)  -intron (Sema5a, intron 3 of 22)          178126             Sema5a        Sema5a       NM_009154             NA              NA
# 207421:                       Intergenic                        Intergenic          -27153             Rnf182        Rnf182       NM_183204             NA              NA
# 207422:              promoter-TSS (Cul2)              -promoter-TSS (Cul2)               3               Cul2          Cul2       NM_029402             NA              NA
# 207423:                       Intergenic                        Intergenic          -82940              Oprm1         Oprm1    NM_001302795             NA              NA
#             Gene Name Gene Alias Gene Description Gene Type
#      1: 9230112J17Rik         NA               NA        NA
#      2:         Nradd         NA               NA        NA
#      3:       Dnajc21         NA               NA        NA
#      4: 6330415B21Rik         NA               NA        NA
#      5:        Atp2a2         NA               NA        NA
#     ---                                                    
# 207419:        Ppp3cc         NA               NA        NA
# 207420:        Sema5a         NA               NA        NA
# 207421:        Rnf182         NA               NA        NA
# 207422:          Cul2         NA               NA        NA
# 207423:         Oprm1         NA               NA        NA

peaks_anno<-peaks_anno[,PeakID:=`PeakID (cmd=annotatePeaks.pl consensus_peaks.mLb.clN.bed genome.fa -gid -gtf genes.gtf -cpu 6)`][,-"PeakID (cmd=annotatePeaks.pl consensus_peaks.mLb.clN.bed genome.fa -gid -gtf genes.gtf -cpu 6)"]
peaks_anno

summary(peaks_anno$`Distance to TSS`)

plot(density(peaks_anno$`Distance to TSS`))

ggplot(peaks_anno[abs(`Distance to TSS`)<50000])+geom_density(aes(x=`Distance to TSS`))

ggplot(peaks_anno[abs(`Distance to TSS`)<5000])+geom_density(aes(x=`Distance to TSS`))

ggplot(peaks_anno[abs(`Distance to TSS`)<5000])+geom_histogram(aes(x=`Distance to TSS`))

table(str_extract(peaks_anno$Annotation,"[A-Za-z]+"))
#   exon Intergenic     intron   promoter        TTS 
#      11143      72575     111483       8875       3347 

peaks_anno[,anno_broad:=str_extract(Annotation,"[A-Za-z]+")]
ggplot(peaks_anno)+geom_bar(aes(x=anno_broad,fill=anno_broad))


#Difference of Accessibilty between PAMH mice (TREATMENT) and PBS mice (Control)####
#1)boolean : Is there Peaks exclusively present in one condition ? i.e. Is there open chromatin region specifically found in one mice group ?
peaks_bool<-fread("~/RUN/Run_669_ATAC/Output/Bulk/bwa/mergedLibrary/macs/broadPeak/consensus/consensus_peaks.mLb.clN.boolean.txt")


peaks_bool[(CONTROL_R1.mLb.clN.bool&CONTROL_R2.mLb.clN.bool&CONTROL_R3.mLb.clN.bool&CONTROL_R4.mLb.clN.bool)&(!TREATMENT_R1.mLb.clN.bool&!TREATMENT_R2.mLb.clN.bool&!TREATMENT_R3.mLb.clN.bool&!TREATMENT_R4.mLb.clN.bool)]

exclu_control<-peaks_bool[(CONTROL_R1.mLb.clN.bool&CONTROL_R2.mLb.clN.bool&CONTROL_R3.mLb.clN.bool&CONTROL_R4.mLb.clN.bool)&(!TREATMENT_R1.mLb.clN.bool&!TREATMENT_R2.mLb.clN.bool&!TREATMENT_R3.mLb.clN.bool&!TREATMENT_R4.mLb.clN.bool)]$interval_id
length(exclu_control) #518



exclu_pamh<-peaks_bool[(!CONTROL_R1.mLb.clN.bool&!CONTROL_R2.mLb.clN.bool&!CONTROL_R3.mLb.clN.bool&!CONTROL_R4.mLb.clN.bool)&(TREATMENT_R1.mLb.clN.bool&TREATMENT_R2.mLb.clN.bool&TREATMENT_R3.mLb.clN.bool&TREATMENT_R4.mLb.clN.bool)]$interval_id
length(exclu_pamh) #8

#without PAMH outlier (TREATMENT_R4 because low  number of peaks detected)
exclu_pamh2<-peaks_bool[(!CONTROL_R1.mLb.clN.bool&!CONTROL_R2.mLb.clN.bool&!CONTROL_R3.mLb.clN.bool&!CONTROL_R4.mLb.clN.bool)&(TREATMENT_R1.mLb.clN.bool&TREATMENT_R2.mLb.clN.bool&TREATMENT_R3.mLb.clN.bool)]$interval_id 

length(exclu_pamh2) #196

peaks_exclu<-peaks_bool[interval_id%in%c(exclu_control,exclu_pamh,exclu_pamh2)]
peaks_exclu[,condition:=ifelse(CONTROL_R1.mLb.clN.bool,"PBS","PAMH")]
peaks_exclu[,All_PAMH_samples:=condition=="PAMH"&TREATMENT_R4.mLb.clN.bool]

ggplot(peaks_exclu[(All_PAMH_samples)|condition=="PBS"])+geom_bar(aes(x=condition,fill=condition))


ggplot(peaks_exclu)+geom_bar(aes(x=condition,fill=condition))


ggplot(peaks_exclu)+geom_bar(aes(x=condition,fill=condition,col=paste(condition,All_PAMH_samples)))



peaks_anno[,interval_id:=PeakID]
peaks_exclu_anno<-merge(peaks_exclu,peaks_anno)

peaks_exclu_anno[,closest_gene:=`Gene Name`]

ggplot(peaks_exclu_anno)+geom_bar(aes(x=anno_broad,fill=anno_broad))



#is it significant ? [todo: permutation test]
#?combn
peaks_bool[(CONTROL_R1.mLb.clN.bool&TREATMENT_R1.mLb.clN.bool&CONTROL_R3.mLb.clN.bool&CONTROL_R4.mLb.clN.bool)&(!CONTROL_R2.mLb.clN.bool&!TREATMENT_R2.mLb.clN.bool&!TREATMENT_R3.mLb.clN.bool&!TREATMENT_R4.mLb.clN.bool)]$interval_id
peaks_bool[(CONTROL_R1.mLb.clN.bool&TREATMENT_R2.mLb.clN.bool&CONTROL_R3.mLb.clN.bool&CONTROL_R4.mLb.clN.bool)&(!CONTROL_R2.mLb.clN.bool&!TREATMENT_R1.mLb.clN.bool&!TREATMENT_R3.mLb.clN.bool&!TREATMENT_R4.mLb.clN.bool)]$interval_id
peaks_bool[(CONTROL_R1.mLb.clN.bool&TREATMENT_R2.mLb.clN.bool&TREATMENT_R1.mLb.clN.bool&CONTROL_R4.mLb.clN.bool)&(!CONTROL_R2.mLb.clN.bool&!CONTROL_R3.mLb.clN.bool&!TREATMENT_R3.mLb.clN.bool&!TREATMENT_R4.mLb.clN.bool)]$interval_id


#characterization of these peaks
#close to gene TSS ?
summary(abs(peaks_anno$`Distance to TSS`))
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  #     0   12713   35982   84612   94369 2203354 

summary(abs(peaks_anno[PeakID%in%exclu_control]$`Distance to TSS`))
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   #   15   15229   37524   68331   86847  705882 

peaks_anno[,exclusivity:=ifelse(PeakID%in%exclu_control,"exclu control",ifelse(PeakID%in%exclu_pamh,"exclu PAMH","found in both group"))]

ggplot(peaks_anno)+
  geom_density(aes(x=`Distance to TSS`,fill=exclusivity,col=exclusivity),alpha=0.3)+
  coord_cartesian(xlim = c(-5e5 ,5e5))

ggplot(peaks_anno)+
  geom_density(aes(x=abs(`Distance to TSS`),fill=exclusivity),alpha=0.5)+
  scale_x_log10()

ggplot(peaks_anno)+
  geom_boxplot(aes(y=abs(`Distance to TSS`),x=exclusivity,fill=exclusivity))+coord_cartesian(ylim = c(0,250000))

genes_ctrl<-peaks_anno[PeakID%in%exclu_control]$`Gene Name`
length(genes_ctrl) #518



genes_pamh<-peaks_anno[PeakID%in%exclu_pamh]$`Gene Name`
genes_pamh #"Dtna"          "1700019C18Rik" "F2rl3"         "Daam2"         "Clca3a1"       "Relt"          "Arhgef28"      "Map3k7"  


genes_pamh2<-peaks_anno[PeakID%in%exclu_pamh2]$`Gene Name`
genes_pamh2 #193

#2) count differences : Is there difference of accessibility between condtion in each peak ? 
#input: the peaks matrix (input of DESeq2)
peaks_mat<-fread("~/RUN/Run_669_ATAC/Output/Bulk/bwa/mergedLibrary/macs/broadPeak/consensus/consensus_peaks.mLb.clN.featureCounts.txt")

#                 Geneid  Chr    Start      End Strand Length TREATMENT_R3.mLb.clN.bam CONTROL_R1.mLb.clN.bam TREATMENT_R2.mLb.clN.bam CONTROL_R4.mLb.clN.bam
#      1:      Interval_1 chr1  3121448  3121586      +    139                       19                     29                       20                     14
#      2:      Interval_2 chr1  3207133  3207409      +    277                       26                     31                       12                     37
#      3:      Interval_3 chr1  3210162  3210639      +    478                       76                     66                       49                     68
#      4:      Interval_4 chr1  3212909  3213108      +    200                       28                     28                       25                     33
#      5:      Interval_5 chr1  3216864  3217723      +    860                       42                     75                       30                     47
#     ---                                                                                                                                                     
# 207419: Interval_207419 chrY 90798705 90800965      +   2261                      986                   1180                     1203                    987
# 207420: Interval_207420 chrY 90803071 90803760      +    690                      344                    460                      488                    343
# 207421: Interval_207421 chrY 90804773 90805255      +    483                      279                    348                      355                    313
# 207422: Interval_207422 chrY 90807097 90809024      +   1928                     1060                   1233                     1229                    901
# 207423: Interval_207423 chrY 90810713 90813008      +   2296                      397                    825                      695                    584


#differential accessibility analysis using DEseq2 [to do: redo DeSeq2 after Variable selection]
res_da<-fread("~/RUN/Run_669_ATAC/Output/Bulk/bwa/mergedLibrary/macs/broadPeak/consensus/deseq2/CONTROLvsTREATMENT/CONTROLvsTREATMENT.mLb.clN.deseq2.results.txt")

res_da<-res_da[,PeakID:=Geneid][,-"Geneid"]
res_da[padj<0.1] #0

res_da_anno<-merge(res_da,)

res_da[pvalue<10^-3&abs(log2FoldChange)>0.25] #109


peak_up_ctrl<-res_da[pvalue<10^-3&log2FoldChange>0.25]$PeakID
length(peak_up_ctrl) #68

peak_up_pamh<-res_da[pvalue<10^-3&log2FoldChange<(-0.25)]$PeakID
length(peak_up_pamh) #41


ggplot(res_da)+geom_point(aes(x=log2FoldChange,y=-log10(pvalue),col=-log10(pvalue)))

genes_peak_up_ctrl<-peaks_anno[PeakID %in%peak_up_ctrl]$`Gene Name`
genes_peak_up_pamh<-peaks_anno[PeakID %in%peak_up_pamh]$`Gene Name`


genes_up_exclu_ctrl<-intersect(genes_peak_up_ctrl,genes_ctrl)
genes_up_exclu_ctrl # 
#  [1] "Ankrd55"       "4930480K15Rik" "Gm833"         "A730082K24Rik" "Sec61a1"       "Stac2"         "Fars2"         "A930012L18Rik" "Col23a1"       "G630093K05Rik"
# [11] "Foxp1" 
#11/68/517
OR(set1 = genes_peak_up_ctrl,set2 = genes_ctrl,size_universe = length(unique(peaks_anno$`Gene Name`)))
#p = 5.048136e-07

# res<-OR3(querys = list(gene_up_ctrl=genes_peak_up_ctrl),terms_list = list(genes_exclu_ctrl=genes_ctrl),background = unique(peaks_anno$`Gene Name`))
# 
# res


intersect(genes_peak_up_pamh,genes_pamh) #0

intersect(genes_peak_up_pamh,genes_pamh2) #6/41





##Can it be integrated with the bulk RNA-seq from same brain regions and condition ?####
#MBH 
res_de_mbh<-fread("/disks/DATATMP/PhD_AlexandrePelletier/atac/twas_res/twas_MBH_untreatedSOPKF3_vs_untreatedControl_gene.csv.gz")


res_de_mbh[fdr<0.05] #13
#fdr0.05
genes_up_pamh<-res_de_mbh[fdr<0.05&log2FoldChange>0.25]$external_gene_name
genes_up_pamh #"Kitl"    "Ramp3"   "Hacd2"   "C4b"     "Slc17a7" "Nsmf"    "Shox2"   "Gm43398" "Rora"    "Impa2"  

genes_up_ctrl<-res_de_mbh[fdr<0.05&log2FoldChange<(-0.25)]$external_gene_name
genes_up_ctrl
# # "Erdr1"   "Bnip5"   "Gm49749"
intersect(genes_up_pamh,genes_pamh)
intersect(genes_up_pamh,genes_pamh2)

intersect(genes_up_ctrl,genes_ctrl) #all 0 overlap

#p 0.05
genes_up_pamh2<-res_de_mbh[pvalue<0.05&log2FoldChange>0.25]$external_gene_name
genes_up_pamh2 #356

genes_up_ctrl2<-res_de_mbh[pvalue<0.05&log2FoldChange<(-0.25)]$external_gene_name
genes_up_ctrl2#300
intersect(genes_up_pamh2,genes_pamh)
intersect(genes_up_pamh2,genes_pamh2)

intersect(genes_up_ctrl2,genes_ctrl) #ns overlap

#overlap with genes differential accessibility and DEGs
library(fgsea)
genes_exclu<-union(genes_ctrl,genes_pamh2)
degs<-union(genes_up_ctrl2,genes_up_pamh2)

intersect(genes_exclu,degs)
OR(genes_exclu,degs,size_universe = length(unique(res_de_mbh$external_gene_name)))
#p=0.25



#fgsea in genes with differential accessibility (closest gene)
#stat = stat of DE test
res_de_mbh[order(stat)]
genes_stats<-res_de_mbh[order(stat)]$stat
names(genes_stats)<-res_de_mbh[order(stat)]$external_gene_name
genes_stats<-genes_stats[!duplicated(names(genes_stats))]


peaksgenes_list<-list(peak_exclu_ctrl=genes_ctrl,
                      peak_exclu_pamh=genes_pamh,
                    peak_exclu_pamh2=genes_pamh2,
                      peak_up_ctrl=genes_peak_up_ctrl,
                      peak_up_pamh=genes_peak_up_pamh)


res_mdh<-fgsea(pathways=peaksgenes_list,stats=genes_stats)
res_mdh #ns

#stat_gsea : sign(log2FC) * -log10(pvalue)
res_de_mbh[,stat_gsea:=sign(log2FoldChange)*-log10(pvalue)] #stat neg = upreg in control

genes_stats<-res_de_mbh[order(stat_gsea)]$stat_gsea
names(genes_stats)<-res_de_mbh[order(stat_gsea)]$external_gene_name
genes_stats<-genes_stats[!duplicated(names(genes_stats))]

res1_mdh<-fgsea(pathways=peaksgenes_list,stats=genes_stats)

res1_mdh  #ns

#others gsea stats ?
{
  #stat_gsea = sign(log2FC) * rank(-log10(pvalue))
res_de_mbh[,stat_gsea:=sign(log2FoldChange)*rank(-log10(pvalue))] #stat_gsea pos= upreg in pamh
res_de_mbh[order(stat_gsea)]

genes_stats<-res_de_mbh[order(stat_gsea)]$stat_gsea
names(genes_stats)<-res_de_mbh[order(stat_gsea)]$external_gene_name
genes_stats<-genes_stats[!duplicated(names(genes_stats))]

res2_mdh<-fgsea(pathways=peaksgenes_list,stats=genes_stats)
res2_mdh

#             pathway       pval      padj    log2err         ES       NES size
# 1:  peak_exclu_ctrl 0.19816514 0.3302752 0.13145761  0.1720516  1.091798  414
# 2:  peak_exclu_pamh 0.28884462 0.3610558 0.11101149  0.4751152  1.159910    7
# 3: peak_exclu_pamh2 0.05736138 0.2868069 0.26166352  0.2362922  1.296284  141
# 4:     peak_up_ctrl 0.38403042 0.3840304 0.09082414  0.2392918  1.048729   49
# 5:     peak_up_pamh 0.16806723 0.3302752 0.15524197 -0.3086734 -1.232257   31
#                                       leadingEdge
# 1:    Sash1,Klf12,Edaradd,Panx2,Nrip1,Rtn4rl1,...
# 2:                            Daam2,Arhgef28,Relt
# 3:         Akr1e1,Fam181b,Gbx2,Vim,Utp6,Ube2h,...
# 4: Fgfr2,Abhd12b,Anapc13,Afap1l2,Ntrk2,Cox6a1,...
# 5:         Panx1,Hdx,Kit,Bdnf,Tmem114,Rundc3b,...



#stat_gsea : simply rank of the gene from least to most DE (p_value)
res_de_mbh[,stat_gsea:=rank(-log10(pvalue))] #stat_gsea pos= upreg in CTRL
res_de_mbh[order(stat_gsea)]

genes_stats<-res_de_mbh[order(stat_gsea)]$stat_gsea
names(genes_stats)<-res_de_mbh[order(stat_gsea)]$external_gene_name
genes_stats<-genes_stats[!duplicated(names(genes_stats))]

res3_mdh<-fgsea(pathways=peaksgenes_list,stats=genes_stats,scoreType="pos")
res3_mdh
#             pathway       pval      padj    log2err        ES       NES size
# 1:  peak_exclu_ctrl 0.28471528 0.3558941 0.07235709 0.2739684 1.0429554  414
# 2:  peak_exclu_pamh 0.56643357 0.5664336 0.03992186 0.3824541 0.9500288    7
# 3: peak_exclu_pamh2 0.03996004 0.0999001 0.22496609 0.3347430 1.2141016  141
# 4:     peak_up_ctrl 0.03696304 0.0999001 0.23439265 0.4010134 1.3343589   49
# 5:     peak_up_pamh 0.12287712 0.2047952 0.12210792 0.3885401 1.2370675   31
#                                         leadingEdge
# 1:      Tmem163,Sash1,Klf12,Edaradd,Tia1,Fndc3b,...
# 2:                       Daam2,Arhgef28,Map3k7,Relt
# 3:      Znrf2,Panx1,Akr1e1,Fam181b,Gbx2,Col27a1,...
# 4: Fgfr2,Abhd12b,Anapc13,Afap1l2,Foxl2os,Hpcal1,...
# 5:                 Panx1,Hdx,Bdh1,Aven,Kit,Bdnf,...



}


#POA####
res_de_poa<-fread("/disks/DATATMP/PhD_AlexandrePelletier/atac/twas_res/twas_POA_untreatedSOPKF3_vs_untreatedControl_gene.csv.gz")


res_de_poa[fdr<0.05] #13
#fdr0.05
genes_up_pamh<-res_de_poa[fdr<0.05&log2FoldChange>0.25]$external_gene_name
genes_up_pamh #"Kitl"    "Ramp3"   "Hacd2"   "C4b"     "Slc17a7" "Nsmf"    "Shox2"   "Gm43398" "Rora"    "Impa2"  

genes_up_ctrl<-res_de_poa[fdr<0.05&log2FoldChange<(-0.25)]$external_gene_name
genes_up_ctrl
# # "Erdr1"   "Bnip5"   "Gm49749"
intersect(genes_up_pamh,genes_pamh)
intersect(genes_up_pamh,genes_pamh2) #Ttr

intersect(genes_up_ctrl,genes_ctrl) 

#p 0.05
genes_up_pamh2<-res_de_poa[pvalue<0.05&log2FoldChange>0.25]$external_gene_name
genes_up_pamh2 #187

genes_up_ctrl2<-res_de_poa[pvalue<0.05&log2FoldChange<(-0.25)]$external_gene_name
genes_up_ctrl2#379
intersect(genes_up_pamh2,genes_pamh)
intersect(genes_up_pamh2,genes_pamh2)# "Ttr"    "Dlx6"   "Panx1"  "Cxcl12"

intersect(genes_up_ctrl2,genes_ctrl)  #"Npas3" "Prrc1" "Mdga2" "Snx29"

#test overlap with genes differential accessibility and DEGs

genes_exclu<-union(genes_ctrl,genes_pamh2)
degs<-union(genes_up_ctrl2,genes_up_pamh2)

intersect(genes_exclu,degs)
OR(genes_exclu,degs,size_universe = length(unique(res_de_poa$external_gene_name)))
#p=0.9743904



#fgsea in genes with differential accessibility (closest gene)
#stat = stat of DE test
library(fgsea)

res_de_poa[order(stat)]
genes_stats<-res_de_poa[order(stat)]$stat
names(genes_stats)<-res_de_poa[order(stat)]$external_gene_name
genes_stats<-genes_stats[!duplicated(names(genes_stats))]


peaksgenes_list<-list(peak_exclu_ctrl=genes_ctrl,
                      peak_exclu_pamh=genes_pamh,
                    peak_exclu_pamh2=genes_pamh2,
                      peak_up_ctrl=genes_peak_up_ctrl,
                      peak_up_pamh=genes_peak_up_pamh)


res_poa<-fgsea(pathways=peaksgenes_list,stats=genes_stats)
res_poa
# 
#            pathway        pval       padj    log2err         ES       NES size
# 1:  peak_exclu_ctrl 0.004417373 0.02208687 0.40701792  0.2305340  1.307852  409
# 2:  peak_exclu_pamh 0.410488246 0.41048825 0.08431443 -0.4679698 -1.042976    7
# 3: peak_exclu_pamh2 0.077334626 0.12889104 0.28780513  0.2492896  1.251653  143
# 4:     peak_up_ctrl 0.036239041 0.09059760 0.32177592  0.3508234  1.412203   47
# 5:     peak_up_pamh 0.316622691 0.39577836 0.12384217  0.3044697  1.081552   31
#                                             leadingEdge
# 1:         Bcl6,Tbc1d2b,Gabra5,Mtss1,Tshz1,Arhgef16,...
# 2:                                 Map3k7,Relt,Arhgef28
# 3:                    Ttr,Dlx6,Panx1,Cxcl12,Dio2,Qk,...
# 4: Fgfr2,Ppp1r15b,Cmklr1,G630093K05Rik,Hpcal1,Sstr2,...
# 5:                Panx1,Cxcl12,Dio2,Xlr,Lrrc4c,Cnr1,...

#stat_gsea : sign(log2FC) * -log10(pvalue)
res_de_poa[,stat_gsea:=sign(log2FoldChange)*-log10(pvalue)] #stat neg = upreg in control

genes_stats<-res_de_poa[order(stat_gsea)]$stat_gsea
names(genes_stats)<-res_de_poa[order(stat_gsea)]$external_gene_name
genes_stats<-genes_stats[!duplicated(names(genes_stats))]

res1_poa<-fgsea(pathways=peaksgenes_list,stats=genes_stats)

res1_poa
#             pathway        pval       padj    log2err         ES       NES size
# 1:  peak_exclu_ctrl 0.010626677 0.02656669 0.38073040  0.2640787  1.260787  409
# 2:  peak_exclu_pamh 0.402097902 0.40209790 0.08359906 -0.5203262 -1.058499    7
# 3: peak_exclu_pamh2 0.008070195 0.02656669 0.38073040  0.3361760  1.422630  143
# 4:     peak_up_ctrl 0.050647325 0.08441221 0.32177592  0.3928663  1.379267   47
# 5:     peak_up_pamh 0.205128205 0.25641026 0.15524197  0.3734796  1.185355   31
#                                             leadingEdge
# 1:         Bcl6,Tbc1d2b,Gabra5,Mtss1,Tshz1,Arhgef16,...
# 2:                                 Map3k7,Relt,Arhgef28
# 3:                    Ttr,Dlx6,Panx1,Cxcl12,Dio2,Qk,...
# 4: Fgfr2,Ppp1r15b,Cmklr1,G630093K05Rik,Hpcal1,Sstr2,...
# 5:                Panx1,Cxcl12,Dio2,Xlr,Lrrc4c,Cnr1,...



#stat_gsea = sign(log2FC) * rank(-log10(pvalue))
res_de_poa[,stat_gsea:=sign(log2FoldChange)*rank(-log10(pvalue))] #stat_gsea pos= upreg in pamh
res_de_poa[order(stat_gsea)]

genes_stats<-res_de_poa[order(stat_gsea)]$stat_gsea
names(genes_stats)<-res_de_poa[order(stat_gsea)]$external_gene_name
genes_stats<-genes_stats[!duplicated(names(genes_stats))]

res2_poa<-fgsea(pathways=peaksgenes_list,stats=genes_stats)
res2_poa

#             pathway        pval       padj    log2err         ES       NES size
# 1:  peak_exclu_ctrl 0.007927863 0.03963931 0.38073040  0.1980668  1.306349  409
# 2:  peak_exclu_pamh 0.388785047 0.45517241 0.08916471 -0.4204460 -1.035575    7
# 3: peak_exclu_pamh2 0.407216495 0.45517241 0.10552094  0.1792585  1.016844  143
# 4:     peak_up_ctrl 0.055517747 0.13879437 0.32177592  0.3170527  1.437482   47
# 5:     peak_up_pamh 0.455172414 0.45517241 0.09196861  0.2501938  1.009477   31
#                                             leadingEdge
# 1:         Bcl6,Tbc1d2b,Gabra5,Mtss1,Tshz1,Arhgef16,...
# 2:                                 Map3k7,Relt,Arhgef28
# 3:                    Ttr,Dlx6,Panx1,Cxcl12,Dio2,Qk,...
# 4: Fgfr2,Ppp1r15b,Cmklr1,G630093K05Rik,Hpcal1,Sstr2,...
# 5:                Panx1,Cxcl12,Dio2,Xlr,Lrrc4c,Cnr1,...

#stat_gsea = rank(-log10(pvalue))
res_de_poa[,stat_gsea:=rank(-log10(pvalue))] #stat_gsea pos= upreg in pamh
res_de_poa[order(stat_gsea)]

genes_stats<-res_de_poa[order(stat_gsea)]$stat_gsea
names(genes_stats)<-res_de_poa[order(stat_gsea)]$external_gene_name
genes_stats<-genes_stats[!duplicated(names(genes_stats))]

res3_poa<-fgsea(pathways=peaksgenes_list,stats=genes_stats,scoreType="pos")
res3_poa#ns



#ccl : GSEA significant only on


#[To do] less stringeant peak-gene linkage: a peak is linked to every genes at (TSS +-100k region ) 



