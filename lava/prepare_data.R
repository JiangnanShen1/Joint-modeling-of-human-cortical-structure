##########LAVA
# 1. Partition files: dir of SUPERGNOVA after reformation /gpfs/gibbs/pi/zhao/yz738/Partition/block/chr
# change the format to the formot of LAVA, sample files of partition in LAVA:/gpfs/gibbs/pi/zhao/cz354/new_local_gc_compare/partition/snp_ldsplit/1KG_eur/lava
# path:/gpfs/gibbs/pi/zhao/yz738/Partition/block/chr@.txt
#
# 2. Reference panel: just put plink file (for different chromsomes)
# path: /gpfs/gibbs/pi/zhao/yz738/1000G/eur_chr@_SNPmaf5
# 
# 3. Summary statistics file:  sample file /gpfs/gibbs/pi/zhao/cz354/new_local_gc_compare/real/code/lava_info.txt 
# 
# 4. sample.overlap file: sample file /gpfs/gibbs/pi/zhao/cz354/new_local_gc_compare/real/code/sample.overlap.txt 
# Code to generate it : /gpfs/gibbs/pi/zhao/cz354/new_local_gc_compare/real/code/overlap.R 
# github code: https://github.com/josefin-werme/LAVA/blob/main/vignettes/sample_overlap.Rmd
# 
# 
# DSQ file to run : /gpfs/gibbs/pi/zhao/cz354/new_local_gc_compare/real/code/lava.sh



###########
python3 /ysm-gpfs/pi/zhao/yz738/software/SUPERGNOVA/supergnova.py /gpfs/ysm/pi/zhao-data/zz375/UKB_analyses/Munged_data/area_aparc_Desikan/649munge.sumstats.gz /gpfs/ysm/pi/zhao/yz738/Munged/AAM.sumstats.gz  --bfile /ysm-gpfs/pi/zhao/yz738/1000G/eur_chr@_SNPmaf5 --partition /ysm-gpfs/pi/zhao/yz738/Partition/block/chr@.txt --out ~/scratch60/IMAGE_SUPERGNOVA/AAM/AAM_area649_supergnova.txt --thread 8;
#########restructured the partition file:
lava_path <- "/gpfs/gibbs/pi/zhao/js4377/UKB/lava/"
LOC<-1
 for(chr in 1:22){
   l1 <- read.table(paste0("/gpfs/gibbs/pi/zhao/yz738/Partition/block/chr",chr,".txt"), header = T)
   colnames(l1)<- c("CHR", "START", "STOP")
   for(i in 1:nrow(l1)){
     l1$LOC[i]<- LOC
     LOC<- LOC+1
   }
   print(paste0("The rows of the chr", chr, ": ", nrow(l1)))
   print(paste0("loc:",LOC))
   file.name <- paste0(lava_path, "Data/partition/partition_chr",chr,".txt")
   write.table(l1,file.name,quote = F, col.names = T, row.names = F)
 }
####################Prepare the sample overlap file for three measures

##########LAVA
# 1. Partition files: dir of SUPERGNOVA after reformation /gpfs/gibbs/pi/zhao/yz738/Partition/block/chr
# change the format to the formot of LAVA, sample files of partition in LAVA:/gpfs/gibbs/pi/zhao/cz354/new_local_gc_compare/partition/snp_ldsplit/1KG_eur/lava
# path:/gpfs/gibbs/pi/zhao/yz738/Partition/block/chr@.txt
#
# 2. Reference panel: just put plink file (for different chromsomes)
# path: /gpfs/gibbs/pi/zhao/yz738/1000G/eur_chr@_SNPmaf5
# 
# 3. Summary statistics file:  sample file /gpfs/gibbs/pi/zhao/cz354/new_local_gc_compare/real/code/lava_info.txt 
# 
# 4. sample.overlap file: sample file /gpfs/gibbs/pi/zhao/cz354/new_local_gc_compare/real/code/sample.overlap.txt 
# Code to generate it : /gpfs/gibbs/pi/zhao/cz354/new_local_gc_compare/real/code/overlap.R 
# github code: https://github.com/josefin-werme/LAVA/blob/main/vignettes/sample_overlap.Rmd
# 
# 
# DSQ file to run : /gpfs/gibbs/pi/zhao/cz354/new_local_gc_compare/real/code/lava.sh




#########restructured the partition file:
lava_path <- "/gpfs/gibbs/pi/zhao/js4377/UKB/lava/"
LOC<-1
for(chr in 1:22){
  l1 <- read.table(paste0("/gpfs/gibbs/pi/zhao/yz738/Partition/block/chr",chr,".txt"), header = T)
  colnames(l1)<- c("CHR", "START", "STOP")
  for(i in 1:nrow(l1)){
    l1$LOC[i]<- LOC
    LOC<- LOC+1
  }
  print(paste0("The rows of the chr", chr, ": ", nrow(l1)))
  print(paste0("loc:",LOC))
  file.name <- paste0(lava_path, "Data/partition/partition_chr",chr,".txt")
  write.table(l1,file.name,quote = F, col.names = T, row.names = F)
}





####################Prepare the sample overlap file for three measures

#####1. run the sh fule to generate the 
#file /gpfs/gibbs/pi/zhao/js4377/UKB/lava/Data/sample_overlap/all.sh
#generate the all_area.sh, all_volume.sh, all_thickness.sh


#####2. generate the covariance matrix
setwd("/gpfs/gibbs/pi/zhao/js4377/UKB/lava/Data/sample_overlap/")




#######Area
scor = read.table("all_area.rg",header=T)              # read in
scor = scor[,c("p1","p2","gcov_int")]             # retain key headers
scor$p1 = gsub("munge.sumstats.gz","",scor$p1)   # assuming the munged files have format [phenotype]_munge.sumstats.gz
scor$p2 = gsub(".sumstats.gz","",scor$p2)   # (adapt as necessary)
scor$p1= gsub("/gpfs/ysm/pi/zhao-data/zz375/UKB_analyses/Munged_data/area_aparc_Desikan//","",scor$p1)
scor$p2= gsub("/gpfs/ysm/pi/zhao-data/zz375/UKB_analyses/Com_Compatible//","",scor$p2)

scor1=scor
scor1$p1=scor$p2
scor1$p2=scor$p1
scor=rbind(scor,scor1)
phen = unique(scor$p1) 
n = length(phen)
mat = matrix(NA,n,n)                    # create matrix
rownames(mat) = colnames(mat) = phen    # set col/rownames
for (i in phen) {
  for (j in phen) {
    if (i==j){
      mat[i,j]=1
    }else{
      if (length(subset(scor, p1==i & p2==j)$gcov_int)==0){
        mat[i,j]= NA
      }else{
        mat[i,j] = subset(scor, p1==i & p2==j)$gcov_int
      } }}
}
#if (!all(t(mat)==mat)) { mat[lower.tri(mat)] = t(mat)[lower.tri(mat)] }  # sometimes there might be small differences in gcov_int depending on which phenotype was analysed as the outcome / predictor
mat = round(cov2cor(mat),5)  
#mat <- mat[1:66, 67:96]
# standardise
write.table(mat, "/gpfs/gibbs/pi/zhao/js4377/UKB/lava/Data/sample_overlap/sample.overlap.area.txt", quote=F)   # save


#######Volume
scor = read.table("all_volume.rg",header=T)              # read in
scor = scor[,c("p1","p2","gcov_int")]             # retain key headers
scor$p1 = gsub("munge.sumstats.gz","",scor$p1)   # assuming the munged files have format [phenotype]_munge.sumstats.gz
scor$p2 = gsub(".sumstats.gz","",scor$p2)   # (adapt as necessary)
scor$p1= gsub("/gpfs/ysm/pi/zhao-data/zz375/UKB_analyses/Munged_data/volume_aparc_Desikan//","",scor$p1)
scor$p2= gsub("/gpfs/ysm/pi/zhao-data/zz375/UKB_analyses/Com_Compatible//","",scor$p2)

scor1=scor
scor1$p1=scor$p2
scor1$p2=scor$p1
scor=rbind(scor,scor1)
phen = unique(scor$p1) 
n = length(phen)
mat = matrix(NA,n,n)                    # create matrix
rownames(mat) = colnames(mat) = phen    # set col/rownames
for (i in phen) {
  for (j in phen) {
    if (i==j){
      mat[i,j]=1
    }else{
      if (length(subset(scor, p1==i & p2==j)$gcov_int)==0){
        mat[i,j]= NA
      }else{
        mat[i,j] = subset(scor, p1==i & p2==j)$gcov_int
      } }}
}
#if (!all(t(mat)==mat)) { mat[lower.tri(mat)] = t(mat)[lower.tri(mat)] }  # sometimes there might be small differences in gcov_int depending on which phenotype was analysed as the outcome / predictor
mat = round(cov2cor(mat),5)  
#mat <- mat[1:66, 67:96]
# standardise
write.table(mat, "/gpfs/gibbs/pi/zhao/js4377/UKB/lava/Data/sample_overlap/sample.overlap.volume.txt", quote=F)   # save



#######Thickness
scor = read.table("all_thickness.rg",header=T)              # read in
scor = scor[,c("p1","p2","gcov_int")]             # retain key headers
scor$p1 = gsub("munge.sumstats.gz","",scor$p1)   # assuming the munged files have format [phenotype]_munge.sumstats.gz
scor$p2 = gsub(".sumstats.gz","",scor$p2)   # (adapt as necessary)
scor$p1= gsub("/gpfs/gibbs/pi/zhao/js4377/UKB/Munge/Munge_data/thickness-aparc-Desikan//","",scor$p1)
scor$p2= gsub("/gpfs/ysm/pi/zhao-data/zz375/UKB_analyses/Com_Compatible//","",scor$p2)

scor1=scor
scor1$p1=scor$p2
scor1$p2=scor$p1
scor=rbind(scor,scor1)
phen = unique(scor$p1) 
n = length(phen)
mat = matrix(NA,n,n)                    # create matrix
rownames(mat) = colnames(mat) = phen    # set col/rownames
for (i in phen) {
  for (j in phen) {
    if (i==j){
      mat[i,j]=1
    }else{
      if (length(subset(scor, p1==i & p2==j)$gcov_int)==0){
        mat[i,j]= NA
      }else{
        mat[i,j] = subset(scor, p1==i & p2==j)$gcov_int
      } }}
}
#if (!all(t(mat)==mat)) { mat[lower.tri(mat)] = t(mat)[lower.tri(mat)] }  # sometimes there might be small differences in gcov_int depending on which phenotype was analysed as the outcome / predictor
mat = round(cov2cor(mat),5)  
#mat <- mat[1:66, 67:96]
# standardise
write.table(mat, "/gpfs/gibbs/pi/zhao/js4377/UKB/lava/Data/sample_overlap/sample.overlap.thickness.txt", quote=F)   # save
############# Prepare info files
library(data.table)
library(readxl)
library(dplyr)
Common_traits_path <- "/gpfs/gibbs/pi/zhao/zhao-data/zz375/UKB_analyses/Com_Compatible/"
Common_info_path <- "/gpfs/gibbs/pi/zhao/js4377/UKB/lava/Data/summary_statistics/pheno_info.xlsx"
common_info <- as.data.frame(read_excel(Common_info_path))
common_info<- common_info[,c("Abbreviation","case","control" )]
colnames(common_info)<- c("phenotype", "cases","controls")
common_info$prevalence <- NA
common_info <- common_info %>%
  rowwise() %>%
  mutate(filename = paste0(Common_traits_path, phenotype, ".sumstats.gz"))


generate_INFO<- function(measure,ROI_measure_path){
file_names<- list.files(ROI_measure_path)
phenotype<- gsub("munge.log","", file_names)
phenotype<- gsub("munge.sumstats.gz","", phenotype)
phenotype<- unique(phenotype)
ROI<- data.frame(phenotype=phenotype, cases=NA, controls= NA, prevalence= NA)
ROI$filename<- paste0(paste0(ROI_measure_path, phenotype, "munge.sumstats.gz"))
return(ROI)
}

ROI_measure_path <- "/gpfs/gibbs/pi/zhao/zhao-data/zz375/UKB_analyses/Munged_data/area_aparc_Desikan/"
ROI_area<- generate_INFO("area",ROI_measure_path)

ROI_measure_path <- "/gpfs/gibbs/pi/zhao/zhao-data/zz375/UKB_analyses/Munged_data/volume_aparc_Desikan/"
ROI_vol<- generate_INFO("volume",ROI_measure_path)
ROI_measure_path <-"/gpfs/gibbs/pi/zhao/js4377/UKB/Munge/Munge_data/thickness-aparc-Desikan/"
ROI_thickness<- generate_INFO("thickness",ROI_measure_path)

ROI <- rbind(common_info, ROI_area,ROI_thickness,ROI_vol)
write.table(ROI,"/gpfs/gibbs/pi/zhao/js4377/UKB/lava/Data/summary_statistics/lava_info.txt", col.names = T, row.names = F, quote = F)