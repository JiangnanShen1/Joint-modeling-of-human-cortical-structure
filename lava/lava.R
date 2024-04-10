library(LAVA)
library(data.table)
library(plyr)

args = commandArgs(trailingOnly=TRUE)
pheno1= args[1]
pheno2= args[2]
#pair=args[2] 
chr=args[3]
measure=args[4]
#measure=args[4]
#output='/gpfs/gibbs/pi/zhao/cz354/new_local_gc_compare/real/result/1KG_eur/lava/'
#output="/gpfs/gibbs/pi/zhao/js4377/UKB/lava/Results/"
#pair=as.numeric(pair)
#phenotype=fread(pheno,header=F)
#pheno1 = phenotype$V1[pair]
#pheno2 = phenotype$V2[pair]
# pheno1="AAM"
# pheno2=1023
# chr=2
# measure="thickness"
#mainDir = '/gpfs/gibbs/pi/zhao/cz354/new_local_gc_compare/real/code/'
mainDir='/gpfs/gibbs/pi/zhao/js4377/UKB/lava/Data/summary_statistics/'
output=paste0("/gpfs/gibbs/pi/zhao/js4377/UKB/lava/Results/", measure, "/")
sampleOverlap="/gpfs/gibbs/pi/zhao/cz354/new_local_gc_compare/real/code/sample.overlap.txt"
sampleOverlap=paste0("/gpfs/gibbs/pi/zhao/js4377/UKB/lava/Data/sample_overlap/sample.overlap.",measure, ".txt")
#REF = paste0('/gpfs/gibbs/pi/zhao/cz354/zhao-data/cz354/genotype/1000G_phase3/EUR/ref_maf1/genetic_map_CHR/1000G_EUR_maf1_chr',chr)
REF = paste0("/gpfs/gibbs/pi/zhao/js4377/UKB/lava/Data/reference_panel/eur_chr", chr, "_SNPmaf5_rmdu")

### Read in summary statistics and related info
input = process.input(input.info.file=paste0(mainDir,"/lava_info.txt"),
                      sample.overlap.file=sampleOverlap,
                      ref.prefix=REF,
                      phenos=c(pheno1,pheno2))   

loci = read.loci(paste0("/gpfs/gibbs/pi/zhao/js4377/UKB/lava/Data/partition/partition_chr",chr,".txt"))



uni=data.frame(matrix(ncol=5,nrow=0))
colnames(uni)=c("phen", "h2.obs", "h2.latent","p","LOC")

# negative variance estimate of phenotype
# when genetic correlation larger or smaller than 1.25, they are excluded and the genetic correlation with -1.25--1 and 1-1.25 are set to be 1 or -1
for(i in 1:nrow(loci)){
  locus = process.locus(loci[i,], input)
  
  tryCatch({
    a=run.univ(locus)
    a$LOC=loci[i,]$LOC
    uni=rbind.fill(uni,a)
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


bi=data.frame(matrix(ncol=10,nrow=0))
colnames(bi)=c("phen1" ,"phen2" ,"rho" ,"rho.lower" ,"rho.upper" ,"r2" ,"r2.lower" ,"r2.upper" ,"p","LOC")


for(i in 1:nrow(loci)){
  locus = process.locus(loci[i,], input)
  
  tryCatch({
    b=run.bivar(locus)
    b$LOC=loci[i,]$LOC
    bi=rbind.fill(bi,b)
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}


write.table(uni,paste0(output,"/uni_",pheno1,"_",pheno2,"_chr",chr,".txt"),row.names=F,col.names = T,quote=F)

write.table(bi,paste0(output,"/bi_",pheno1,"_",pheno2,"_chr",chr,".txt"),row.names=F,col.names = T,quote=F)
