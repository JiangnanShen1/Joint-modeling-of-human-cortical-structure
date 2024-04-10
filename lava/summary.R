library(data.table)







############Summary lava results
library(data.table)

measure <-"area"
pheno_pair<- read.table(paste0("/gpfs/gibbs/pi/zhao/js4377/UKB/lava/Scipt/joblist_", measure,".txt"))
pheno_pair<- unique(pheno_pair[,c("V3","V4")])
  for( i in 1740:nrow(pheno_pair)){
    print(i)
    pair=paste0(pheno_pair$V3[i],"_",pheno_pair$V4[i])
    d=fread(paste0("/gpfs/gibbs/pi/zhao/js4377/UKB/lava/Results/", measure,"/uni_",pair,"_chr1.txt"))
    e=fread(paste0("/gpfs/gibbs/pi/zhao/js4377/UKB/lava/Results/", measure,"/bi_",pair,"_chr1.txt"))
    partition=fread(paste0("/gpfs/gibbs/pi/zhao/js4377/UKB/lava/Data/partition/partition_chr1.txt"))
    d=merge(d,partition,by="LOC")
    e=merge(e,partition,by="LOC") 
    for(j in 2:22){
      d1=fread(paste0("/gpfs/gibbs/pi/zhao/js4377/UKB/lava/Results/", measure,"/uni_",pair,"_chr",j,".txt"))
      e1=fread(paste0("/gpfs/gibbs/pi/zhao/js4377/UKB/lava/Results/", measure,"/bi_",pair,"_chr",j,".txt"))
      partition=fread(paste0("/gpfs/gibbs/pi/zhao/js4377/UKB/lava/Data/partition/partition_chr",j,".txt"))
      d1=merge(d1,partition,by="LOC")
      e1=merge(e1,partition,by="LOC")
      d=rbind(d,d1)
      e=rbind(e,e1)
    }
    write.table(d,paste0("/gpfs/gibbs/pi/zhao/js4377/UKB/lava/Results/", measure,"/summary/uni_",pair,".txt"),row.names=F,col.names=T,quote=F)
    write.table(e,paste0("/gpfs/gibbs/pi/zhao/js4377/UKB/lava/Results/", measure,"/summary/bi_",pair,".txt"),row.names=F,col.names=T,quote=F)
  }

######################################









measure <-"volume"
library(data.table)
dir<-paste0("/gpfs/gibbs/pi/zhao/js4377/UKB/lava/Results/", measure)
pheno_pair<- read.table(paste0("/gpfs/gibbs/pi/zhao/js4377/UKB/lava/Scipt/joblist_", measure,".txt"))
pheno_pair<- unique(pheno_pair[,c("V3","V4")])
print(nrow(pheno_pair))
nblocks=c()
data=as.data.frame(matrix(,nrow=0,ncol=13))
colnames(data)=c("CHR","START","STOP","phen1","phen2","rho","rho.lower","rho.upper","p" ,"h2_1","p_h2_1" ,"h2_2" ,"p_h2_2")
for(i in 1:nrow(pheno_pair)){
  print(i)
  phen1=pheno_pair$V3[i]
  phen2=pheno_pair$V4[i]
  bi=fread(paste0(dir,"/summary/bi_",phen1,"_",phen2,".txt"))[,c("CHR","START","STOP","phen1","phen2","rho","rho.lower","rho.upper","p")]
  bi$phen2<- as.character(bi$phen2)
  uni1=fread(paste0(dir,"/summary/uni_",phen1,"_",phen2,".txt"))[,c("CHR","START","STOP","phen","h2.obs","p")]
  colnames(uni1)[c(4,5,6)]=c("phen1","h2_1","p_h2_1")
  uni2=fread(paste0(dir,"/summary/uni_",phen1,"_",phen2,".txt"))[,c("CHR","START","STOP","phen","h2.obs","p")]
  colnames(uni2)[c(4,5,6)]=c("phen2","h2_2","p_h2_2")
  bi=merge(bi,uni1,by=c("CHR","START","STOP","phen1"))
  bi=merge(bi,uni2,by=c("CHR","START","STOP","phen2"))
  bi=bi[!is.na(bi$rho)]
  data=rbind(data,bi)
  nblocks=c(nblocks,nrow(bi))
}

nrow(data)
nrow(data[data$h2_1>0 & data$h2_2>0,])
mu=data[data$h2_1>0 & data$h2_2>0,]
nrow(mu[abs(mu$rho)<=1.25,])
lava=data
fwrite(lava,paste0("/gpfs/gibbs/pi/zhao/js4377/UKB/lava/Results/summary/summary_",measure,"_lava.txt"),row.names = F, col.names = T,quote = F,na="NA", sep=" ")
