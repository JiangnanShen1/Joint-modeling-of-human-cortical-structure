trait_list = 649:715
remov=which(trait_list==682)
Trait1 = trait_list[-remov]

traits=c("ASD", "ADHD", "AD", "AN", "ADs", "BD", "epilepsy", "MDD", "OCD","SCZ", "PD", "EA", "NSM", "SmkInit", 
         "DrnkWk", "CP", "SD", "Insomnia","RA", "Crohn", "BC", "AAM", "ANM", "T2D", "Height", "BMI", "HDL", "LDL","IS", "CAD")
x <- expand.grid(Trait1,traits)


area <- as.data.frame(matrix(ncol=6,nrow=1980))
names(area) <- c("Trait1","Trait2","rg","se","z","p")
area[,1] <- x[,1]
area[,2] <- x[,2]

txtfile <- paste0('/gpfs/ysm/pi/zhao-data/zz375/UKB_analyses/area_Com30/',x[,1],'_',x[,2],'.log')

for (i in 1:1980) {
  dta <- read.delim(txtfile[i])
  dta <- strsplit(as.character(dta[54,]),split = " ")
  dta <- as.data.frame(dta)
  dta = dta[!dta[,1]=="",]
  area$rg[i] <- as.numeric(as.character(dta[3]))
  area$se[i] <- as.numeric(as.character(dta[4]))
  area$z[i] <- as.numeric(as.character(dta[5]))
  area$p[i] <- as.numeric(as.character(dta[6]))
}
write.table(area,"/gpfs/ysm/pi/zhao-data/zz375/UKB_analyses/area_Com30/areaCom30_aparc_Desikan_correlation.csv",quote=F,row.names = F,sep=",")