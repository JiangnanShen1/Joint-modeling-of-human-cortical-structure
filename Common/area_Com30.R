trait_list = 649:715
remov=which(trait_list==682)
Trait1 = trait_list[-remov]
trait_list2 <-c("ASD", "ADHD", "AD", "AN", "ADs", "BD", "epilepsy", "MDD", "OCD","SCZ", "PD", "EA", "NSM", "SmkInit", 
                "DrnkWk", "CP", "SD", "Insomnia","RA", "Crohn", "BC", "AAM", "ANM", "T2D", "Height", "BMI", "HDL", "LDL","IS", "CAD")

gwas_path = "/gpfs/ysm/pi/zhao-data/zz375/UKB_analyses/Munged_data/area_aparc_Desikan/"
gwas_path2 <- "/gpfs/ysm/pi/zhao-data/zz375/UKB_analyses/Com_Compatible/"
# You can use the ldscore calculated by me
ldscore = "/gpfs/ysm/pi/zhao-data/yz738/Local-Trans/HYPERGNOVA/Reference_Eur_hapmap3/chr@hapmap3"

# LDSC software
ldsc_software = "/ysm-gpfs/pi/zhao-data/zz375/IMAGE/ldsc/ldsc/ldsc.py"

# Generate ldsc for each pair
joblist = c()
for(i in 1:length(Trait1)){
    trait1 = Trait1[i]
  for(j in 1:length(trait_list2)){
    trait2 = trait_list2[j]
    # please change the python here to your python
    joblist = c(joblist, paste0("python2 ", ldsc_software, " --rg ", gwas_path, "/", trait1,
                                "munge.sumstats.gz,", gwas_path2, "/", trait2,
                                ".sumstats.gz --ref-ld-chr ", ldscore, " --w-ld-chr ",
                                ldscore, " --out ./", trait1, "_", trait2))
  }
}

write.table(joblist, "./ldsc.sh", quote=F, row.names = F, col.names=F)

# change the partition if necessary
partition = "scavenge"
# generate joblist to submit the jobs
joblist = c("module load GCC;", "partition=scavenge;", "module load dSQ",
            paste0("/gpfs/loomis/apps/avx/software/dSQ/1.05/dSQ --jobfile ldsc.sh -p ",
                   partition, " -n 1 --mem-per-cpu=10g", " -t 24:00:00", " --mail-type=ALL",
                   " --batch-file ldsc.pbs"),
            "sbatch ldsc.pbs")
write.table(joblist, "./runCom_ldsc.sh", quote=F, row.names = F, col.names=F)