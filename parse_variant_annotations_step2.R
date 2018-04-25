options(stringsAsFactors=FALSE)

library(dplyr)
library(data.table)

#anno_base <- read.csv("/home/murch/DATA/project/pheno_geno_correlation/geno_data/annotations/all_var_annotated.hg38_multianno_full.csv", header=TRUE)
anno_base <- data.table::fread("/home/murch/DATA/project/pheno_geno_correlation/geno_data/annotations/all_var_annotated.hg38_multianno_full.csv",header=TRUE, sep=",")

#anno_base <- as.data.frame(anno_base)
data.table::setDF(anno_base)

print(nrow(anno_base))

freq_col <- c("gnomAD_genome_ALL","1000g2015aug_all","Kaviar_AF","esp6500siv2_all","ExAC_ALL","abraom_freq","GME_AF")

print(anno_base)
#anno_base[,..freq_col] <- as.numeric(unlist(anno_base[,..freq_col]))

anno_base[,freq_col] <- apply(anno_base[,freq_col],2,function(y){
  
  return(as.numeric(y))
  
})
print(unlist(lapply(anno_base, class)))
print(nrow(anno_base))

## Calculate aggregate MAF for each site

colnames(anno_base)[which(colnames(anno_base)=="1000g2015aug_all")] <- "X1000g2015aug_all"

#gnomAD_genome_ALL,X1000g2015aug_all,Kaviar_AF,esp6500siv2_all,ExAC_ALL,HRC_AF,HRC_non1000G_AF,abroam_freq,GME_AF

anno_base <- anno_base %>%
  dplyr::rowwise() %>%
  dplyr::mutate(all_MAF_max = pmax(gnomAD_genome_ALL,X1000g2015aug_all,Kaviar_AF,esp6500siv2_all,ExAC_ALL,abraom_freq,GME_AF,na.rm=TRUE)) %>% 
  dplyr::mutate(all_MAF_min = pmin(gnomAD_genome_ALL,X1000g2015aug_all,Kaviar_AF,esp6500siv2_all,ExAC_ALL,abraom_freq,GME_AF, na.rm=TRUE)) %>%
  dplyr::mutate(all_MAF_mean = mean(c(gnomAD_genome_ALL,X1000g2015aug_all,Kaviar_AF,esp6500siv2_all,ExAC_ALL,abraom_freq,GME_AF), na.rm=TRUE)) %>%
dplyr::mutate(all_MAF_median = median(c(gnomAD_genome_ALL,X1000g2015aug_all,Kaviar_AF,esp6500siv2_all,ExAC_ALL,abraom_freq,GME_AF), na.rm=TRUE))
### START - All var

print(nrow(anno_base))
con<-file("/home/murch/DATA/project/pheno_geno_correlation/geno_data/annotations/all_var_annotated.hg38_multianno_full_aggregMAF.csv",encoding="UTF-8")
write.table(anno_base,file=con,sep=",",quote=TRUE,row.names = FALSE,col.names=TRUE)