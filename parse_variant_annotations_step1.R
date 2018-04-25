options(stringsAsFactors=FALSE)

library(dplyr)
library(data.table)


#file_path_1 <- "/home/murch/DATA/project/pheno_geno_correlation/geno_data/annotations/all_var_annotated_sample.csv"
#file_path_2 <- "/home/murch/DATA/project/pheno_geno_correlation/geno_data/annotations/all_var_annotated_sample_set2.csv"

file_path_1 <- "/home/murch/DATA/project/pheno_geno_correlation/geno_data/annotations/all_var_annotated.hg38_multianno.csv"
file_path_2 <- "/home/murch/DATA/project/pheno_geno_correlation/geno_data/annotations/all_var_annotated_set2.hg38_multianno.csv"

#anno_base_set1 <- read.table(file_path_1, sep=",", header=TRUE)
anno_base_set1 <- data.table::fread(file_path_1, header=TRUE)
anno_base_set2 <- data.table::fread(file_path_2, header=TRUE, fill=TRUE)

#anno_base_set1[,"Start"] <- as.numeric(anno_base_set1[,"Start"])
#anno_base_set1[,"End"] <- as.numeric(anno_base_set1[,"End"])

anno_base <- dplyr::left_join(anno_base_set1, anno_base_set2, by= c("Chr","Start","End","Ref","Alt"))

#anno_base <- cbind(anno_base_set1,anno_base_set2)

anno_base[anno_base=="."] <- NA

con<-file("/home/murch/DATA/project/pheno_geno_correlation/geno_data/annotations/all_var_annotated.hg38_multianno_full.csv",encoding="UTF-8")
write.table(anno_base,file=con,sep=",",quote=TRUE,row.names = FALSE,col.names=TRUE)