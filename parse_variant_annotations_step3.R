options(stringsAsFactors=FALSE)

library(dplyr)
library(data.table)

anno_base <- data.table::fread("/home/murch/DATA/project/pheno_geno_correlation/geno_data/annotations/all_var_annotated.hg38_multianno_full_aggregMAF.csv", header=TRUE)

print(unlist(lapply(anno_base, class)))

step <- "START: All variants with annotations"
var_total <- nrow(anno_base)

filter_steps <- data.frame(step, var_total)

##Filter by genotype quality/read depth etc



### Select only exonic/spling var

coding_opts <- c("exonic", "splicing")
anno_coding <- anno_base %>%
  dplyr::filter(Func.refGene %in% coding_opts)

filter_steps <- rbind(filter_steps, c("Exon and splicing variants only", nrow(anno_coding)))

con<-file("/home/murch/DATA/project/pheno_geno_correlation/geno_data/annotations/all_var_annotated.hg38_multianno_exonic.csv",encoding="UTF-8")
write.table(anno_coding,file=con,sep=",",quote=TRUE,row.names = FALSE,col.names=TRUE)
### Select non-synonymous

exon_coding_opts <- c("synonymous SNV")

anno_coding_nonsynon <- anno_coding %>%
  dplyr::filter(!(ExonicFunc.refGene %in% exon_coding_opts))

filter_steps <- rbind(filter_steps, c("Non-synonymous variants", nrow(anno_coding_nonsynon)))

con<-file("/home/murch/DATA/project/pheno_geno_correlation/geno_data/annotations/all_var_annotated.hg38_multianno_exonic_nonsynon.csv",encoding="UTF-8")
write.table(anno_coding_nonsynon,file=con,sep=",",quote=TRUE,row.names = FALSE,col.names=TRUE)

### Select common alleles

anno_coding_nonsynon_common <- anno_coding_nonsynon  %>%
  dplyr::filter(all_MAF_mean > 0.1)

con<-file("/home/murch/DATA/project/pheno_geno_correlation/geno_data/annotations/all_var_annotated.hg38_multianno_common.csv",encoding="UTF-8")
write.table(anno_coding_nonsynon_common,file=con,sep=",",quote=TRUE,row.names = FALSE,col.names=TRUE)

filter_steps <- rbind(filter_steps, c("Common Alleles [>0.1 MAF] ->", nrow(anno_coding_nonsynon_common)))

### Select rare alleles

anno_coding_nonsynon_rare <- anno_coding_nonsynon  %>%
  dplyr::filter(all_MAF_mean < 0.05)

filter_steps <- rbind(filter_steps, c("Rare Alleles [<0.05 MAF] ->", nrow(anno_coding_nonsynon_rare)))

con<-file("/home/murch/DATA/project/pheno_geno_correlation/geno_data/annotations/all_var_annotated.hg38_multianno_rare.csv",encoding="UTF-8")
write.table(anno_coding_nonsynon_rare,file=con,sep=",",quote=TRUE,row.names = FALSE,col.names=TRUE)

### Select very rare alleles

anno_coding_nonsynon_vrare <- anno_coding_nonsynon  %>%
  dplyr::filter(all_MAF_mean < 0.01)

filter_steps <- rbind(filter_steps, c("Very Rare Alleles [<0.01 MAF] ->", nrow(anno_coding_nonsynon_vrare)))

con<-file("/home/murch/DATA/project/pheno_geno_correlation/geno_data/annotations/all_var_annotated.hg38_multianno_vrare.csv",encoding="UTF-8")
write.table(anno_coding_nonsynon_vrare,file=con,sep=",",quote=TRUE,row.names = FALSE,col.names=TRUE)

print(filter_steps)



con<-file("/home/murch/DATA/project/pheno_geno_correlation/geno_data/annotations/all_var_annotated.hg38_multianno_filter_steps.csv",encoding="UTF-8")
write.table(filter_steps,file=con,sep=",",quote=TRUE,row.names = FALSE,col.names=TRUE)

