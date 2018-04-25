library(gdsfmt)
library(SNPRelate)
library(dplyr)
library(data.table)

anno_base_common <- data.table::fread("/home/murch/DATA/project/pheno_geno_correlation/geno_data/annotations/all_var_annotated.hg38_multianno_common.csv", header=TRUE)

anno_base_common_regions <- anno_base_common[,c("Chr","Start")]

con<-file("/home/murch/DATA/project/pheno_geno_correlation/geno_data/annotations/all_var_annotated.hg38_multianno_common_regions.tsv",encoding="UTF-8")
write.table(anno_base_common_regions,file=con,sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)

anno_base_exonic_nonsynon <- data.table::fread("/home/murch/DATA/project/pheno_geno_correlation/geno_data/annotations/all_var_annotated.hg38_multianno_exonic_nonsynon.csv", header=TRUE)

anno_base_exonic_nonsynon_regions <- anno_base_exonic_nonsynon[,c("Chr","Start")]

con<-file("/home/murch/DATA/project/pheno_geno_correlation/geno_data/annotations/all_var_annotated.hg38_multianno_exonic_nonsynon_regions.tsv",encoding="UTF-8")
write.table(anno_base_exonic_nonsynon,file=con,sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)

##COMMON

common_variant_vcf <- "/home/murch/DATA/project/pheno_geno_correlation/geno_data/full_set_common.vcf"

common_gds_path <- snpgdsVCF2GDS(common_variant_vcf, "common.gds", method="biallelic.only")


common_gds <- snpgdsOpen(common_gds_path)

#common_ldpruned <- snpgdsLDpruning(common_gds , ld.threshold=0.2)
#common_ldpruned_id <- unlist(common_ldpruned)

# <- snpgdsPCA(common_gds, snp.id=common_ldpruned_id, num.thread=5)
common_pca <- snpgdsPCA(common_gds, num.thread=5)

common_pca_eigen <- data.frame(sample.id = common_pca$sample.id,
                  EV1 = common_pca$eigenvect[,1],    # the first eigenvector
                  EV2 = common_pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)

plot(common_pca_eigen$EV2, common_pca_eigen$EV1, xlab="eigenvector 2", ylab="eigenvector 1")

con<-file("/home/murch/DATA/project/pheno_geno_correlation/output/pca_common_R.csv",encoding="UTF-8")
write.table(common_pca_eigen,file=con,sep="\t",quote=FALSE,row.names = FALSE,col.names=TRUE)

##ALL PCA

exonic_variant_vcf <- "/home/murch/DATA/project/pheno_geno_correlation/geno_data/full_set_exonic.vcf"

exonic_gds_path <- snpgdsVCF2GDS(exonic_variant_vcf, "common.gds", method="biallelic.only")


common_gds <- snpgdsOpen(exonic_gds_path)

#common_ldpruned <- snpgdsLDpruning(common_gds , ld.threshold=0.2)
#common_ldpruned_id <- unlist(common_ldpruned)

# <- snpgdsPCA(common_gds, snp.id=common_ldpruned_id, num.thread=5)
exonic_pca <- snpgdsPCA(exonic_gds, num.thread=5)

exonic_pca_eigen <- data.frame(sample.id = exonic_pca$sample.id,
                               EV1 = exonic_pca$eigenvect[,1],    # the first eigenvector
                               EV2 = exonic_pca$eigenvect[,2],    # the second eigenvector
                               stringsAsFactors = FALSE)

plot(exonic_pca_eigen$EV2, exonic_pca_eigen$EV1, xlab="eigenvector 2", ylab="eigenvector 1")

con<-file("/home/murch/DATA/project/pheno_geno_correlation/output/pca_exonic_R.csv",encoding="UTF-8")
write.table(exonic_pca_eigen,file=con,sep="\t",quote=FALSE,row.names = FALSE,col.names=TRUE)
