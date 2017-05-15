#Homology Mapper

library(dplyr)

test <- read.table("Nur77GeneSpeciesComparison.csv",sep=",",quote = "",stringsAsFactors=FALSE)
testMus <- read.table("expMusMusculusFilt100000.csv",sep=",",quote = "",header=TRUE,stringsAsFactors=FALSE)
musRNAseq <- read.table("mouseRNAseq.csv",sep=",",quote = "",header=TRUE,stringsAsFactors=FALSE)


homologyMapRaw <- read.table("homologene.data",sep="\t",quote = "",stringsAsFactors=FALSE)
homologyTaxonID <- read.table("taxid_taxname.txt",sep="\t",quote = "")

colnames(homologyTaxonID) <- c("TaxonID","Species")
homologyColName <- read.table("HomoloGene_ColNames.txt",sep="\t",quote = "",stringsAsFactors=FALSE)

colnames(homologyMapRaw) <- unlist(homologyColName)

EnsemblMapRaw <- read.table("EnsembleToGeneName.txt",sep=",",quote = "",header=TRUE,stringsAsFactors=FALSE)

EnsemblMap <- unique(EnsemblMapRaw[,c(1,3)])


mapGeneNameToHID <- function(geneList,species)
{
  
  speciesTaxon <- homologyTaxonID[which(homologyTaxonID[,"Species"]==species),"TaxonID"]
  
  
  homologyMapBySpecies <- homologyMapRaw %>%
    filter(TaxonomyID==speciesTaxon)
  
  homologyID <- lapply(geneList,function(x)
    {
    
    return(homologyMapBySpecies[which(homologyMapBySpecies[,"GeneSymbol"]==x),"HomologyID"])
    
    
  })
  
  return(homologyID)
}

mapHIDToGeneName <- function(geneList,species)
{
  
  speciesTaxon <- homologyTaxonID[which(homologyTaxonID[,"Species"]==species),"TaxonID"]
  
  homologyMapBySpecies <- homologyMapRaw %>%
    filter(TaxonomyID==speciesTaxon)
  
  geneName <- lapply(geneList,function(x)
  {
    
    return(homologyMapBySpecies[which(homologyMapBySpecies[,"HomologyID"]==x),"GeneSymbol"])
    
    
  })
  
  return(geneName)
}

mapEnsemblToGeneName <- function(ensemblList,species)
{
  
 # speciesTaxon <- homologyTaxonID[which(homologyTaxonID[,"Species"]==species),"TaxonID"]
  
  #homologyMapBySpecies <- homologyMapRaw %>%
    #filter(TaxonomyID==speciesTaxon)
  
  geneName <- lapply(ensemblList,function(x)
  {
    
    return(EnsemblMap[which(EnsemblMap[,"Gene.stable.ID"]==x),"Gene.name"])
    
    
    
  })
  
  return(geneName)
}

test1 <- mapGeneNameToHID(test[-1,1],"Homo sapiens")
test2 <- mapGeneNameToHID(test[-1,2],"Mus musculus")

testMus2 <- mapEnsemblToGeneName(testMus[,1],"Homo sapiens")


testMus[,1] <- unlist(testMus2)

testRNAseq <- mapEnsemblToGeneName(musRNAseq[,1],"Homo sapiens")

#musRNAseq <- cbind(musRNAseq,testRNAseq)
testRNAseq <- as.matrix(as.character(testRNAseq))

musRNAseq <- cbind(musRNAseq,testRNAseq)

#musRNAseq[,1] <- testRNAseq'



dataSet <- cbind(test1,test2)

overlapDataSet <- dataSet[which(dataSet[,"test1"] %in% dataSet[,"test2"]),]

colnames(dataSet) <- c("Homo sapiens","Mus musculus")

overlapDataSet <- as.data.frame(dataSet[which(dataSet[,1] %in% dataSet[,2]),])

overlap <- mapHIDToGeneName(overlapDataSet[,1],"Homo sapiens")

overlap <- as.data.frame(unlist(overlap),stringsAsFactors=FALSE)

colnames(overlap) <- "GeneName"

con<-file("RNAfixedGeneNAme.csv",encoding="UTF-8")
write.table(musRNAseq,file=con,sep=",",quote=FALSE,row.names = FALSE,col.names=TRUE)

con<-file("ComparisonHomologyGene.csv",encoding="UTF-8")
write.table(dataSet,file=con,sep=",",quote=FALSE,row.names = FALSE,col.names=TRUE)

con<-file("overlap.csv",encoding="UTF-8")
write.table(overlap,file=con,sep=",",quote=FALSE,row.names = FALSE,col.names=TRUE)


con<-file("Ly6CMouseByGeneName.csv",encoding="UTF-8")
write.table(testMus,file=con,sep=",",quote=FALSE,row.names = FALSE,col.names=TRUE)

