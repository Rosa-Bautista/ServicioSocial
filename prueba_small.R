####Prueba para ver si corre bien y qué es lo que obtengo

# es con el mgnify de menos muestras para que no sea tan pesado

# If you dont have installed this package

  #para tener la version de bioconductor que necesitan
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.20")

  #para los paquetes
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("MGnifyR")
library(MGnifyR)
library(phyloseq)
library(mia)


# Set up the MGnify client instance
mgclnt <- MgnifyClient(usecache = TRUE, cache_dir = "/tmp/MGnify_cache")

# aqui ponemos el mgnify que queremos para que lo busque
accession_list_3 <- searchAnalysis(mgclnt, "studies","MGYS00002309", usecache = TRUE)

# descargar el mgnify que puse arriba
meta_dataframe_3 <- getMetadata(mgclnt, accession_list_3, usecache = TRUE)

# Convert analyses outputs to a single `MultiAssayExperiment` object
#Esto también puede tardar mucho
mae_3 <- getResult(mgclnt, meta_dataframe_3$analysis_accession, usecache = TRUE)
mae_3
if (mae_3 > 1){  
  mae_phyloseq_3_ <- c()
  tax_table_3_ <- c()
  otu_table_3_ <- c()
  mae3_phyloseq_ <- c()
  for(i in 1: length(mae_3)){
    mae_phyloseq_3_[i] <- makePhyloseqFromTreeSummarizedExperiment(mae_3[[i]])
    save(mae_phyloseq_3_[i],file="03_Results/Phyloseq_objects/mae_phyloseq_3_[i].rds")
    library(microbiome)
    View(abundances(mae_phyloseq_3,transform = "identity")) 
    abundances_3<-abundances(mae_phyloseq_3,transform = "identity")
    View(abundances_3)
    meta_data_3_[i]<-as.data.frame(sample_data(mae_phyloseq_[i]))
    write.csv(meta_dataframe_3_[i],file="03_Results/meta_data_[i]<.csv")
    taxtable_3_[i]<-as.data.frame(tax_table(mae_phyloseq_3_[i]))
    write.csv(taxtable_3_[i],file="03_Results/taxtable_3_[i].csv")
    otutable_3_[i]<-as.data.frame(otu_table(mae_phyloseq_3_[i]))
    write.csv(otutable_3_[i],file="03_Results/otutable_3_[i].csv")
    mae3_phyloseq_[i] <-makePhyloseqFromTreeSummarizedExperiment(mae_)
    View(sample_data(mae_phyloseq_3_[i]))
} else {
  mae_phyloseq_3 <- makePhyloseqFromTreeSummarizedExperiment(mae_3)
  save(mae_phyloseq_3,file="03_Results/Phyloseq_objects/mae_phyloseq_3.rds")
  library(microbiome)
  View(abundances(mae_phyloseq_3,transform = "identity")) 
  abundances_3<-abundances(mae_phyloseq_3,transform = "identity")
  View(abundances_3)
  meta_data_3<-as.data.frame(sample_data(mae_phyloseq_3))
  write.csv(meta_dataframe_3,file="03_Results/meta_data_3.csv")
  taxtable_3<-as.data.frame(tax_table(mae_phyloseq_3))
  write.csv(taxtable_3, file = "03_Results/taxtable_3.csv")
  otutable_3<-as.data.frame(otu_table(mae_phyloseq_3))
  write.csv(otutable_3,file="03_Results/otutable_3.csv")
  mae3_phyloseq <- makePhyloseqFromTreeSummarizedExperiment(mae3)
  View(sample_data(mae_phyloseq3)) 
}
