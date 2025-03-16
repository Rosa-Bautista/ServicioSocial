### funcion final con Mgnify 

#lo que hace la función es ordenar una lista de codigos de mgnify segun la cantidad de muestras
#que tiene el proyecto, toma cada uno de los codigo y corre lo necesario para obtener un objeto phyloseq
#con metadatos para apartir de ahí hacer redes/otros análisis

#need meter un objeto donde se tenga SOLO los codigos de Mgnify y la cantidad de muestras que tiene
mgn_phyloseq <- function (lista) {
  id_mgn <- c()
num_muestra <- c()

for(i in 1: length(solonames1)){
  id_mgn[i]<-solonames[solonames1[i],1]
  num_muestra[i]<-solonames[solonames1[i],2]  
}
 print(id_mgn)
  accession_list_ <- c()
meta_dataframe_ <- c()
mae_ <- c()
mae_phyloseq_ <- c()
abundances_ <- c()
meta_data_ <- c()
taxtable_ <- c()
otutable_ <- c()

for(i in 1:id_mgn){
  accession_list_[i]<- searchAnalysis(mgclnt, "studies",as.character([i]), usecache = TRUE) ##aqui checarlo que si funcione el as.character
  meta_dataframe_[i] <- getMetadata(mgclnt, accession_list_[i], usecache = TRUE)
  mae_[i] <- getResult(mgclnt, meta_dataframe_[i]$analysis_accession, usecache = TRUE)
  mae_[i]
  if (mae_[i] > 1){ #aqui tengo que mejorar lo de los objetos que crea,  o hago uno para cada uno de loas posibles numeros de objetos o hago uno que haga x cantidad de objetos segun se necesite
    mae_phyloseq_[i]_ <- c()
    tax_table_[i]_ <- c()
    otu_table_[i]_ <- c()
    mae[i]_phyloseq_ <- c()
    for(i in 1: length(mae_[i])){
      mae_phyloseq_[i]_[i] <- makePhyloseqFromTreeSummarizedExperiment(mae_[i][[i]])
      save(mae_phyloseq_[i]_[i],file="03_Results/Phyloseq_objects/mae_phyloseq_[i]_[i].rds")
      library(microbiome)
      View(abundances(mae_phyloseq_[i],transform = "identity")) ## checar si es que influye y c+omo influye que tenga varios
      abundances_[i]<-abundances(mae_phyloseq_[i],transform = "identity")
      View(abundances_[i])
      meta_data_3_[i]<-as.data.frame(sample_data(mae_phyloseq_[i]))
      write.csv(meta_dataframe_3_[i],file="03_Results/meta_data_[i]<.csv")
      taxtable_3_[i]<-as.data.frame(tax_table(mae_phyloseq_3_[i]))
      write.csv(taxtable_3_[i],file="03_Results/taxtable_3_[i].csv")
      otutable_3_[i]<-as.data.frame(otu_table(mae_phyloseq_3_[i]))
      write.csv(otutable_3_[i],file="03_Results/otutable_3_[i].csv")
      mae3_phyloseq_[i] <-makePhyloseqFromTreeSummarizedExperiment(mae_)
      View(sample_data(mae_phyloseq_3_[i]))
    }
  } else {
    mae_phyloseq_3 <- makePhyloseqFromTreeSummarizedExperiment(mae_3)
    save(mae_phyloseq_3,file="03_Results/Phyloseq_objects/mae_phyloseq_3.rds")
    library(microbiome)
    View(abundances(mae_phyloseq_3,transform = "identity")) ## checar si es que influye y c+omo influye que tenga varios
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
}
}


