###PARA QUE ORDENE LOS ARCHIVOS DE ACUERDO A QUÉ TAN PESADOS SON PARA HACER LO DEMÁS

solonames <- read.csv("pruebaservicio - Hoja 1.csv")
solonames
solonames <- as.data.frame(solonames)
solonames

id_mgn <- c()
num_muestra <- c()

for(i in 1: length(solonames1)){
  id_mgn[i]<-solonames[solonames1[i],1]
  num_muestra[i]<-solonames[solonames1[i],2]  
}



###ciclo para agilizar el mgnify 

  #de forma individual 
# aqui ponemos el mgnify que queremos para que lo busque
accession_list_3 <- searchAnalysis(mgclnt, "studies","MGYS00002309", usecache = TRUE)

# descargar el mgnify que puse arriba
meta_dataframe_3 <- getMetadata(mgclnt, accession_list_3, usecache = TRUE)

# Convert analyses outputs to a single `MultiAssayExperiment` object
#Esto también puede tardar mucho
mae_3 <- getResult(mgclnt, meta_dataframe_3$analysis_accession, usecache = TRUE)
mae_3
if (mae_3 > 1){ #aqui tengo que mejorar lo de los objetos que crea,  o hago uno para cada uno de loas posibles numeros de objetos o hago uno que haga x cantidad de objetos segun se necesite
  mae_phyloseq_3_ <- c()
  tax_table_3_ <- c()
  otu_table_3_ <- c()
  mae3_phyloseq_ <- c()
  for(i in 1: length(mae_3)){
    mae_phyloseq_3_[i] <- makePhyloseqFromTreeSummarizedExperiment(mae_3[[i]])
    save(mae_phyloseq_3_[i],file="03_Results/Phyloseq_objects/mae_phyloseq_3_[i].rds")
    library(microbiome)
    View(abundances(mae_phyloseq_3,transform = "identity")) ## checar si es que influye y c+omo influye que tenga varios
    abundances_3<-abundances(mae_phyloseq_3,transform = "identity")
    View(abundances_3)
    meta_data_1<-as.data.frame(sample_data(mae_phyloseq_1))
    write.csv(meta_dataframe_1,file="03_Results/meta_data_1.csv")
    taxtable_3_[i]<-as.data.frame(tax_table(mae_phyloseq_3_[i]))
    write.csv(taxtable_3_[i],file="03_Results/taxtable_3_[i].csv")
    otutable_3_[i]<-as.data.frame(otu_table(mae_phyloseq_3_[i]))
    write.csv(otutable_3_[i],file="03_Results/otutable_3_[i].csv")
    mae3_phyloseq_[i] <-makePhyloseqFromTreeSummarizedExperiment(mae_)
    View(sample_data(mae_phyloseq))
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

  # con varios en desorden     ### crear antes los objetos que tienen la [i] pero vacios
for(i in 1:solonames1){
  accession_list_[i]<- searchAnalysis(mgclnt, "studies",as.character([i]), usecache = TRUE) ##aqui error,checarlo
  meta_dataframe_[i] <- getMetadata(mgclnt, accession_list_[i], usecache = TRUE)
  mae_[i] <- getResult(mgclnt, meta_dataframe_[i]$analysis_accession, usecache = TRUE)
  mae_[i]
  if (mae_[i] > 1){ #aqui tengo que mejorar lo de los objetos que crea,  o hago uno para cada uno de loas posibles numeros de objetos o hago uno que haga x cantidad de objetos segun se necesite
    mae_phyloseq_[i]_ <- c()
    tax_table_[i]_ <- c()
    otu_table_[i]_ <- c()
    mae[i]_phyloseq_ <- c()
    for(i in 1: length(mae_3)){
      mae_phyloseq_3_[i] <- makePhyloseqFromTreeSummarizedExperiment(mae_3[[i]])
      save(mae_phyloseq_3_[i],file="03_Results/Phyloseq_objects/mae_phyloseq_3_[i].rds")
      library(microbiome)
      View(abundances(mae_phyloseq_3,transform = "identity")) ## checar si es que influye y c+omo influye que tenga varios
      abundances_3<-abundances(mae_phyloseq_3,transform = "identity")
      View(abundances_3)
      meta_data_1<-as.data.frame(sample_data(mae_phyloseq_1))
      write.csv(meta_dataframe_1,file="03_Results/meta_data_1.csv")
      taxtable_3_[i]<-as.data.frame(tax_table(mae_phyloseq_3_[i]))
      write.csv(taxtable_3_[i],file="03_Results/taxtable_3_[i].csv")
      otutable_3_[i]<-as.data.frame(otu_table(mae_phyloseq_3_[i]))
      write.csv(otutable_3_[i],file="03_Results/otutable_3_[i].csv")
      mae3_phyloseq_[i] <-makePhyloseqFromTreeSummarizedExperiment(mae_)
      View(sample_data(mae_phyloseq))
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

  # con varios ordenandolos+
for(i in 1: length(solonames1)){
  id_mgn[i]<-solonames[solonames1[i],1]
  num_muestra[i]<-solonames[solonames1[i],2]  
}