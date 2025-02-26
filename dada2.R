## DADA" original que cualquier persona podría usar solo ocupa meter sus datos 
#y mover algunos parámetros


###paquete

install.packages("dada2")
library(dada2)

## direccion donde esta el archivo con los datos
path <- "carpeta" # aqui ponemos la carpeta donde estan los fastq ya unzipped
list.files(path)


##obtener los fastq que son forward y los que son reverse
forward_fn <- sort(list.files(path, pattern="fastq con los forwards", full.names = TRUE))
reverse_fn <- sort(list.files(path, pattern="fastq con los reversa", full.names = TRUE))

## Extraer nombres de las muestras
#esto es asuminedo que tiene el formato SAMPLENAME_XXX.fastq; pero este podría variar entre archivos
sample.names <- sapply(strsplit(basename(forward_fn), "_"), `[`, 1)

## ver calidad de las secuencias
plotQualityProfile(forward_fn[1:2])
plotQualityProfile(reverse_fn[1:2])

# filtrar y corte

## crear el espacio para guardar los datos filtrados
filtFs <- file.path(path, "filtrado", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtrado", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

## guardar las secuencias cortadas e cada uno de los 
out <- filterAndTrim(forward_fn, filtFs, reverse_fn, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # en windows este ultimo es false
head(out) #no moví nada de los valores que venian en la pagina

# Error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ = TRUE)

# Sample inference a las muestras ya filtradas y cortados
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]
dadaRs[[1]]

#fusionar las lecturas F y R alineadas
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
head(mergers[[1]])

#Construct sequence table// ASV table // amplicon sequence variant table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

table(nchar(getSequences(seqtab))) #Inspeccionar la distribución de longitudes de secuencia

# quitar quimeras: easier porque ya se quito cierto ruido
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)
    ##la frecuencias de secuencias quimericas va a variar mucho entre conjuntos de datos


#Número de leccturas hechas para cada paso
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# en caso de que sea un SINGLE SAMPLE?????? quitar el sapply y poner solo getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#Asignar taxonomia
taxa <- assignTaxonomy(seqtab.nochim, ".fa.gz", multithread=TRUE)
  #revisar las asignaciones taxonomicas que se hicieron
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


#Evaluar occurracy: en el ejemplo usan mock porque una de las muestras incluidas
#en los datos que usaron fue una “comunidad simulada”, en la que se secuenció una mezcla de 20 cepas 
#conocidas. Las secuencias de referencia correspondientes 
#a estas cepas se proporcionaron en el archivo zip descargado. Volvemos a esa 
#muestra y comparamos las variantes de secuencia inferidas por DADA2 con la 
#composición esperada de la comunidad.

unqs.mock <- seqtab.nochim["Mock",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")

mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")


###PHYLOSEQ

library(phyloseq)
library(Biostrings)
library(ggplot2)
theme_set(theme_bw()) ##???

#Construir un data frame a partir de los datos del archivo
samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out

#Ahora si ya podemos hacerlo un phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps)

#cambiar por nombres más ccortos pero quedarnos con la secuencias
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

#Gráficos
  #diversidad alfa
plot_richness(ps, x="Day", measures=c("Shannon", "Simpson"), color="When")

  #Bray-Curtis
#transformar a proporciones
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
#gráfica
plot_ordination(ps.prop, ord.nmds.bray, color="When", title="Bray NMDS")

  #Barplot para ver abundancias
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Day", fill="Family") + facet_wrap(~When, scales="free_x")

