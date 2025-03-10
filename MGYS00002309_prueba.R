####prueba: https://shiny-portal.embl.de/shinyapps/app/06_mgnify-notebook-lab?jlpath=mgnify-examples/R%20Examples/Fetch%20Analyses%20metadata%20for%20a%20Study.ipynb&jlvar_MGYS=MGYS00002309

source("lib/variable_utils.r")
mgnify_study_accesion <- get_variable_from_link_or_input("MGYS", "Study Accession", "MGYS00002309")

library(vegan)
library(ggplot2)
library(phyloseq)

library(MGnifyR)

analyses_accessions <- mgnify_analyses_from_studies

metadata <- mgnify_get_analyses_metadata (mg, )
