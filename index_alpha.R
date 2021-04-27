# AUTHOR  Hugo Castelan Sanchez 
# CREATED (2019)
# USAGE  Rscript index_alpha.R archive_table.txt 
# DESCRIPTION
# this script is used to make grahps of index alpha


library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")

args <- commandArgs(TRUE)
file_in <- as.character(args[1])

datos<-read.table(file_in,sep=",",header=T, row.names=1)
OTU <- otu_table(datos, taxa_are_rows = TRUE)
GP <- prune_species(speciesSums(OTU) > 0, OTU)
#Todos los indices
plot_richness(GP)
#Solo Cho1 y Shannon 
plot_richness(GP, measures=c("Chao1", "Shannon"))

