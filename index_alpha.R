


library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")

datos<-read.table("~/Desktop/generos_termofilos.csv",sep=",",header=T, row.names=1)
OTU <- otu_table(datos, taxa_are_rows = TRUE)
GP <- prune_species(speciesSums(OTU) > 0, OTU)
#Todos los indices
plot_richness(GP)
#Solo Cho1 y Shannon 
plot_richness(GP, measures=c("Chao1", "Shannon"))

