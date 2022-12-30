############################################################################################################################
######### Análisis de metaprofiling ###############################################################
#############################################################################################################################

######### Hugo Castelan ################################################################################################


#DADA2 para el procesamiento de lecturas
#DADA2 es un paquete de software que modela y corrige errores de amplicón secuenciados por Illumina. DADA2 infiere secuencias de muestra con exactitud, sin granularidad gruesa en Unidad Taxonómica Operativa (OTU), y resuelve diferencias de tan solo un nucleótido. En varias comunidades simuladas. DADA2 ha destacado de otras herramientas similares en identificar más variantes reales produciendo menos secuencias espurias.


#Instalación de librerias
#Antes de la instalación de DADA2 es importante contar con la libreria devtools. El objetivo de devtools es facilitar su vida como desarrollador de paquetes al proporcionar funciones de R que simplifican muchas tareas comunes.
install.packages("devtools")

#Posteriormente instalamos DADA2.
#El pipeline DADA2 tiene como punto de partida un conjunto de archivos fastq secuenciados por Illumina ("desmultiplexado"), que realiza un control de calidad que elimina quimeras y adaptadores. El producto final es una tabla de "variantes de secuencia de amplicón" (ASV), que es un análogo de mayor resolución de la tabla OTU tradicional que ofrecen otros programas.


library("devtools")
library(dada2)


#Instalamos phyloseq
#El paquete phyloseq es una herramienta para importar, almacenar, analizar y mostrar gráficamente datos de secuenciación filogenética complejos que ya se han agrupado en unidades taxonómicas operativas (OTU), especialmente cuando hay datos de muestra asociados, árbol filogenético y/o asignación taxonómica de las OTUs.

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("phyloseq")

library(phyloseq)

#Instalamos microbiome

install_github("microbiome/microbiome")

library(microbiome)


#Instalamos Fantaxtic
devtools::install_github("gmteunisse/Fantaxtic")

library(fantaxtic)

library(ggplot2)

library(dplyr)

library (vegan)



######### Read sequence data ###################################################################################


#Ruta de las lecturas
#Indicamos la ruta de las lecturas
#Cargar los archivos de la corrida de secuenciación.


PATH = "/content/"
list.files(PATH)


fnFs <- sort(list.files(PATH, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(PATH, pattern="_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#Visualizamos control de calidad de las lecturas de las lecturas Forward

plotQualityProfile(fnFs[1:3])

#Visualizamos control de calidad de las lecturas de las lecturas Reverse


plotQualityProfile(fnRs[1:3])

#Filtrado de las lecturas
#Indicamos que generamos una subcaperta o subdirectorio llamado filtered donde se van almacenar los datos


filtFs <- file.path(PATH, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(PATH, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#Filtrado de las lecturas por tamaño y con una calidad de Q20 

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(150,100),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE)

#Creating output directory: /content//filtered

#Visualización de las lecturas después del trimming
#Visualizamos el filtrado de lecturas Forward, aquí se esta visualizando unicamente 3 lecturas  

plotQualityProfile(filtFs[1:3])

#Visualizamos el filtrado de lecturas Reverse

plotQualityProfile(filtRs[1:3])


head(out)

#Filtrado de errores en las lecturas
#Para este paso se utiliza la función lernErrors que es un modelo de error paramétrico (err), donde cada conjunto de datos de lecturas tiene un conjunto diferente de tasas de error. El método learnErrors aprende este modelo de error de los datos, alternando la estimación de las tasas de error y la inferencia de la composición de la muestra hasta que convergen en una solución conjunta consistente.


errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

#Merge de las lecturas


mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)


#Remoción de chimeras

seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
write.table(seqtab.nochim, file = "AmpliconSequenceVariableTable.txt", sep = "\t")

#Resumen de los filtrados de control de calidad

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track
write.table(track,"Track_reads.tsv", sep="\t", quote=F, col.names=NA)


#Asignación taxonómica con la base de datos de SILVA
#Puedes descargar la base de datos desde la siguiente liga recuerda que tienes que descargar el archivo silva_nr_v132_train_set.fa.gz y subirlo a este colab


taxa<-assignTaxonomy(seqtab.nochim, "/content/silva_nr_v132_train_set.fa.gz", multithread=TRUE)


#Removemos las secuencias y visualizamos la tabla de asignación taxonómica

taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)
write.table(taxa, file = "TaxonomyTable.txt", sep = "\t")
saveRDS(seqtab.nochim, "seqtab_final1.rds")
saveRDS(taxa, "tax_final1.rds")

seqtab <- readRDS("seqtab_final1.rds") 
taxtab <- readRDS("tax_final1.rds")


#Leemos los metadaos
metadata<-read.table("/content/sample_data/sampledata.txt",sep="\t", header=T, row.names=1)
samples = sample_data(metadata)



#Guardamos y creamos un objeto Phyloseq

asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)


ps <- phyloseq(otu_table(seqtab_nochim, taxa_are_rows=FALSE), 
               sample_data(metadata), 
               tax_table(taxa))


#Gráficas apiladas a nivel de Phylum 
ps1_phylum <- tax_glom(ps, "Phylum", NArm = TRUE)

ps1_phylum_relabun <- transform_sample_counts(ps1_phylum, function(OTU) OTU/sum(OTU) * 100)

taxa_abundance_table_phylum <- psmelt(ps1_phylum_relabun)

StackedBarPlot_phylum <- taxa_abundance_table_phylum %>% 
  ggplot(aes(x =SampleID, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  labs(x = "Samples",
       y = "Relative Abundance",
       title = "Phylum Relative Abundance") +
  facet_grid(~ Diet, scales = "free") +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12)
  )

StackedBarPlot_phylum


#Gráficas apiladas a nivel de Genero
ps1_genus <- tax_glom(ps, "Genus", NArm = TRUE)

ps1_genus_relabun <- transform_sample_counts(ps1_phylum, function(OTU) OTU/sum(OTU) * 100)

taxa_abundance_table_genus <- psmelt(ps1_genus_relabun)

StackedBarPlot_genus <- taxa_abundance_table_genus %>% 
  ggplot(aes(x =SampleID, y = Abundance, fill = genus)) +
  geom_bar(stat = "identity") +
  labs(x = "Samples",
       y = "Relative Abundance",
       title = "Genus Relative Abundance") +
  facet_grid(~ Diet, scales = "free") +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12)
  )

StackedBarPlot_phylum



#Gráficas apiladas para un top 20 

top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, fill="Family")


#Graficamos con Fantaxtic a nivel de genero
fantaxtic_bar(ps.top20, color_by = "Family", label_by = "Genus", other_label = "Other")



#Boxplot Phylum 
BoxPlot_phylum <- taxa_abundance_table_phylum %>% 
  ggplot(aes(x =Phylum, y = Abundance, fill = Phylum)) +
  geom_boxplot() +
  labs(x = "Phylum",
       y = "Relative Abundance",
       title = "Phylum Relative Abundance") +
  facet_grid(~ Diet, scales = "free") +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12)
  )



#Indices Alfa

plot_richness(objetophy, color="Site") 

#Indices Chao1, Shannon y Simpson 

plot_richness(objetophy, measures=c("Chao1", "Shannon", "Simpson"),color="Site") + geom_point(size=5, alpha=0.7)

#Calculo de diversidad beta 

pcoa_bc = ordinate(ps, "PCoA", "bray") 

plot_ordination(ps, pcoa_bc, color = "Site") + geom_point(size=5, alpha=0.7)

#Curva de rafefaction
ggrare(ps, step = 10, label = NULL, color = NULL,
  plot = TRUE, parallel = FALSE, se = TRUE)
