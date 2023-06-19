隆############################################################################################################################
######### An谩lisis de metaprofiling ###############################################################
#############################################################################################################################

######### Hugo Castelan ################################################################################################


#DADA2 para el procesamiento de lecturas
#DADA2 es un paquete de software que modela y corrige errores de amplic贸n secuenciados por Illumina. DADA2 infiere secuencias de muestra con exactitud, sin granularidad gruesa en Unidad Taxon贸mica Operativa (OTU), y resuelve diferencias de tan solo un nucle贸tido. En varias comunidades simuladas. DADA2 ha destacado de otras herramientas similares en identificar m谩s variantes reales produciendo menos secuencias espurias.


#Instalaci贸n de librerias
#Antes de la instalaci贸n de DADA2 es importante contar con la libreria devtools. El objetivo de devtools es facilitar su vida como desarrollador de paquetes al proporcionar funciones de R que simplifican muchas tareas comunes.
install.packages("devtools")

#Posteriormente instalamos DADA2.
#El pipeline DADA2 tiene como punto de partida un conjunto de archivos fastq secuenciados por Illumina ("desmultiplexado"), que realiza un control de calidad que elimina quimeras y adaptadores. El producto final es una tabla de "variantes de secuencia de amplic贸n" (ASV), que es un an谩logo de mayor resoluci贸n de la tabla OTU tradicional que ofrecen otros programas.


library("devtools")
library(dada2)


#Instalamos phyloseq
#El paquete phyloseq es una herramienta para importar, almacenar, analizar y mostrar gr谩ficamente datos de secuenciaci贸n filogen茅tica complejos que ya se han agrupado en unidades taxon贸micas operativas (OTU), especialmente cuando hay datos de muestra asociados, 谩rbol filogen茅tico y/o asignaci贸n taxon贸mica de las OTUs.

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("phyloseq")

library(phyloseq)

#Instalamos microbiome

install_github("microbiome/microbiome")

library(microbiome)

devtools::install_github("microsud/microbiomeutilities")

library(microbiomeutilities)


#Instalamos Fantaxtic
devtools::install_github("gmteunisse/Fantaxtic")

library(fantaxtic)


library(ggplot2)

devtools::install_github("david-barnett/microViz")

library(microViz)

library(dplyr)

library (vegan)

library(compositions)

library (svglite)

library(patchwork)

library(ggpubr)

library(RColorBrewer)


######### Read sequence data ###################################################################################


#Ruta de las lecturas
#Indicamos la ruta de las lecturas
#Cargar los archivos de la corrida de secuenciaci贸n.


PATH = "D:/UNAM/Semestre_9/Prueba_Maleki/16S_Maleki"
list.files(PATH)


fnFs <- sort(list.files(PATH, pattern="_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(PATH, pattern="_2.fastq", full.names = TRUE))
#sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`,1)
sample.names <- sapply(strsplit(basename(fnFs), "\\."), `[`, 1)

#Visualizamos control de calidad de las lecturas de las lecturas Forward

plotQualityProfile(fnFs[1:99])

#Visualizamos control de calidad de las lecturas de las lecturas Reverse


plotQualityProfile(fnRs[1:99])

#Filtrado de las lecturas
#Indicamos que generamos una subcaperta o subdirectorio llamado filtered donde se van almacenar los datos


filtFs <- file.path(PATH, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(PATH, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#Filtrado de las lecturas por tama帽o y con una calidad de Q20 

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250,200),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE)
head(out)

#Creating output directory: /content//filtered

#Visualizacion de las lecturas despues del trimming
#Visualizamos el filtrado de lecturas Forward, aqu铆 se esta visualizando unicamente 3 lecturas  

plotQualityProfile(filtFs[1:3])

#Visualizamos el filtrado de lecturas Reverse

plotQualityProfile(filtRs[1:3])


#Filtrado de errores en las lecturas
#Para este paso se utiliza la funci贸n lernErrors que es un modelo de error param茅trico (err), donde cada conjunto de datos de lecturas tiene un conjunto diferente de tasas de error. El m茅todo learnErrors aprende este modelo de error de los datos, alternando la estimaci贸n de las tasas de error y la inferencia de la composici贸n de la muestra hasta que convergen en una soluci贸n conjunta consistente.


errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

#Denoising
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

#Merge de las lecturas


mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)


#Remocion de quimeras

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


#Asignaci贸n taxon贸mica con la base de datos de SILVA
#Puedes descargar la base de datos desde la siguiente liga recuerda que tienes que descargar el archivo silva_nr_v132_train_set.fa.gz y subirlo a este colab


taxa<-assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)


#Removemos las secuencias y visualizamos la tabla de asignaci贸n taxon贸mica

taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)
write.table(taxa, file = "TaxonomyTable.txt", sep = "\t")
saveRDS(seqtab.nochim, "seqtab_final1.rds")
saveRDS(taxa, "tax_final1.rds")

seqtab <- readRDS("seqtab_final1.rds") 
taxtab <- readRDS("tax_final1.rds")


#Leemos los metadatos
metadata<-read.table("sampledata.txt",sep="\t", header=T, row.names=1)
samples = sample_data(metadata)

#En caso de iniciar con archivos RDS, creamos un objeto Phyloseq

asv_seqs <- colnames(seqtab)
asv_headers <- vector(dim(seqtab)[2], mode="character")
asv_tab <- t(seqtab)
row.names(asv_tab) <- sub(">", "", asv_headers)
asv_tax <- taxtab
row.names(asv_tax) <- sub(">", "", asv_headers)

write.table(seqtab.nochim,file="seqtab.nochim.txt", sep="\t", quote = F)
write.table(taxa,file="taxa.txt", sep="\t", quote = F)

ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), 
               sample_data(metadata), 
               tax_table(taxtab))


#Continuando sin RDS, guardamos y creamos un objeto Phyloseq

asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)

write.table(seqtab.nochim,file="seqtab.nochim.txt", sep="\t", quote = F)
write.table(taxa,file="taxa.txt", sep="\t", quote = F)

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(metadata), 
               tax_table(taxa))

# TOP Taxa

mycols <- c("brown3", "steelblue", "grey50", "green", "red")
top_tax <- plot_taxa_boxplot(ps,
                             taxonomic.level="Phylum",
                             top.otu = 3,
                             group = "Time",
                             group.order = c("S0", "S1", "S2", "S3", "S4"),
                             group.colors = mycols,
                             add.violin = TRUE,
                             title = "Top Three Phylum",
                             keep.other = FALSE,
                             dot.size=1)

svglite("Top_phylum_abd.svg", width=15, height=10)
top_tax
dev.off()

# Revisar la funci髇 plot_listed_taxa para gr醘icas similares a la anterior


# An醠isis Beta-diversidad con distancia Aitchison entre 2 grupos

phy_sub <- subset_samples(ps, Time %in% c("S0", "S1"))
ibd <- phy_sub
ibd <- tax_fix(ibd)

plot_aitchison <- ibd %>%
  tax_transform("identity", rank = "Genus") %>% # don't transform!
  dist_calc("aitchison") %>%
  ord_calc("PCoA") %>%
  ord_plot(color = "Type", shape = "Time", plot_taxa = 1:5, size = 4) +
  scale_colour_brewer(palette = "Dark2")
svglite("Beta-Diversity_Aitchison.svg", width=10, height=10)
plot_aitchison
dev.off()

# An醠isis Beta-diversidad con distancia CLR entre 2 grupos

plot_clr <- ibd %>%
  tax_transform("clr", rank = "Genus") %>% # don't transform!
  ord_calc("PCA") %>%
  ord_plot(color = "Type", shape = "Time", plot_taxa = 1:5, size = 4) +
  scale_colour_brewer(palette = "Dark2")
svglite("Beta-Diversity_PCAclr.svg", width=10, height=10)
plot_clr
dev.off()

phy_sub2 <- subset_samples(ps, Time %in% c("S0", "S4"))
ibd2 <- phy_sub2
ibd2 <- tax_fix(ibd2)

ibd2 %>%
  tax_transform("clr", rank = "Genus") %>% # don't transform!
  ord_calc("PCA") %>%
  ord_plot(color = "Type", shape = "Time", plot_taxa = 1:10, size = 4) +
  scale_colour_brewer(palette = "Dark2")


########## CORE MICROBIOME ANALYSIS ######################

# Keep only taxa with positive sums
ps.filt <- prune_taxa(taxa_sums(ps)>0, ps)
print(ps.filt)
# Calculate compositional version of the data (relative abundances)
ps.filt.rel <- microbiome::transform(ps.filt, "compositional")
print(ps.filt.rel)
#Change it to OTUs name
dna <- Biostrings::DNAStringSetList(taxa_names(ps.filt.rel))
taxa_names(dna) <- taxa_names(ps.filt.rel)
ps.filt.rel2 <- merge_phyloseq(ps.filt.rel, dna)
taxa_names(ps.filt.rel2) <- paste0("ASV", seq(ntaxa(ps.filt.rel2)))
taxa_names(ps.filt.rel2)[1:3]

# Core microbiota analysis
core.taxa.standard <- core_members(ps.filt.rel2, detection = 0, prevalence = 5/100)
ps.core <- core(ps.filt.rel2, detection = 0, prevalence = .05)
core.taxa <- taxa(ps.core)

# Get the taxonomy data
#tax.mat <- tax_table(ps.core)
#tax.df <- as.data.frame(tax.mat)

# Add the OTUs to last column
#tax.df$OTU <- rownames(tax.df)

# Select taxonomy of only those OTUs that are core members based on the thresholds that were used
#core.

#Core visualization - genus level
ps.core.gen <- aggregate_taxa(ps.core, "Genus")
prevalences <- seq(.5, 1, .5)
detections <- 10^seq(log10(1e-4), log10(.2), length = 10)
p1 <- plot_core(ps.core.gen,
                plot.type = "heatmap",
                colours = rev(brewer.pal(5,"RdBu")),
                prevalences = prevalences,
                detections = detections, min.prevalence = 0.05) +
  xlab("Detection Threshold (Relative Abundance (%)")
p1 <- p1 + theme_bw() + ylab("ASVs")
p1


#Gr谩ficas apiladas a nivel de Genero
ps1_genus_core <- tax_glom(ps.core, "Genus", NArm = TRUE)

ps1_genus_relabun_core <- transform_sample_counts(ps1_genus_core, function(OTU) OTU/sum(OTU) * 100)

taxa_abundance_table_genus_core <- psmelt(ps1_genus_relabun_core)

write.table(taxa_abundance_table_genus_core,file="taxa_abundance_table_genus_core.txt", sep="\t", quote = F)

StackedBarPlot_genus_core <- taxa_abundance_table_genus_core %>% 
  ggplot(aes(x =Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  labs(x = "Samples",
       y = "Relative Abundance",
       title = "Genus Relative Abundance") +
  facet_grid(~ Time, scales = "free") +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12)
  )


StackedBarPlot_genus_core

#Indices Chao1, Shannon y Simpson

p.shannon <- boxplot_alpha(ps.core, index = "shannon", x_var = "Time")
p.shannon <- p.shannon + theme_minimal() + labs(x="Time", y="Shannon diversity") + 
  theme(
    axis.text = element_text(size = 12), 
    axis.title = element_text(size=16),
    legend.text=element_text(size=12),
    legend.title=element_text(size=16)
  )
p.shannon

p.simpson <- boxplot_alpha(ps.core, index = "simpson", x_var = "Time")
p.simpson <- p.simpson + theme_minimal() + labs(x="Time", y="Simpson diversity") + 
  theme(
    axis.text = element_text(size = 12), 
    axis.title = element_text(size=16),
    legend.text=element_text(size=12),
    legend.title=element_text(size=16)
  )
p.simpson

# Guardar con patchwork library
svglite("Alpha_diversity_core_microbiome.svg", width=20, height=10)
(p.shannon | p.simpson)
dev.off()

#otra representaci髇 con ggpubr library
ggarrange(p.shannon,p.simpson, legend="right")



######### END OF CORE MICROBIOME ANALYSIS ###################3








#Gr谩ficas apiladas a nivel de Phylum 
ps1_phylum <- tax_glom(ps, "Phylum", NArm = TRUE)

ps1_phylum_relabun <- transform_sample_counts(ps1_phylum, function(OTU) OTU/sum(OTU) * 100)

taxa_abundance_table_phylum <- psmelt(ps1_phylum_relabun)

write.table(taxa_abundance_table_phylum,file="taxa_abundance_table_phylum.txt", sep="\t", quote = F)

StackedBarPlot_phylum <- taxa_abundance_table_phylum %>% 
  ggplot(aes(x =Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  labs(x = "Samples",
       y = "Relative Abundance",
       title = "Phylum Relative Abundance") +
  facet_grid(~ Time, scales = "free") +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12)
  )


svglite("Phylum_socks.svg", width=20, height=10)
StackedBarPlot_phylum
dev.off()


#Gr谩ficas apiladas a nivel de Genero
ps1_genus <- tax_glom(ps, "Genus", NArm = TRUE)

ps1_genus_relabun <- transform_sample_counts(ps1_genus, function(OTU) OTU/sum(OTU) * 100)

taxa_abundance_table_genus <- psmelt(ps1_genus_relabun)

write.table(taxa_abundance_table_genus,file="taxa_abundance_table_genus.txt", sep="\t", quote = F)

StackedBarPlot_genus <- taxa_abundance_table_genus %>% 
  ggplot(aes(x =Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  labs(x = "Samples",
       y = "Relative Abundance",
       title = "Genus Relative Abundance") +
  facet_grid(~ Time, scales = "free") +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12) 
  ) + guides(fill = guide_legend(nrow = 34))


svglite("Genus_socks.svg", width=25, height=10)
StackedBarPlot_genus
dev.off()


#Gr谩ficas apiladas para un top 20 

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
  facet_grid(~ Time, scales = "free") +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12)
  )



#Indices Alfa

plot_richness(ps, color="Time") 

#Indices Chao1, Shannon y Simpson 

plot_richness(ps, measures=c("Chao1", "Shannon"),color="Time") + geom_point(size=5, alpha=0.7)
plot_richness(ps, measures=c("Chao1", "Shannon"),color="Time") + geom_boxplot(aes(fill = "Time"), alpha=0.7)

#Calculo de diversidad beta 

pcoa_bc = ordinate(objetophy, "PCoA", "bray") 

plot_ordination(objetophy, pcoa_bc, color = "Site") + geom_point(size=5, alpha=0.7)

#Curva de rafefaction
ggrare(objetophy, step = 10, label = NULL, color = NULL,
  plot = TRUE, parallel = FALSE, se = TRUE)
