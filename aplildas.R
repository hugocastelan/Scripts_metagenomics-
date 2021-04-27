


# AUTHOR  Hugo Castelan Sanchez 
# CREATED (2019)
# USAGE  aplildada.R archive_table.txt 
# DESCRIPTION
# this script is used to make stacked graphs 


library(ggplot2)
require(reshape2)
#Argumentos
args <- commandArgs(TRUE)
file_in <- as.character(args[1])
grafica<-read.table(file_in, sep=",", header=T, row.names=1)

#Proporcion
prop<-prop.table(data.matrix(grafica), 2)
dat_m <- melt(prop)
colnames(dat_m)<-c("Genus", "Metagenome", "Abundance")

#Grafica
 ggplot(dat_m, aes(Metagenome, Abundance, fill =Genus) + 
   geom_bar(stat = "identity")+xlab("Metagenome")+ 
   ylab("Relative abundance")+theme_bw()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+
   scale_fill_manual(values = c("darkslategrey", "darkseagreen4", "seagreen3", 
                                "darkolivegreen1","forestgreen","mediumspringgreen",
                                "limegreen","magenta3","maroon4","purple","orchid4",
                                "palevioletred4","salmon","peru","sienna1","sienna3",
                                "tan3","tomato","skyblue","turquoise","turquoise3","skyblue4"
                                ,"slateblue4","steelblue","steelblue","deepskyblue3",
                                "dodgerblue4","cornflowerblue","antiquewhite","antiquewhite4"
                                ,"azure3","azure4","gainsboro","gray","gray56","honeydew","honeydew3"
                                ,"lightcyan3","lightcyan4","lightgoldenrod4","lightgoldenrod","khaki"
                                ,"khaki4","gold3","goldenrod1","gold","lightgoldenrod4","midnightblue"
                                ,"navyblue","royalblue","skyblue1","steelblue4","lightcyan1","cadetblue"
                                ,"cadetblue3","aquamarine4","aquamarine"))
