# Cargar las bibliotecas necesarias
library(ggplot2)
library(reshape2)

# Leer los datos
virus <- read.table("Matriz_virus.csv", header = TRUE, sep = ",", row.names = 1)
virus.t <- t(virus)  # Transponer los datos
virus.matrix <- as.matrix(virus.t)  # Convertir a matriz

# Convertir la matriz en un data frame largo para ggplot
virus_melt <- melt(virus.matrix)

# Crear una paleta de colores
my_palette <- colorRampPalette(c("white", "blue"))(100)

# Crear el heatmap con ggplot2 sin dendrograma
ggplot(data = virus_melt, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "black") +  # Crear las celdas del heatmap con borde negro
  scale_fill_gradientn(colours = my_palette) +  # Aplicar la paleta de colores
  theme_minimal() +  # Usar un tema simple
  labs(title = "Metagenomes", x = "Variables", y = "Samples") +  # Etiquetas
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  # Rotar los nombres en el eje X
