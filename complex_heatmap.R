# Cargar la librería ComplexHeatmap y otras necesarias
library(ComplexHeatmap)
library(circlize)  # Para control de paletas de colores

# Leer los datos
virus <- read.table("Matriz_virus.csv", header = TRUE, sep = ",", row.names = 1)
virus.t <- t(virus)  # Transponer los datos
virus.matrix <- as.matrix(virus.t)  # Convertir a matriz

# Definir una paleta de colores para el heatmap
my_palette <- colorRamp2(c(min(virus.matrix), mean(virus.matrix), max(virus.matrix)), 
                         c("white", "blue", "darkblue"))

# Crear el heatmap con ComplexHeatmap
heatmap <- Heatmap(virus.matrix, 
                   name = "Metagenomes",  # Nombre de la leyenda
                   col = my_palette,  # Usar la paleta definida
                   row_dend_side = "left",  # Dendrograma en el lado izquierdo
                   column_dend_side = "top",  # Dendrograma en la parte superior
                   cluster_rows = TRUE,  # Clustering por filas
                   cluster_columns = TRUE,  # Clustering por columnas
                   row_names_side = "left",  # Etiquetas de las filas a la izquierda
                   column_names_rot = 90,  # Rotar nombres de las columnas
                   heatmap_legend_param = list(title = "Expression", 
                                               title_position = "topcenter", 
                                               legend_direction = "horizontal"))

# Guardar la imagen en formato PNG
png("virus_matrix_heatmap_complex.png", width = 5*600, height = 5*400, res = 300)
draw(heatmap)  # Dibujar el heatmap
dev.off()  # Cerrar el archivo PNG
