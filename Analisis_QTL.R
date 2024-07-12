##Este código requerirá de los siguientes 4 paquetes
##Instalar en caso que sea necesario

#install.packages("Matrix")
#install.packages("lme4")
#install.packages("qtl")
#install.packages("qtlcharts")

##Cargar las bibliotecas
library(Matrix)
library(lme4)
library(qtl)
library(qtlcharts)

rm(list = ls()) # Limpiando el ambiente global

# Guardando la ubicación del directorio de trabajo,
# una vez que fuese definida manualmente
directorio <- getwd()
head(directorio)

## Importando datos de genotipo y fenotipo
qtlfile <- read.cross("csvs", ".", "Cortijo_2014-Epigenotypes.csv",
                      "Cortijo_2014-Phenotypes.csv", crosstype = "riself" )

##Resumen datos originales
summary.data.frame(qtlfile)

## Convertir valores extremos (+-3 Dev. Est.) en valores no disponibles (NA)

qtlfile$pheno[c(which(
  qtlfile$pheno[,2] > 
    (mean(qtlfile$pheno[,2], na.rm = TRUE) + 3*(sd(qtlfile$pheno[,2], na.rm = TRUE))))
  ,
  
  which(
    qtlfile$pheno[,2] < 
      (mean(qtlfile$pheno[,2], na.rm = TRUE) - 3*(sd(qtlfile$pheno[,2], na.rm = TRUE))))),2] <- NA

#Resumen luego de eliminar valores extremos
summary.data.frame(qtlfile)

## Cálculo de probabilidades en base al modelo hidden Markov
qtlfile <- calc.genoprob(qtlfile, step=2, error.prob = 0.001)

## Escaneo del genoma mediante regresión Haley-Knott
bloom.scan1 <- scanone(qtlfile, pheno.col = 2, method = "hk")
summary.data.frame(bloom.scan1)

## Realizando el análisis de QTL y las permutaciones para la estimación de
## umbrales significativos
bloom.scan1.perm <- scanone(qtlfile, pheno.col = 2, method = "hk", n.perm = 1000)
summary(bloom.scan1.perm, alpha=c(0.6, 0.10, 0.05, 0.01))

## Generando un gráfico interactivo de curva de valores LOD
iplotScanone(bloom.scan1)


## Función PDF para guardar gráficos con resultados del análisis
pdf(file = "lod_graph.pdf")

##Visualización
plot(bloom.scan1, col = "blue",
     xlab = "Posición en el cromosoma",
     ylab = "LOD score/Threshold",
     ylim=c(0, max(as.matrix(bloom.scan1[,-c(1:2)]))),
     main="Mapeo de QTLs asociados a tiempo de floración en Arabidopsis")

add.threshold(bloom.scan1, perms=bloom.scan1.perm, alpha = 0.05,
              lty="dashed", col = "blue")

summary(bloom.scan1, perm=bloom.scan1.perm, lodcolum=1, alpha=0.05)

dev.off()







