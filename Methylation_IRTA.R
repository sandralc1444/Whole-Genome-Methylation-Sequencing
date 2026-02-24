### Global methylation paper - IRTA ######

### BARPLOT DE DMRs HIPER / HIPO / TOTAL ###

library(ggplot2)
library(tidyverse)
install.packages("readxl", type = "source") #no me iba la función readxl pq no encontraba funciones básicas del sistema (que raro pq antes no me ha pasado esto)
library(readxl)

# Leer datos desde Excel
data <- read_excel("DMR_Methylation_IRTA.xlsx", sheet = 1)

# Reorganizar datos a formato largo para ggplot
data_long <- data %>%
  pivot_longer(cols = c(`Hyper DMRs`, `Hypo DMRs`, `All DMRs`),
               names_to = "DMR_type",
               values_to = "Count")

# Orden de las comparaciones
data_long$Comparison <- factor(data_long$Comparison,
                               levels = c("INF_T vs CTRL_T",
                                          "INF_O vs CTRL_O",
                                          "CTRL_O vs CTRL_T",
                                          "INF_O vs INF_T"))

# Etiquetas más bonitas
data_long$DMR_type <- recode(data_long$DMR_type,
                             "Hyper DMRs" = "Hypermethylated",
                             "Hypo DMRs" = "Hypomethylated",
                             "All DMRs" = "Total")

# Crear el barplot (DMRs)
ggplot(data_long, aes(x = Comparison, y = Count, fill = DMR_type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_text(aes(label = Count), 
            vjust = -0.3,
            position = position_dodge(width = 0.8),
            size = 3) +
  scale_fill_brewer(palette = "Paired", name = "DMRs category") +
  labs(x = "Treatment and sex comparisons", y = "Differentially Methylated Regions (DMRs)", 
       title = "Differentially Methylated Regions (DMRs)") +
  theme_minimal(base_size = 12)

# Aplicar etiquetas, colores y tema
ggplot(data_long, aes(x = Comparison, y = Count, fill = DMR_type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_text(aes(label = Count),
            vjust = -0.3,
            position = position_dodge(width = 0.8),
            size = 3) +
  labs(x = "Comparison of gonadal tissues treatments", 
       y = "Differentially Methylated Regions (DMRs)") +
  scale_fill_manual(values = c(
    "Hypermethylated" = "#CABCEF", 
    "Hypomethylated" = "#A3CBA3", 
    "Total" = "#D3D3D3"),
    name = "DMRs category") +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggplot(data_long, aes(x = Comparison, y = Count, fill = DMR_type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
  geom_text(aes(label = Count),
            vjust = -0.3,
            position = position_dodge(width = 0.9),
            size = 3) +
  labs(x = "Comparison", 
       y = "Differentially Methylated Regions (DMRs)") +
  scale_fill_manual(values = c(
    "Hypermethylated" = "#CABCEF", 
    "Hypomethylated" = "#A3CBA3", 
    "Total" = "#D3D3D3"),
    name = "DMRs category") +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10)
  )

##All CpGs
# Leer datos desde Excel
data <- read_excel("DMR_Methylation_IRTA.xlsx", sheet = 2)

# Reorganizar datos a formato largo para ggplot
data_long <- data %>%
  pivot_longer(cols = c(`Hyper CpGs`, `Hypo CpGs`, `All CpGs`),
               names_to = "CpGs_type",
               values_to = "Count")

# Orden de las comparaciones
data_long$Comparison <- factor(data_long$Comparison,
                               levels = c("INF_O vs CTRL_O",
                                          "INF_T vs CTRL_T",
                                          "INF_O vs INF_T",
                                          "CTRL_O vs CTRL_T"))

# Etiquetas más bonitas
data_long$CpG_type <- recode(data_long$CpGs_type,
                             "Hyper CpGs" = "Hypermethylated",
                             "Hypo CpGs" = "Hypomethylated",
                             "All CpGs" = "Total")

# Crear el barplot
ggplot(data_long, aes(x = Comparison, y = Count, fill = CpG_type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_text(aes(label = Count), 
            vjust = -0.3,
            position = position_dodge(width = 0.8),
            size = 3) +
  scale_fill_brewer(palette = "Paired", name = "CpGs category") +
  labs(x = "Treatment and sex comparisons", y = "Differentially CpG sites", 
       title = "Differentially CpG sites") +
  theme_minimal(base_size = 12)

# Aplicar etiquetas, colores y tema
ggplot(data_long, aes(x = Comparison, y = Count, fill = CpG_type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_text(aes(label = Count),
            vjust = -0.3,
            position = position_dodge(width = 0.8),
            size = 3) +
  labs(x = "Comparison of gonadal tissues treatments", 
       y = "Differentially CpG sites") +
  scale_fill_manual(values = c(
    "Hypermethylated" = "#548B54", 
    "Hypomethylated" = "#8B4789", 
    "Total" = "#EE9A00"),
    name = "CpGs category") +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggplot(data_long, aes(x = Comparison, y = Count, fill = CpG_type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
  geom_text(aes(label = Count),
            vjust = -0.3,
            position = position_dodge(width = 0.9),
            size = 3) +
  labs(x = "Comparison", 
       y = "Differentially CpG sites") +
  scale_fill_manual(values = c(
    "Hypermethylated" = "#548B54", 
    "Hypomethylated" = "#8B4789", 
    "Total" = "#EE9A00"),
    name = "CpGs category") +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10)
  )

ggplot(data_long, aes(x = Comparison, y = Count, fill = CpG_type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
  geom_text(aes(label = Count),
            vjust = -0.3,
            position = position_dodge(width = 0.9),
            size = 3) +
  labs(x = "Comparison", 
       y = "Differentially CpG sites") +
  scale_y_continuous(labels = scales::label_comma()) +
  scale_fill_manual(values = c(
    "Hypermethylated" = "#548B54", 
    "Hypomethylated" = "#8B4789", 
    "Total" = "#EE9A00"),
    name = "CpGs category") +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10)
  )


###Separated graphs###
### --------- Gráfico 1: Comparaciones infección DMRs ---------
data_infection <- data_long %>%
  filter(Comparison %in% c("INF_T vs CTRL_T", "INF_O vs CTRL_O"))

ggplot(data_infection, aes(x = Comparison, y = Count, fill = DMR_type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
  geom_text(aes(label = Count),
            vjust = -0.3,
            position = position_dodge(width = 0.9),
            size = 3) +
  labs(x = "Comparison between treatments", 
       y = "Differentially Methylated Regions (DMRs)") +
  scale_fill_manual(values = c(
    "Hypermethylated" ="#548B54", 
    "Hypomethylated" = "#8B4789", 
    "Total" = "#EE9A00"),
    name = "DMRs category") +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10)
  )

### --------- Gráfico 2: Comparaciones de sexo (DMRs) ---------
data_sex <- data_long %>%
  filter(Comparison %in% c("CTRL_F vs CTRL_M", "INF_F vs INF_M"))

ggplot(data_sex, aes(x = Comparison, y = Count, fill = DMR_type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
  geom_text(aes(label = Count),
            vjust = -0.3,
            position = position_dodge(width = 0.9),
            size = 3) +
  labs(x = "Comparison between ovaries and testes", 
       y = "Differentially Methylated Regions (DMRs)") +
  scale_fill_manual(values = c(
    "Hypermethylated" = "#87B6D9", 
    "Hypomethylated" = "#D09EBF", 
    "Total" = "#C6C6C6"),
    name = "DMRs category") +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10)
  )


########Differentially methylated regions distribution (DMRs)#########
####only the columns of INF_O vs CTRL_O and INF_T vs CTRL_T####
# Cargar librerías
library(ggplot2)
library(dplyr)

# Crear el dataframe con las dos últimas comparaciones del gráfico
df <- data.frame(
  comparison = rep(c("INF_T vs. CTRL_T", "INF_O vs. CTRL_O"), each = 4),
  region = rep(c("promoter", "intron", "intergenic", "exon"), times = 2),
  percentage = c(30.59, 22.9, 23.8, 22.7,
                 15.42, 41.45, 26.28, 16.85)
)
df$comparison <- factor(df$comparison, levels = c("INF_T vs. CTRL_T", "INF_O vs. CTRL_O"))
df

# Crear el gráfico de barras apilado con etiquetas
ggplot(df, aes(x = comparison, y = percentage, fill = region)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  geom_text(aes(label = paste0(round(percentage, 2), "%")),
            position = position_fill(vjust = 0.5),
            color = "white", size = 4.2) +
  labs(y = "%", x = NULL, fill = "Region") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0)) +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "lightgray", color = NA), # Fondo de la zona de dibujo
  )


########DMRs and DMLs overlapping genes#########
####INF_O vs CTRL_O and INF_T vs CTRL_T####
# Leer datos desde Excel
data2 <- read_excel("DMRs and DMLs overlapping genes_barplot.xlsx", sheet = 1)

# Reorganizar datos a formato largo para ggplot
data_long2 <- data2 %>%
  pivot_longer(cols = c(`Hypermethylated`, `Hypomethylated`, `Total`),
               names_to = "gene_type",
               values_to = "Count")

# Orden de las comparaciones
data_long2$Comparison <- factor(data_long2$Comparison,
                               levels = c("INF_O vs CTRL_O",
                                          "INF_T vs CTRL_T"
                                          ))

ggplot(data_long2, aes(x = Comparison, y = Count, fill = gene_type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
  geom_text(aes(label = Count),
            vjust = -0.3,
            position = position_dodge(width = 0.9),
            size = 3) +
  labs(x = "Comparison between treatments", 
       y = "Differentially Methylated Genes (DMGs)") +
  scale_fill_manual(values = c(
    "Hypermethylated" ="#548B54", 
    "Hypomethylated" = "#8B4789", 
    "Total" = "#EE9A00"),
    name = "Genes category") +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10)
  )

##############Gene annotation - biomart#############
# Instalar solo si no lo tienes
if (!requireNamespace("biomaRt", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("biomaRt")
}

library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "dlabrax_gene_ensembl")  # European seabass

#Subo el txt (archivo en formato txt)
gene_ids <- read.table("INF_F_vs_CTRL_F_genes_DMRS.txt", header = FALSE)$V1
head(gene_ids)  # Verifica que se vean tipo ENSDLAG0000...

#Gene annotation
annotated_genes <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "description"),
  filters = "ensembl_gene_id",
  values = gene_ids,
  mart = ensembl
)

#Exportar el CSV con los reusltados
write.csv(annotated_genes, "Annotated_DMRS_genes.csv", row.names = FALSE)

##No se me cargaba la libreria de biomaRt asi que he cerrado R y he hecho:
install.packages("tidyselect")
library(biomaRt)

#European sea bass no esta disponible en la versión principal de Ensambl
library(biomaRt)

ensembl_metazoa <- useMart("metazoa_mart", host = "https://metazoa.ensembl.org")

listDatasets(ensembl_metazoa) #no esta E sea bass

ensembl_metazoa <- useMart(
  biomart = "metazoa_mart",
  dataset = "dlabrax_eg_gene",
  host = "https://metazoa.ensembl.org"
)


#ALternativa --> 2021 del ensamblaje genómico de Dicentrarchus labrax, con anotaciones de Ensembl (https://projects.ensembl.org/aqua-faang/)

install.packages("BiocManager")
BiocManager::install("rtracklayer")

library(rtracklayer)

# Abrir gtf
gtf_path <- "Dicentrarchus_labrax.dlabrax2021.106.gtf.gz"
gtf_data <- import(gtf_path)

#Subo el txt con mis genes (archivo en formato txt)
gene_ids <- read.table("INF_F_vs_CTRL_F_genes_DMRS.txt", header = FALSE)$V1
head(gene_ids)  # Verifica que se vean tipo ENSDLAG0000...

#Anotar genes
gene_annotations <- gtf_data[gtf_data$type == "gene"]

# Asegúrate de que las IDs coinciden con el campo correspondiente en el GTF
annotated_genes <- gene_annotations[gene_annotations$gene_id %in% gene_list$gene_id]

# Visualizar
head(annotated_genes)




###Definitivo --> Funciona!!
# Instala si no lo tienes
install.packages("readr")
install.packages("dplyr")
install.packages("openxlsx")
install.packages("writexl")

library(readr)
library(dplyr)
library(openxlsx)
library(writexl)
library(openxlsx)

# Cargar lista de genes
gene_list <- read_tsv("INF_F_vs_INF_M_genes_DMRS.txt", col_names = "gene_id")
gene_list <- read_excel("M_bien exportado.xlsx", sheet = 4) #Hypomethylated (significantlly) genes _ t
gene_list <- read_excel("M_bien exportado.xlsx", sheet = 6) #Hypermethylated (significantlly) genes _ t
gene_list <- read_excel("F_bienexportado.xlsx", sheet = 4) #Hypomethylated (significantlly) genes _ o
gene_list <- read_excel("F_bienexportado.xlsx", sheet = 6) #Hypermethylated (significantlly) genes _ o

gene_list <- read_excel("Dataset 5.1_DEGs.xlsx", sheet = 2) #down T
gene_list <- read_excel("Dataset 5.1_DEGs.xlsx", sheet = 3) #up T
gene_list <- read_excel("Dataset 5.1_DEGs.xlsx", sheet = 5) #down O
gene_list <- read_excel("Dataset 5.1_DEGs.xlsx", sheet = 6) #up O



# Leer el archivo GTF
gtf_lines <- readLines(gzfile("Dicentrarchus_labrax.dlabrax2021.114.gtf.gz"))

# Filtrar solo líneas que contienen "gene"
gtf_genes <- gtf_lines[grepl("\tgene\t", gtf_lines)]

# Función robusta para parsear atributos
parse_attributes <- function(attr_string) {
  fields <- strsplit(attr_string, ";\\s*")[[1]]
  attrs <- list()
  for (field in fields) {
    parts <- strsplit(field, " +")[[1]]
    if (length(parts) >= 2) {
      key <- parts[1]
      value <- gsub('"', '', parts[2])
      attrs[[key]] <- value
    }
  }
  return(attrs)
}

# Parsear el GTF a tabla con manejo de errores
gtf_parsed <- lapply(gtf_genes, function(line) {
  parts <- strsplit(line, "\t")[[1]]
  if (length(parts) < 9) return(NULL)
  attrs <- tryCatch(parse_attributes(parts[9]), error = function(e) return(NULL))
  
  if (!is.null(attrs[["gene_id"]])) {
    return(data.frame(
      gene_id = attrs[["gene_id"]],
      gene_name = ifelse(!is.null(attrs[["gene_name"]]), attrs[["gene_name"]], NA),
      gene_biotype = ifelse(!is.null(attrs[["gene_biotype"]]), attrs[["gene_biotype"]], NA),
      chromosome = parts[1],
      start = as.numeric(parts[4]),
      end = as.numeric(parts[5]),
      strand = parts[7],
      stringsAsFactors = FALSE
    ))
  } else {
    return(NULL)
  }
})

# Filtrar resultados válidos
gtf_df <- do.call(rbind, Filter(Negate(is.null), gtf_parsed))

# Verificar cuántos genes hay en el GTF
cat("Genes únicos en el GTF:", length(unique(gtf_df$gene_id)), "\n")

# Hacer join con tu lista de genes
annotated_genes <- left_join(gene_list, gtf_df, by = "gene_id")

# Mostrar cuántos genes se han podido anotar
cat("Genes anotados:", sum(!is.na(annotated_genes$gene_name.y)), "de", nrow(annotated_genes), "\n")

# Guardar resultado
write_xlsx(annotated_genes, "DEGs_up_O_v114.xlsx")

###Diagrama de sectores que compara el tipo de genes encontrados en los DMRs
# Instalar y cargar paquetes si es necesario
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
if (!require("cowplot")) install.packages("cowplot")

library(ggplot2)
library(dplyr)
library(cowplot)
library(patchwork)

##INF_FvsCTROL_F
# ------------------------------
# DATOS desde la figura
# ------------------------------

# Total de genes (todos)
all_genes <- data.frame(
  biotype = c("protein_coding", "lncRNA", "miRNA", "pseudogene", "snoRNA", "snRNA"),
  count = c(6164, 118, 8, 11, 3, 2)
)

# Total de genes anotados
annotated_genes <- data.frame(
  biotype = c("protein_coding", "lncRNA", "miRNA", "pseudogene", "snoRNA", "snRNA"),
  count = c(3532, 0, 0, 0, 1, 2)
)

# ------------------------------
# COLORES pastel personalizados
# ------------------------------

# Paleta de colores pastel
pastel_colors <- c(
  "protein_coding" = "#f4a6b5", # rosa pastel
  "lncRNA"         = "#d9b3ff", # lila pastel
  "miRNA"          = "#b3d9ff", # azul pastel
  "pseudogene"     = "#b3f0b3", # verde pastel
  "snoRNA"         = "#ffffb3", # amarillo pastel
  "snRNA"          = "#b3ffff"  # azul claro pastel
)

# Función de gráfico sin % (sin geom_text)
plot_pie <- function(data, title) {
  ggplot(data, aes(x = "", y = count, fill = biotype)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = pastel_colors) +
    labs(title = title, fill = "Biotype") +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 11)
    )
}

# Crear gráficos
p1 <- plot_pie(all_genes, "All DMRs-associated genes")
p2 <- plot_pie(annotated_genes, "Annotated DMRs-associated genes")

# Mostrar juntos con patchwork
p1 + p2

all_genes %>%
  mutate(percentage = round(count / sum(count) * 100, 2)) %>%
  arrange(desc(percentage))

annotated_genes %>%
  mutate(percentage = round(count / sum(count) * 100, 2)) %>%
  arrange(desc(percentage))

#INF_MvsCTRL_M
# Total de genes (todos)
all_genes <- data.frame(
  biotype = c("protein_coding", "lncRNA", "miRNA", "pseudogene", "snoRNA", "snRNA", "rRNA", "misc_RNA", "IG_V_gene"),
  count = c(3293, 57, 4, 10, 2, 9, 4, 2, 2)
)

# Total de genes anotados
annotated_genes <- data.frame(
  biotype = c("protein_coding", "lncRNA", "miRNA", "pseudogene", "snoRNA", "snRNA", "rRNA", "misc_RNA", "IG_V_gene"),
  count = c(1720, 0, 0, 0, 0, 9, 4, 2, 0)
)

# ------------------------------
# COLORES pastel personalizados
# ------------------------------

# Paleta de colores pastel extendida
pastel_colors <- c(
  "protein_coding" = "#f4a6b5", # rosa pastel
  "lncRNA"         = "#d9b3ff", # lila pastel
  "miRNA"          = "#b3d9ff", # azul pastel
  "pseudogene"     = "#b3f0b3", # verde pastel
  "snoRNA"         = "#ffffb3", # amarillo pastel
  "snRNA"          = "#b3ffff", # azul claro pastel
  "rRNA"           = "#f9d1ff", # rosa-lila pastel
  "misc_RNA"       = "#c9ffe5", # verde-menta pastel
  "IG_V_gene"      = "#ffdfbf"  # naranja suave pastel
)

# Función de gráfico sin % (sin geom_text)
plot_pie <- function(data, title) {
  ggplot(data, aes(x = "", y = count, fill = biotype)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = pastel_colors) +
    labs(title = title, fill = "Biotype") +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 11)
    )
}

# Crear gráficos
p1 <- plot_pie(all_genes, "All DMRs-associated genes")
p2 <- plot_pie(annotated_genes, "Annotated DMRs-associated genes")

# Mostrar juntos con patchwork
p1 + p2

# Tablas resumen con porcentajes
all_summary <- all_genes %>%
  mutate(percentage = round(count / sum(count) * 100, 2)) %>%
  arrange(desc(percentage))

annotated_summary <- annotated_genes %>%
  mutate(percentage = round(count / sum(count) * 100, 2)) %>%
  arrange(desc(percentage))

# Mostrar tablas en consola
print(all_summary)
print(annotated_summary)

#CTRL_FvsCTRL_M
# Total de genes (todos)
all_genes <- data.frame(
  biotype = c("protein_coding", "lncRNA", "miRNA", "pseudogene", "snoRNA", "snRNA", "scaRNA", "rRNA", "misc_RNA", "vault_RNA", "Y_RNA", "IG_V_gene", "ribozyme"),
  count = c(17897, 533, 44, 39, 21, 21, 1, 4, 2, 3, 1, 7, 1)
)

# Total de genes anotados
annotated_genes <- data.frame(
  biotype = c("protein_coding", "lncRNA", "miRNA", "pseudogene", "snoRNA", "snRNA", "scaRNA", "rRNA", "misc_RNA", "vault_RNA", "Y_RNA", "IG_V_gene", "ribozyme"),
  count = c(9686, 0, 6, 0, 12, 20, 0, 4, 2, 3, 1, 0, 1)
)

# ------------------------------
# COLORES pastel personalizados
# ------------------------------

# Paleta de colores pastel extendida
pastel_colors <- c(
  "protein_coding" = "#f4a6b5", # rosa pastel
  "lncRNA"         = "#d9b3ff", # lila pastel
  "miRNA"          = "#b3d9ff", # azul pastel
  "pseudogene"     = "#b3f0b3", # verde pastel
  "snoRNA"         = "#ffffb3", # amarillo pastel
  "snRNA"          = "#b3ffff", # azul claro pastel
  "rRNA"           = "#f9d1ff", # rosa-lila pastel
  "misc_RNA"       = "#c9ffe5", # verde-menta pastel
  "IG_V_gene"      = "#ffdfbf"  # naranja suave pastel
)

# Función de gráfico sin % (sin geom_text)
plot_pie <- function(data, title) {
  ggplot(data, aes(x = "", y = count, fill = biotype)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = pastel_colors) +
    labs(title = title, fill = "Biotype") +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 11)
    )
}

# Crear gráficos
p1 <- plot_pie(all_genes, "All DMRs-associated genes")
p2 <- plot_pie(annotated_genes, "Annotated DMRs-associated genes")

# Mostrar juntos con patchwork
p1 + p2

# Tablas resumen con porcentajes
all_summary <- all_genes %>%
  mutate(percentage = round(count / sum(count) * 100, 2)) %>%
  arrange(desc(percentage))

annotated_summary <- annotated_genes %>%
  mutate(percentage = round(count / sum(count) * 100, 2)) %>%
  arrange(desc(percentage))

# Mostrar tablas en consola
print(all_summary)
print(annotated_summary)


#INF_FvsINF_M
# Total de genes (todos)
all_genes <- data.frame(
  biotype = c("protein_coding", "lncRNA", "miRNA", "sca_RNA", "rRNA", "misc_RNA", "pseudogene", "snoRNA", "snRNA", "Y_RNA", "IG_V_gene", "ribozyme", "vault_RNA"),
  count = c(17984, 545, 49, 1, 3, 3, 37, 20, 23, 1, 7, 1, 2)
)

# Total de genes anotados
annotated_genes <- data.frame(
  biotype = c("protein_coding", "lncRNA", "miRNA", "sca_RNA", "rRNA", "misc_RNA", "pseudogene", "snoRNA", "snRNA", "Y_RNA", "IG_V_gene", "ribozyme", "vault_RNA"),
  count = c(9696, 0, 8, 0, 3, 3, 0, 11, 22, 1, 0, 1, 2)
)

# ------------------------------
# COLORES pastel personalizados
# ------------------------------

# Paleta de colores pastel extendida
pastel_colors <- c(
  "protein_coding" = "#f4a6b5", # rosa pastel
  "lncRNA"         = "#d9b3ff", # lila pastel
  "miRNA"          = "#b3d9ff", # azul pastel
  "pseudogene"     = "#b3f0b3", # verde pastel
  "snoRNA"         = "#ffffb3", # amarillo pastel
  "snRNA"          = "#b3ffff", # azul claro pastel
  "rRNA"           = "#f9d1ff", # rosa-lila pastel
  "misc_RNA"       = "#c9ffe5", # verde-menta pastel
  "IG_V_gene"      = "#ffdfbf"  # naranja suave pastel
)

# Función de gráfico sin % (sin geom_text)
plot_pie <- function(data, title) {
  ggplot(data, aes(x = "", y = count, fill = biotype)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = pastel_colors) +
    labs(title = title, fill = "Biotype") +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 11)
    )
}

# Crear gráficos
p1 <- plot_pie(all_genes, "All DMRs-associated genes")
p2 <- plot_pie(annotated_genes, "Annotated DMRs-associated genes")

# Mostrar juntos con patchwork
p1 + p2

# Tablas resumen con porcentajes
all_summary <- all_genes %>%
  mutate(percentage = round(count / sum(count) * 100, 2)) %>%
  arrange(desc(percentage))

annotated_summary <- annotated_genes %>%
  mutate(percentage = round(count / sum(count) * 100, 2)) %>%
  arrange(desc(percentage))

# Mostrar tablas en consola
print(all_summary)
print(annotated_summary)

###Normalización annotated genes (para que sean todos en mayúsculas)
library(readxl)
library(dplyr)
library(stringr)
library(writexl)

# Leer archivo anotado
annot <- read_xlsx("annotated_genes_with_RStudio_INF_FvsCTRL_F.xlsx")

# Ver qué valores raros hay en gene_name
unique(annot$gene_name)

# Filtrar si hay algún gene_name raro como "chr"
annot_clean <- annot %>% filter(gene_name != "chr")

# Normalizar nombres de genes (opcional):
# Por ejemplo, poner todos los nombres en mayúsculas
annot_clean <- annot_clean %>%
  mutate(gene_name = ifelse(!is.na(gene_name), toupper(gene_name), NA))

# Guardar archivo limpio
write_xlsx(annot_clean, "annotated_genes_cleaned.xlsx")


##Finding the differential methylation interval in specific interesting genes
# Leer el archivo correctamente (reemplaza el nombre si es diferente)
data <- read.delim("DSS_diff_CpG_INF_F_vs_CTRL_F_with_overlap_genes copia.txt", header = TRUE, stringsAsFactors = FALSE)

# Verifica los nombres de las columnas
print(colnames(data))

# Opcional: si hay muchos espacios como separadores
# data <- read.table("archivo.txt", header = TRUE, sep = "", stringsAsFactors = FALSE)

# Verifica el tipo de dato de la columna diff
str(data$diff)

# Filtrar por FDR y pval
filtered <- subset(data, fdr <= 0.001 & pval <= 0.001)

# Filtrar por el gene_id de interés
vash2 <- subset(filtered, gene_id == "ENSDLAG00005014412")

# Asegúrate que 'diff' sea numérica
vash2$diff <- as.numeric(vash2$diff)

# Verifica los valores únicos
summary(vash2$diff)

# Número de CpGs y rango
cat("Número de CpGs:", length(vash2$diff), "\n")
cat("Rango de diff:", min(vash2$diff), "a", max(vash2$diff), "\n")



library(readr)

# Leer bien el archivo separando por tabulación
# Esto evita que se lea como una sola columna
correct_data <- read_delim("DSS_diff_CpG_INF_F_vs_CTRL_F_with_overlap_genes copia.txt", delim = " ")
#DSS_diff_CpG_INF_F_vs_CTRL_F_with_overlap_genes copia.txt - F
#DSS_diff_CpG_INF_M_vs_CTRL_M_with_overlap_genes.txt - M


# Verifica estructura y tipos
str(correct_data$diff)  # Debe mostrar numérico

# Si diff no es numérica (por ejemplo, "character"), conviértela:
correct_data$diff <- as.numeric(correct_data$diff)

# Repite los pasos de filtrado
filtered <- subset(correct_data, fdr <= 0.001 & pval <= 0.001)
vash2 <- subset(filtered, gene_id == "ENSDLAG00005025730")

# Resultado final
cat("N CpGs:", nrow(vash2), "\n")
cat("Rango diff:", min(vash2$diff), "a", max(vash2$diff), "\n")

# Hipometilados (diff < 0)
hipo <- subset(vash2, diff < 0)
cat("\n--- HIPOMETILADOS ---\n")
cat("N CpGs:", nrow(hipo), "\n")
cat("Rango de diff:", min(hipo$diff), "a", max(hipo$diff), "\n")

# Hipermetilados (diff > 0)
hiper <- subset(vash2, diff > 0)
cat("\n--- HIPERMETILADOS ---\n")
cat("N CpGs:", nrow(hiper), "\n")
cat("Rango de diff:", min(hiper$diff), "a", max(hiper$diff), "\n")

#####FUNCTIONAL ANALYSIS#####
#############Bubble plot 15 most enriched GOs per category####################

# Cargar las librerías necesarias
library(ggplot2)
library(dplyr)
library(forcats)
library(magrittr)

##############PLOT ALL###########
# Leer el archivo CSV o Excel
df <- read.xlsx("GOterms_dlabrax_INFFvsCTRLF_allgenes.xlsx", sheet = 2)


# Agrupar por 'source' y obtener los 15 términos más enriquecidos para cada grupo
df_filtered <- df %>%
  group_by(source) %>%
  top_n(10, wt = -adjusted_p_value) %>%
  ungroup()

# Ordenar y reordenar 'term_name' por 'source'
df_filtered %<>%
  group_by(source) %>%
  arrange(term_name, .by_group = TRUE) %>%
  ungroup() %>%
  mutate(term_name = forcats::fct_reorder(term_name, source))

# Escala de colores para el gráfico
color_scale <- scale_color_gradient(low = "lightblue", high = "black")

# Crear el gráfico de burbujas
bubble_plot1 <- ggplot(df_filtered, aes(x=intersection_size, y=term_name, size = intersection_size)) + 
  geom_point(aes(color = adjusted_p_value), alpha = 0.8) + 
  geom_tile(aes(width = Inf, y=term_name, fill=source), alpha = 0.2) + 
  scale_fill_manual(values = c("#104E8B", "#8B814C", "#328BA2")) +
  color_scale +
  xlab("Gene number") +
  ylab("Top 15 enriched GO terms per category") +
  theme_bw() +  
  theme(panel.background = element_blank(), 
        axis.text = element_text(colour = "black", size = rel(1.5)),
        axis.line = element_line(size = 0.8), 
        axis.title = element_text(size = 12), 
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 10), 
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10))

# Mostrar el gráfico
print(bubble_plot1)


####GOs con barplot######
# Cargar librerías
library(ggplot2)
library(dplyr)
library(readr)

# Leer el CSV
go_df <- read_csv("gProfiler_dlabrax_19-9-2025_12-33-03__intersectionsOvaries.csv")

# Ordenar los términos por número de genes (intersection_size)
go_df <- go_df %>%
  arrange(intersection_size) %>%
  mutate(term_name = factor(term_name, levels = term_name))

# Crear gráfico
ggplot(go_df, aes(x = intersection_size, y = term_name, fill = adjusted_p_value)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "red3", high = "blue3") +
  labs(
    x = "Gene Count",
    y = NULL,
    fill = "Adjusted p-value",
    title = "Enriched MF in testes Hypermethylated Genes"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(color = "black"),
    axis.text.y = element_text(size = 10)
  )


####Only with BP
# Filtrar solo términos GO:BP y los más significativos (por ejemplo, top 15)
bp_df <- go_df %>%
  filter(source == "GO:BP") %>%
  arrange(adjusted_p_value) %>%
  slice_head(n = 10) %>%
  mutate(term_name = factor(term_name, levels = rev(term_name)))


# Crear gráfico
ggplot(bp_df, aes(x = intersection_size, y = term_name, fill = adjusted_p_value)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "red3", high = "blue3", trans = "log10", guide = guide_colorbar(reverse = TRUE)) +
  labs(
    x = "Gene Count",
    y = NULL,
    fill = "Adjusted p-value",
    title = "Top enriched Biological Processes\nHypermethylated genes in testes"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(color = "black"),
    axis.text.y = element_text(size = 10)
  )

##FINAL VERSION
ggplot(bp_df, aes(x = intersection_size, y = term_name, fill = adjusted_p_value)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(
    low = "red3",
    high = "blue3",
    trans = "log10",
    labels = scales::label_number(accuracy = 0.0001),
    guide = guide_colorbar(reverse = TRUE)
  ) +
  labs(
    x = "Gene Count",
    y = NULL,
    fill = "Adjusted p-value",
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(color = "black"),
    axis.text.y = element_text(size = 10)
  )


library(scales)

ggplot(bp_df, aes(x = intersection_size, y = term_name, fill = adjusted_p_value)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = intersection_size), 
            hjust = -0.1, size = 4.5, color = "black") +  # Ajusta el tamaño y la posición del texto
  scale_fill_gradient(
    low = "#a06",
    high = "blue3",
    trans = "log10",
    labels = label_number(accuracy = 0.0001),
    guide = guide_colorbar(reverse = TRUE)
  ) +
  labs(
    x = "Gene count",
    y = NULL,
    fill = "adjusted p-value"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(color = "black", size = 13),  # Tamaño eje X
    axis.text.y = element_text(color = "black", size = 13),  # Tamaño eje Y (términos GO)
    axis.title.x = element_text(color = "black", size = 14),
    axis.title.y = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(color = "black")
  ) +
  coord_cartesian(xlim = c(0, max(bp_df$intersection_size) * 1.1))  # Para que haya espacio para los números


