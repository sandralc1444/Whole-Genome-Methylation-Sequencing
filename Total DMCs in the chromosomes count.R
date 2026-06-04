######Total DMCs in the chromosomes count#####
####All DMCs

library(dplyr)
library(stringr)

# Leer archivo
df <- read.table(
  "DSS_diff_CpG_INF_M_vs_CTRL_M_with_overlap_genes.txt",
  header = TRUE,
  fill = TRUE,
  stringsAsFactors = FALSE
)

# Extraer número de cromosoma
df$chromosome <- as.numeric(
  str_extract(df$chr, "(?<=CAJNNU010000)\\d+(?=\\.1)")
)

# Contar DMCs por cromosoma
chrom_counts <- data.frame(chromosome = 1:24) %>%
  left_join(
    df %>%
      count(chromosome, name = "n_DMCs"),
    by = "chromosome"
  ) %>%
  mutate(n_DMCs = ifelse(is.na(n_DMCs), 0, n_DMCs))

print(chrom_counts)



##Calcular %
chrom_counts <- chrom_counts %>%
  mutate(
    percentage = round(100 * n_DMCs / sum(n_DMCs), 2)
  )

print(chrom_counts)

# Guardar tabla
write.csv(
  chrom_counts,
  "Number of DMCs_testes_all.csv",
  row.names = FALSE
)


##Representation
library(ggplot2)

ggplot(chrom_counts,
       aes(x = factor(chromosome), y = n_DMCs)) +
  geom_bar(stat = "identity") +
  xlab("Chromosome") +
  ylab("Number of DMCs_testes_all") +
  theme_bw()



##Significant DMCs (pval ≤ 0.001 & fdr ≤ 0.001)
###### Total DMCs in the chromosomes count #####

library(dplyr)
library(stringr)
library(ggplot2)

# Leer archivo
df <- read.table(
  "DSS_diff_CpG_INF_M_vs_CTRL_M_with_overlap_genes.txt",
  header = TRUE,
  fill = TRUE,
  stringsAsFactors = FALSE
)

# Filtrar DMCs significativos
df_filtered <- df %>%
  filter(pval <= 0.001 & fdr <= 0.001)

# Extraer número de cromosoma
df_filtered$chromosome <- as.numeric(
  str_extract(df_filtered$chr, "(?<=CAJNNU010000)\\d+(?=\\.1)")
)

# Contar DMCs por cromosoma
chrom_counts <- data.frame(chromosome = 1:24) %>%
  left_join(
    df_filtered %>%
      count(chromosome, name = "n_DMCs"),
    by = "chromosome"
  ) %>%
  mutate(n_DMCs = ifelse(is.na(n_DMCs), 0, n_DMCs))

print(chrom_counts)

# Calcular porcentajes
chrom_counts <- chrom_counts %>%
  mutate(
    percentage = round(100 * n_DMCs / sum(n_DMCs), 2)
  )

print(chrom_counts)

# Guardar tabla
write.csv(
  chrom_counts,
  "Number_of_DMCs_testes_p001_fdr001.csv",
  row.names = FALSE
)

# Representación
ggplot(chrom_counts,
       aes(x = factor(chromosome), y = n_DMCs)) +
  geom_bar(stat = "identity") +
  xlab("Chromosome") +
  ylab("Number of DMCs (p ≤ 0.001, FDR ≤ 0.001)") +
  theme_bw()








