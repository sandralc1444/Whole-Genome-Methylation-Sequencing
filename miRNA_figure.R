# ---- Paquetes ----
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(readxl)
library(forcats)

# ==== 1) LECTURA ====
df_raw <- read_excel("miRNA_figure.xlsx", sheet = 2)

# Normaliza nombres de columnas (por si en el Excel pone Testicles/Testes, etc.)
names(df_raw) <- names(df_raw) %>%
  str_trim() %>%
  str_replace_all("(?i)^mirna$", "miRNA") %>%
  str_replace_all("(?i)^ovar(y|ies)$", "Ovaries") %>%
  str_replace_all("(?i)^test(e|i)s$|(?i)^testicles$", "Testes")

# Verifica columnas mínimas
stopifnot(all(c("miRNA", "Ovaries", "Testes") %in% names(df_raw)))

# ==== 2) FORMATO LARGO ====
df_long <- df_raw %>%
  pivot_longer(
    cols = c(Ovaries, Testes),
    names_to = "Tissue",
    values_to = "State"
  ) %>%
  mutate(
    Tissue = factor(Tissue, levels = c("Ovaries", "Testes")),
    State = str_to_lower(str_trim(as.character(State))),
    State = case_when(
      State %in% c("hypo", "hypomethylated", "hypomethylation") ~ "hypo",
      State %in% c("hyper", "hypermethylated", "hypermethylation") ~ "hyper",
      State %in% c("no", "none", "ns", "no diff", "sin diferencia") ~ "no",
      TRUE ~ "no"  # por seguridad
    )
  )

# Mantener el orden de miRNA tal como aparece en el Excel (primero arriba)
mirna_order <- df_raw$miRNA %>% as.character()
df_long <- df_long %>%
  mutate(miRNA = factor(miRNA, levels = rev(unique(mirna_order))))

# ==== 3) COLORES ====
fill_vals <- c(
  "hyper" = "#548B54",  # verde
  "hypo"  = "#8B4789",  # lila
  "no"    = "#FFFFFF"   # blanco
)

label_vals <- c(
  "hyper" = "Hypermethylation",
  "hypo"  = "Hypomethylation",
  "no"    = "No difference"
)

# ==== 4) PLOT ====
p <- ggplot(df_long, aes(x = Tissue, y = miRNA, fill = State)) +
  geom_tile(color = "grey40", width = 0.92, height = 0.92) +
  scale_fill_manual(
    values = fill_vals,
    breaks = c("hypo", "hyper", "no"),  # orden de leyenda (puedes cambiarlo)
    labels = label_vals,
    name = "Methylation"
  ) +
  labs(
    x = NULL, y = NULL,
    title = "miRNA methylation status by sex"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

print(p)

####Basal conditions####
library(dplyr)
library(stringr)
library(ggplot2)
library(readxl)

# ==== 1) LECTURA ====
df_raw <- read_excel("miRNA_figure.xlsx", sheet = 2)

# Normaliza nombres de columnas
names(df_raw) <- names(df_raw) %>%
  str_trim() %>%
  str_replace_all("(?i)^mirna$", "miRNA") %>%
  str_replace_all("(?i)^em[-_ ]?seq$", "Em-seq")

# Verifica columnas mínimas
stopifnot(all(c("miRNA", "Em-seq") %in% names(df_raw)))

# ==== 2) LIMPIEZA / ESTADOS ====
df_plot <- df_raw %>%
  mutate(
    `Em-seq` = str_to_lower(str_trim(as.character(`Em-seq`))),
    State = case_when(
      `Em-seq` %in% c("hypo", "hypomethylated", "hypomethylation") ~ "hypo",
      `Em-seq` %in% c("hyper", "hypermethylated", "hypermethylation") ~ "hyper",
      `Em-seq` %in% c("hypo and hyper", "hyper and hypo", "both", "mixed") ~ "both",
      `Em-seq` %in% c("no", "none", "ns", "no diff", "no difference", "sin diferencia") ~ "no",
      TRUE ~ "no"
    )
  )

# Mantener el orden del Excel (primer miRNA arriba)
mirna_order <- df_raw$miRNA %>% as.character()
df_plot <- df_plot %>%
  mutate(
    miRNA = factor(miRNA, levels = rev(unique(mirna_order))),
    X = "Em-seq"
  )

# ==== 3) COLORES ====
# Mantengo tus colores y añado uno intermedio/daltónico-friendly para "both"
fill_vals <- c(
  "hyper" = "#548B54",  # verde
  "both"  = "#4C566A",  # gris azulado oscuro (muy visible y apto para daltónicos)
  "hypo"  = "#8B4789",  # lila
  "no"    = "#FFFFFF"   # blanco
)

label_vals <- c(
  "hyper" = "Hypermethylation",
  "both"  = "Hypo and hyper",
  "hypo"  = "Hypomethylation",
  "no"    = "No difference"
)

# ==== 4) PLOT ====
p <- ggplot(df_plot, aes(x = X, y = miRNA, fill = State)) +
  geom_tile(color = "grey40", width = 0.2, height = 0.90) +  # columna fina
  scale_fill_manual(
    values = fill_vals,
    breaks = c("hypo", "both", "hyper", "no"),
    labels = label_vals,
    name = "EM-seq"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(
    x = NULL, y = NULL,
    title = "Differential methylation in miRNA-associated regions (EM-seq)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

print(p)

# ==== 5) GUARDAR (opcional) ====
ggsave("miRNA_EMseq_one_column.png", p, width = 4.5, height = 7, dpi = 300)





