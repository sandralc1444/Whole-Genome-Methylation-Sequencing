# ---- Paquetes ----
library(dplyr)
library(tidyr)
library(stringr)
library(forcats)
library(ggplot2)
library(readxl)

# ==== 1) LECTURA ====
df_raw <- read_excel("summary_figure.xlsx", sheet = 1)

# Normaliza nombres de columnas
names(df_raw) <- names(df_raw) |>
  str_replace_all("(?i)^em[-_ ]?seq$", "Em-seq") |>
  str_replace_all("(?i)^rna[-_ ]?seq?$", "RNA") |>
  str_replace_all("(?i)^qpcr$", "qPCR") |>
  str_replace_all("(?i)^gene$", "Gene") |>
  str_replace_all("(?i)^tissue$", "Tissue")

stopifnot(all(c("Gene","Em-seq","RNA","qPCR","Tissue") %in% names(df_raw)))

# Homogeneiza valores
df <- df_raw %>%
  mutate(
    Tissue = case_when(
      str_detect(Tissue, "(?i)^ovar") ~ "Ovaries",
      str_detect(Tissue, "(?i)^test") ~ "Testes",
      TRUE ~ Tissue
    ),
    `Em-seq` = str_to_lower(`Em-seq`),
    RNA      = str_to_lower(RNA),
    qPCR     = str_to_lower(qPCR)
  )

# ==== 2) LARGO + COMPLETAR SEXO FALTANTE ====
long <- df %>%
  pivot_longer(cols = c(`Em-seq`, RNA, qPCR),
               names_to = "Layer", values_to = "State") %>%
  mutate(Layer = recode(Layer, "Em-seq" = "EM-seq", "RNA" = "RNA-seq"))

# Crea contraparte en el otro tejido como "no" (para tener siempre 6 columnas)
mirror <- long %>%
  mutate(Tissue = if_else(Tissue == "Ovaries", "Testes", "Ovaries"),
         State  = "no")

long2 <- bind_rows(long, mirror) %>%
  distinct(Gene, Tissue, Layer, .keep_all = TRUE)

# ==== 3) ORDEN EJE X Y GENES ====
long2 <- long2 %>%
  mutate(
    Tissue  = factor(Tissue, levels = c("Ovaries","Testes")),
    Layer   = factor(Layer,  levels = c("EM-seq","RNA-seq","qPCR")),
    X       = interaction(Tissue, Layer, sep = " | ", lex.order = TRUE)
  )
x_levels <- c("Ovaries | EM-seq","Ovaries | RNA-seq","Ovaries | qPCR",
              "Testes | EM-seq","Testes | RNA-seq","Testes | qPCR")
long2$X <- factor(long2$X, levels = x_levels)

# Orden de genes: primero los que venían de Ovarios, luego Testes; respeta el orden de tu hoja
genes_order <- df %>%
  group_by(Gene) %>% summarise(source = first(Tissue), .groups = "drop") %>%
  arrange(factor(source, levels = c("Ovaries","Testes"))) %>%
  pull(Gene)
long2 <- long2 %>% mutate(Gene = factor(Gene, levels = rev(unique(genes_order)))) # rev para que el primero quede arriba

# ==== 4) LEYENDA (CATEGORÍAS) Y COLORES ====
# Definimos categorías para poder tener leyenda
long2 <- long2 %>%
  mutate(
    Cat = case_when(
      Layer == "EM-seq" & State == "hyper" ~ "EM: hyper",
      Layer == "EM-seq" & State == "hypo"  ~ "EM: hypo",
      Layer == "EM-seq" & State == "no"    ~ "no",
      Layer != "EM-seq" & State == "up"    ~ "RNA/qPCR: up",
      Layer != "EM-seq" & State == "down"  ~ "RNA/qPCR: down",
      TRUE                                 ~ "no"
    )
  )

fill_vals <- c(
  "EM: hyper"      = "#548B54",  # verde
  "EM: hypo"       = "#8B4789",  # lila
  "RNA/qPCR: up"   = "#E9C46A",  # naranja
  "RNA/qPCR: down" = "#2F8DBA",  # azul
  "no"             = "#FFFFFF"   # blanco
)

label_vals <- c(
  "EM: hyper"      = "EM-seq: hyper",
  "EM: hypo"       = "EM-seq: hypo",
  "RNA/qPCR: up"   = "RNA/qPCR: up",
  "RNA/qPCR: down" = "RNA/qPCR: down",
  "no"             = "Sin diferencia"
)

# ==== 5) CONCORDANCIA RNA vs qPCR (✓/✕ sobre qPCR) ====
# Construimos tabla ancha por Gene+Tissue para comparar
wide_cmp <- long2 %>%
  select(Gene, Tissue, Layer, State) %>%
  pivot_wider(names_from = Layer, values_from = State) %>%
  mutate(
    agree = case_when(
      `RNA-seq` %in% c("up","down") & qPCR %in% c("up","down") &
        (`RNA-seq` == qPCR) ~ "ok",
      `RNA-seq` %in% c("up","down") & qPCR %in% c("up","down") &
        (`RNA-seq` != qPCR) ~ "bad",
      TRUE ~ NA_character_
    )
  )

# Nos quedamos con la marca solo para colocarla sobre las celdas qPCR
marks <- wide_cmp %>%
  mutate(X = factor(paste(Tissue, "qPCR", sep=" | "), levels = x_levels)) %>%
  select(Gene, X, agree) %>%
  mutate(mark = ifelse(agree == "ok", "✓",
                       ifelse(agree == "bad", "✕", NA)))

# ==== 6) PLOT ====
p <- ggplot(long2, aes(x = X, y = Gene)) +
  geom_tile(aes(fill = Cat), color = "grey40", width = 0.92, height = 0.92) +  # columnas algo más estrechas
  scale_fill_manual(values = fill_vals, breaks = names(label_vals), labels = label_vals, name = "Dirección") +
  # marcas ✓/✕ sobre qPCR
  geom_text(
    data = marks,
    aes(x = X, y = Gene, label = mark),
    na.rm = TRUE, size = 3.5, fontface = "bold"
  ) +
  # separador y rótulos de sexo
  geom_vline(xintercept = 3.5, color = "black", linewidth = 0.7) +
  annotate("text", x = 2,   y = Inf, label = "Ovaries", vjust = 1.4, fontface = "bold") +
  annotate("text", x = 5,   y = Inf, label = "Testes",  vjust = 1.4, fontface = "bold") +
  scale_x_discrete(labels = c("EM-seq","RNA-seq","qPCR","EM-seq","RNA-seq","qPCR")) +
  labs(x = NULL, y = NULL,
       title = "EM-seq, RNA-seq y qPCR por sexo",
       subtitle = "Colores por técnica y dirección; ✓ = concordancia RNA vs qPCR; ✕ = discrepancia") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(face = "italic"),  # genes en cursiva
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

print(p)

# ==== 7) EXPORTAR (columnas más estrechas) ====
n_genes <- length(levels(long2$Gene))
# Altura en función del nº de genes; ancho contenido para que las columnas no “engorden”
ggsave("fig_sexo_methyl_rna_qpcr.pdf", p, width = 5.2, height = 0.35*n_genes + 1.6, units = "in")
ggsave("fig_sexo_methyl_rna_qpcr.png", p, width = 5.2, height = 0.35*n_genes + 1.6, units = "in", dpi = 300)


##Sexual dimorphism comparison
# ---- Paquetes ----
library(dplyr)
library(tidyr)
library(stringr)
library(forcats)
library(ggplot2)
library(readxl)

# ==== 1) LECTURA ====
df_raw <- read_excel("summary_figure.xlsx", sheet = 2)

# Normaliza nombres de columnas
names(df_raw) <- names(df_raw) |>
  str_replace_all("(?i)^em[-_ ]?seq$", "Em-seq") |>
  str_replace_all("(?i)^rna[-_ ]?seq?$", "RNA") |>
  str_replace_all("(?i)^qpcr$", "qPCR") |>
  str_replace_all("(?i)^gene$", "Gene")

stopifnot(all(c("Gene","Em-seq","RNA","qPCR") %in% names(df_raw)))

# Homogeneiza valores
df <- df_raw %>%
  transmute(
    Gene,
    `Em-seq` = str_to_lower(`Em-seq`),
    RNA      = str_to_lower(RNA),
    qPCR     = str_to_lower(qPCR)
  )

# ==== 2) ORDEN MANUAL DE GENES ====
orden_manual <- c(
  "akna",
  "aqp1a",
  "impdh",
  "map4k",
  "nod1",
  "ptprz",
  "pxdn",
  "sting1",
  "vash2",
  "zfp64"
)

# ==== 3) FORMATO LARGO (3 columnas) + aplicar orden ====
long2 <- df %>%
  pivot_longer(cols = c(`Em-seq`, RNA, qPCR),
               names_to = "Layer", values_to = "State") %>%
  mutate(
    Layer = recode(Layer, "Em-seq" = "EM-seq", "RNA" = "RNA-seq"),
    Layer = factor(Layer, levels = c("EM-seq","RNA-seq","qPCR")),
    Gene  = factor(Gene, levels = rev(orden_manual))  # primero de tu lista arriba
  )

# ==== 3) COLORES Y LEYENDA ====
long2 <- long2 %>%
  mutate(
    Cat = case_when(
      Layer == "EM-seq" & State == "hyper" ~ "EM: hyper",
      Layer == "EM-seq" & State == "hypo"  ~ "EM: hypo",
      Layer == "EM-seq" & State == "no"    ~ "no",
      Layer != "EM-seq" & State == "up"    ~ "RNA/qPCR: up",
      Layer != "EM-seq" & State == "down"  ~ "RNA/qPCR: down",
      TRUE                                 ~ "no"
    )
  )

fill_vals <- c(
  "EM: hyper"      = "#548B54",  # verde
  "EM: hypo"       = "#8B4789",  # lila
  "RNA/qPCR: up"   = "#E9C46A",  # amarillo mostaza
  "RNA/qPCR: down" = "#2F8DBA",  # azul
  "no"             = "#FFFFFF"   # blanco
)

label_vals <- c(
  "EM: hyper"      = "EM-seq: hyper",
  "EM: hypo"       = "EM-seq: hypo",
  "RNA/qPCR: up"   = "RNA/qPCR: up",
  "RNA/qPCR: down" = "RNA/qPCR: down",
  "no"             = "Sin diferencia"
)

# ==== 4) CONCORDANCIA RNA vs qPCR (✓/✕ SOLO en qPCR) ====
wide_cmp <- df %>%
  mutate(
    RNA = recode(str_to_lower(RNA), "rna-seq"="RNA-seq"),
    qPCR = str_to_lower(qPCR)
  ) %>%
  mutate(
    agree = case_when(
      RNA %in% c("up","down") & qPCR %in% c("up","down") & RNA == qPCR ~ "ok",
      RNA %in% c("up","down") & qPCR %in% c("up","down") & RNA != qPCR ~ "bad",
      TRUE ~ NA_character_
    )
  )

marks <- wide_cmp %>%
  transmute(
    Gene = factor(Gene, levels = levels(long2$Gene)),
    Layer = factor("qPCR", levels = c("EM-seq","RNA-seq","qPCR")),
    mark = ifelse(agree == "ok", "✓",
                  ifelse(agree == "bad", "✕", NA))
  )

# ==== 5) PLOT (3 columnas, mismo estilo) ====
p <- ggplot(long2, aes(x = Layer, y = Gene)) +
  geom_tile(aes(fill = Cat), color = "grey40", linewidth = 0.25,
            width = 0.92, height = 0.92) +
  scale_fill_manual(values = fill_vals,
                    breaks = names(label_vals),
                    labels = label_vals,
                    name = "Dirección") +
  geom_text(
    data = marks,
    aes(x = Layer, y = Gene, label = mark),
    na.rm = TRUE, size = 3.5, fontface = "bold"
  ) +
  scale_x_discrete(labels = c("EM-seq","RNA-seq","qPCR")) +
  labs(x = NULL, y = NULL,
       title = "EM-seq, RNA-seq y qPCR",
       subtitle = "Colores por técnica y dirección; ✓ = concordancia RNA vs qPCR; ✕ = discrepancia") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(face = "italic"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

print(p)

# ==== 6) EXPORTAR ====
n_genes <- length(levels(long2$Gene))
ggsave("fig_emseq_rnaseq_qpcr_singlecomp.pdf", p,
       width = 4.8, height = 0.35*n_genes + 1.6, units = "in")
ggsave("fig_emseq_rnaseq_qpcr_singlecomp.png", p,
       width = 4.8, height = 0.35*n_genes + 1.6, units = "in", dpi = 300)



##Without RNA-seq columns##
# ---- Paquetes ----
library(dplyr)
library(tidyr)
library(stringr)
library(forcats)
library(ggplot2)
library(readxl)

# ==== 1) LECTURA ====
df_raw <- read_excel("summary_figure.xlsx", sheet = 1)

# Normaliza nombres de columnas
names(df_raw) <- names(df_raw) |>
  str_replace_all("(?i)^em[-_ ]?seq$", "Em-seq") |>
  str_replace_all("(?i)^rna[-_ ]?seq?$", "RNA") |>
  str_replace_all("(?i)^qpcr$", "qPCR") |>
  str_replace_all("(?i)^gene$", "Gene") |>
  str_replace_all("(?i)^tissue$", "Tissue")

stopifnot(all(c("Gene","Em-seq","RNA","qPCR","Tissue") %in% names(df_raw)))

# Homogeneiza valores
df <- df_raw %>%
  mutate(
    Tissue = case_when(
      str_detect(Tissue, "(?i)^ovar") ~ "Ovaries",
      str_detect(Tissue, "(?i)^test") ~ "Testes",
      TRUE ~ Tissue
    ),
    `Em-seq` = str_to_lower(`Em-seq`),
    qPCR     = str_to_lower(qPCR)
  )

# ==== 2) LARGO (solo EM-seq y qPCR) + COMPLETAR SEXO FALTANTE ====
long <- df %>%
  pivot_longer(cols = c(`Em-seq`, qPCR),
               names_to = "Layer", values_to = "State") %>%
  mutate(Layer = recode(Layer, "Em-seq" = "EM-seq"))

# Crea contraparte en el otro tejido como "no" (para tener siempre 4 columnas)
mirror <- long %>%
  mutate(Tissue = if_else(Tissue == "Ovaries", "Testes", "Ovaries"),
         State  = "no")

long2 <- bind_rows(long, mirror) %>%
  distinct(Gene, Tissue, Layer, .keep_all = TRUE)

# ==== 3) ORDEN EJE X Y GENES ====
long2 <- long2 %>%
  mutate(
    Tissue  = factor(Tissue, levels = c("Ovaries","Testes")),
    Layer   = factor(Layer,  levels = c("EM-seq","qPCR")),
    X       = interaction(Tissue, Layer, sep = " | ", lex.order = TRUE)
  )

x_levels <- c("Ovaries | EM-seq","Ovaries | qPCR",
              "Testes | EM-seq","Testes | qPCR")
long2$X <- factor(long2$X, levels = x_levels)

# Orden de genes (respeta orden de tu hoja; ovarios primero)
genes_order <- df %>%
  group_by(Gene) %>% summarise(source = first(Tissue), .groups = "drop") %>%
  arrange(factor(source, levels = c("Ovaries","Testes"))) %>%
  pull(Gene)
long2 <- long2 %>% mutate(Gene = factor(Gene, levels = rev(unique(genes_order))))

# ==== 4) LEYENDA Y COLORES ====
long2 <- long2 %>%
  mutate(
    Cat = case_when(
      Layer == "EM-seq" & State == "hyper" ~ "EM: hyper",
      Layer == "EM-seq" & State == "hypo"  ~ "EM: hypo",
      Layer == "EM-seq" & State == "no"    ~ "no",
      Layer == "qPCR"   & State == "up"    ~ "qPCR: up",
      Layer == "qPCR"   & State == "down"  ~ "qPCR: down",
      TRUE                                 ~ "no"
    )
  )

fill_vals <- c(
  "EM: hyper"  = "#548B54",  # verde
  "EM: hypo"   = "#8B4789",  # lila
  "qPCR: up"   = "#E9C46A",  # naranja
  "qPCR: down" = "#2F8DBA",  # azul
  "no"         = "#FFFFFF"   # blanco
)

label_vals <- c(
  "EM: hyper"  = "EM-seq: hyper",
  "EM: hypo"   = "EM-seq: hypo",
  "qPCR: up"   = "qPCR: up",
  "qPCR: down" = "qPCR: down",
  "no"         = "Sin diferencia"
)

# ==== 5) PLOT (sin rótulos de sexo ni marcas internas) ====
p <- ggplot(long2, aes(x = X, y = Gene)) +
  geom_tile(aes(fill = Cat), color = "grey90", width = 0.85, height = 0.92) +
  scale_fill_manual(values = fill_vals, breaks = names(label_vals),
                    labels = label_vals, name = "Dirección") +
  geom_vline(xintercept = 2.5, color = "grey55", linewidth = 0.7) +  # separador entre sexos
  scale_x_discrete(labels = c("EM-seq","qPCR","EM-seq","qPCR")) +
  labs(x = NULL, y = NULL,
       title = "EM-seq y qPCR por sexo",
       subtitle = "Colores por técnica y dirección") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(face = "italic"),  # genes en cursiva
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

print(p)

# ==== 6) EXPORTAR ====
n_genes <- length(levels(long2$Gene))
ggsave("fig_emseq_qpcr_sexo.pdf", p, width = 4.6, height = 0.35*n_genes + 1.4, units = "in")
ggsave("fig_emseq_qpcr_sexo.png", p, width = 4.6, height = 0.35*n_genes + 1.4, units = "in", dpi = 300)


##Ordenando los genes
# ---- Paquetes ----
library(dplyr)
library(tidyr)
library(stringr)
library(forcats)
library(ggplot2)
library(readxl)

# ==== 1) LECTURA ====
df_raw <- read_excel("summary_figure.xlsx", sheet = 1)

# Normaliza nombres de columnas
names(df_raw) <- names(df_raw) |>
  str_replace_all("(?i)^em[-_ ]?seq$", "Em-seq") |>
  str_replace_all("(?i)^rna[-_ ]?seq?$", "RNA") |>
  str_replace_all("(?i)^qpcr$", "qPCR") |>
  str_replace_all("(?i)^gene$", "Gene") |>
  str_replace_all("(?i)^tissue$", "Tissue")

stopifnot(all(c("Gene","Em-seq","RNA","qPCR","Tissue") %in% names(df_raw)))

# Homogeneiza valores
df <- df_raw %>%
  mutate(
    Tissue = case_when(
      str_detect(Tissue, "(?i)^ovar") ~ "Ovaries",
      str_detect(Tissue, "(?i)^test") ~ "Testes",
      TRUE ~ Tissue
    ),
    `Em-seq` = str_to_lower(`Em-seq`),
    qPCR     = str_to_lower(qPCR)
  )

# ==== 2) ORDEN MANUAL DE GENES (mismo que antes) ====
orden_manual <- c(
  "akna",
  "aqp1a",
  "impdh",
  "map4k",
  "nod1",
  "ptprz",
  "pxdn",
  "sting1",
  "vash2",
  "zfp64"
)

# ==== 3) LARGO (solo EM-seq y qPCR) + COMPLETAR SEXO FALTANTE ====
long <- df %>%
  pivot_longer(cols = c(`Em-seq`, qPCR),
               names_to = "Layer", values_to = "State") %>%
  mutate(Layer = recode(Layer, "Em-seq" = "EM-seq"))

# Crea contraparte en el otro tejido como "no" (para tener siempre 4 columnas)
mirror <- long %>%
  mutate(Tissue = if_else(Tissue == "Ovaries", "Testes", "Ovaries"),
         State  = "no")

long2 <- bind_rows(long, mirror) %>%
  distinct(Gene, Tissue, Layer, .keep_all = TRUE)

# ==== 4) ORDEN EJE X + aplicar orden manual de genes ====
long2 <- long2 %>%
  mutate(
    Tissue  = factor(Tissue, levels = c("Ovaries","Testes")),
    Layer   = factor(Layer,  levels = c("EM-seq","qPCR")),
    X       = interaction(Tissue, Layer, sep = " | ", lex.order = TRUE),
    Gene    = factor(Gene, levels = rev(orden_manual))  # primero de tu lista arriba
  )

x_levels <- c("Ovaries | EM-seq","Ovaries | qPCR",
              "Testes | EM-seq","Testes | qPCR")
long2$X <- factor(long2$X, levels = x_levels)

# ==== 5) LEYENDA Y COLORES ====
long2 <- long2 %>%
  mutate(
    Cat = case_when(
      Layer == "EM-seq" & State == "hyper" ~ "EM: hyper",
      Layer == "EM-seq" & State == "hypo"  ~ "EM: hypo",
      Layer == "EM-seq" & State == "no"    ~ "no",
      Layer == "qPCR"   & State == "up"    ~ "qPCR: up",
      Layer == "qPCR"   & State == "down"  ~ "qPCR: down",
      TRUE                                 ~ "no"
    )
  )

fill_vals <- c(
  "EM: hyper"  = "#548B54",  # verde
  "EM: hypo"   = "#8B4789",  # lila
  "qPCR: up"   = "#E9C46A",  # amarillo mostaza
  "qPCR: down" = "#2F8DBA",  # azul
  "no"         = "#FFFFFF"   # blanco
)

label_vals <- c(
  "EM: hyper"  = "EM-seq: hyper",
  "EM: hypo"   = "EM-seq: hypo",
  "qPCR: up"   = "qPCR: up",
  "qPCR: down" = "qPCR: down",
  "no"         = "Sin diferencia"
)

# ==== 6) PLOT ====
p <- ggplot(long2, aes(x = X, y = Gene)) +
  geom_tile(aes(fill = Cat), color = "black", width = 0.85, height = 0.92) +
  scale_fill_manual(values = fill_vals, breaks = names(label_vals),
                    labels = label_vals, name = "Dirección") +
  geom_vline(xintercept = 2.5, color = "grey55", linewidth = 0.7) +
  scale_x_discrete(labels = c("EM-seq","qPCR","EM-seq","qPCR")) +
  labs(x = NULL, y = NULL,
       title = "EM-seq y qPCR por sexo",
       subtitle = "Colores por técnica y dirección") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(face = "italic"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

print(p)

# ==== 7) EXPORTAR ====
n_genes <- length(levels(long2$Gene))
ggsave("fig_emseq_qpcr_sexo.pdf", p, width = 4.6, height = 0.35*n_genes + 1.4, units = "in")
ggsave("fig_emseq_qpcr_sexo.png", p, width = 4.6, height = 0.35*n_genes + 1.4, units = "in", dpi = 300)




##Sexual dimorphism --> only Em-seq y qPCR
# ---- Paquetes ----
library(dplyr)
library(tidyr)
library(stringr)
library(forcats)
library(ggplot2)
library(readxl)

# ==== 1) LECTURA (sin Tissue) ====
df_raw <- read_excel("summary_figure.xlsx", sheet = 2)

# Normaliza nombres de columnas (por si vienen con variantes)
names(df_raw) <- names(df_raw) |>
  str_replace_all("(?i)^em[-_ ]?seq$", "Em-seq") |>
  str_replace_all("(?i)^rna[-_ ]?seq?$", "RNA") |>
  str_replace_all("(?i)^qpcr$", "qPCR") |>
  str_replace_all("(?i)^gene$", "Gene")

# Comprobación mínima
stopifnot(all(c("Gene","Em-seq","qPCR") %in% names(df_raw)))

# Homogeneiza valores
df <- df_raw %>%
  transmute(
    Gene,
    `Em-seq` = str_to_lower(`Em-seq`),
    qPCR     = str_to_lower(qPCR)
  )

# ==== 2) A FORMATO LARGO (solo EM-seq y qPCR) ====
long2 <- df %>%
  pivot_longer(cols = c(`Em-seq`, qPCR),
               names_to = "Layer", values_to = "State") %>%
  mutate(
    Layer = recode(Layer, "Em-seq" = "EM-seq"),
    # orden de los ejes
    Layer = factor(Layer, levels = c("EM-seq","qPCR")),
    Gene  = factor(Gene, levels = rev(unique(df$Gene))) # primero arriba
  )

# ==== 3) CATEGORÍAS Y COLORES ====
long2 <- long2 %>%
  mutate(
    Cat = case_when(
      Layer == "EM-seq" & State == "hyper" ~ "EM: hyper",
      Layer == "EM-seq" & State == "hypo"  ~ "EM: hypo",
      Layer == "EM-seq" & State == "no"    ~ "no",
      Layer == "qPCR"   & State == "up"    ~ "qPCR: up",
      Layer == "qPCR"   & State == "down"  ~ "qPCR: down",
      TRUE                                 ~ "no"
    )
  )

fill_vals <- c(
  "EM: hyper"  = "#548B54",  # verde
  "EM: hypo"   = "#8B4789",  # lila
  "qPCR: up"   = "#E9C46A",  # amarillo mostaza
  "qPCR: down" = "#2F8DBA",  # azul
  "no"         = "#FFFFFF"   # blanco
)



label_vals <- c(
  "EM: hyper"  = "EM-seq: hyper",
  "EM: hypo"   = "EM-seq: hypo",
  "qPCR: up"   = "qPCR: up",
  "qPCR: down" = "qPCR: down",
  "no"         = "Sin diferencia"
)

# ==== 4) PLOT (dos columnas) ====
p <- ggplot(long2, aes(x = Layer, y = Gene)) +
  geom_tile(aes(fill = Cat), color = "grey90", width = 0.8, height = 0.9) +
  scale_fill_manual(values = fill_vals, breaks = names(label_vals),
                    labels = label_vals, name = "Dirección") +
  labs(x = NULL, y = NULL,
       title = "EM-seq y qPCR (♀ control vs ♂ control)",
       subtitle = "Colores por técnica y dirección") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(face = "italic"),  # genes en cursiva
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

print(p)

# ==== 5) EXPORTAR ====
n_genes <- length(levels(long2$Gene))
ggsave("fig_emseq_qpcr_sex_ctrls.pdf", p, width = 3.8, height = 0.35*n_genes + 1.2, units = "in")
ggsave("fig_emseq_qpcr_sex_ctrls.png", p, width = 3.8, height = 0.35*n_genes + 1.2, units = "in", dpi = 300)

#Ordenando los genes
##Sexual dimorphism --> only EM-seq y qPCR (con orden manual y misma paleta)

# ---- Paquetes ----
library(dplyr)
library(tidyr)
library(stringr)
library(forcats)
library(ggplot2)
library(readxl)

# ==== 1) LECTURA (sin Tissue) ====
df_raw <- read_excel("summary_figure.xlsx", sheet = 2)

# Normaliza nombres de columnas
names(df_raw) <- names(df_raw) |>
  str_replace_all("(?i)^em[-_ ]?seq$", "Em-seq") |>
  str_replace_all("(?i)^rna[-_ ]?seq?$", "RNA") |>
  str_replace_all("(?i)^qpcr$", "qPCR") |>
  str_replace_all("(?i)^gene$", "Gene")

stopifnot(all(c("Gene","Em-seq","qPCR") %in% names(df_raw)))

# Homogeneiza valores
df <- df_raw %>%
  transmute(
    Gene,
    `Em-seq` = str_to_lower(`Em-seq`),
    qPCR     = str_to_lower(qPCR)
  )

# ==== 2) ORDEN MANUAL DE GENES (mismo que antes) ====
orden_manual <- c(
  "akna",
  "aqp1a",
  "impdh",
  "map4k",
  "nod1",
  "ptprz",
  "pxdn",
  "sting1",
  "vash2",
  "zfp64"
)

# ==== 3) A FORMATO LARGO (solo EM-seq y qPCR) + aplicar orden ====
long2 <- df %>%
  pivot_longer(cols = c(`Em-seq`, qPCR),
               names_to = "Layer", values_to = "State") %>%
  mutate(
    Layer = recode(Layer, "Em-seq" = "EM-seq"),
    Layer = factor(Layer, levels = c("EM-seq","qPCR")),
    Gene  = factor(Gene, levels = rev(orden_manual))  # primero arriba
  )

# ==== 4) CATEGORÍAS Y COLORES ====
long2 <- long2 %>%
  mutate(
    Cat = case_when(
      Layer == "EM-seq" & State == "hyper" ~ "EM: hyper",
      Layer == "EM-seq" & State == "hypo"  ~ "EM: hypo",
      Layer == "EM-seq" & State == "no"    ~ "no",
      Layer == "qPCR"   & State == "up"    ~ "qPCR: up",
      Layer == "qPCR"   & State == "down"  ~ "qPCR: down",
      TRUE                                 ~ "no"
    )
  )

fill_vals <- c(
  "EM: hyper"  = "#548B54",  # verde
  "EM: hypo"   = "#8B4789",  # lila
  "qPCR: up"   = "#E9C46A",  # amarillo mostaza
  "qPCR: down" = "#2F8DBA",  # azul
  "no"         = "#FFFFFF"   # blanco
)

label_vals <- c(
  "EM: hyper"  = "EM-seq: hyper",
  "EM: hypo"   = "EM-seq: hypo",
  "qPCR: up"   = "qPCR: up",
  "qPCR: down" = "qPCR: down",
  "no"         = "Sin diferencia"
)

# ==== 5) PLOT (dos columnas) ====
p <- ggplot(long2, aes(x = Layer, y = Gene)) +
  geom_tile(aes(fill = Cat), color = "grey40", linewidth = 0.25,
            width = 0.8, height = 0.9) +
  scale_fill_manual(values = fill_vals, breaks = names(label_vals),
                    labels = label_vals, name = "Dirección") +
  labs(x = NULL, y = NULL,
       title = "EM-seq y qPCR (♀ control vs ♂ control)",
       subtitle = "Colores por técnica y dirección") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(face = "italic"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

print(p)

# ==== 6) EXPORTAR ====
n_genes <- length(levels(long2$Gene))
ggsave("fig_emseq_qpcr_sex_ctrls.pdf", p,
       width = 3.8, height = 0.35*n_genes + 1.2, units = "in")
ggsave("fig_emseq_qpcr_sex_ctrls.png", p,
       width = 3.8, height = 0.35*n_genes + 1.2, units = "in", dpi = 300)


