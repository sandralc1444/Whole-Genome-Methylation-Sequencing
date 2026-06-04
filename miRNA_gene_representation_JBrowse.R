library(readr)
library(dplyr)
library(stringr)

make_gene_beds <- function(infile,
                           gene_id,
                           out_prefix = gene_id,
                           fdr_thr = 0.001,
                           pval_thr = 0.001,
                           name_tag = gene_id,
                           strand_out = "+") {
  
  # 1) Leer archivo (separa bien con espacios múltiples o tab)
  df <- read_table2(infile, col_types = cols(.default = col_guess()))
  
  # 2) Comprobaciones mínimas
  needed <- c("chr", "pos", "diff", "pval", "fdr", "gene_id")
  miss <- setdiff(needed, names(df))
  if (length(miss) > 0) stop("Faltan columnas en el DSS: ", paste(miss, collapse = ", "))
  
  # 3) Tipos y filtrado global
  df <- df %>%
    mutate(
      pos  = as.integer(pos),
      diff = as.numeric(diff)
    ) %>%
    filter(!is.na(pos), !is.na(diff)) %>%
    filter(fdr <= fdr_thr, pval <= pval_thr) %>%
    filter(gene_id == !!gene_id)
  
  if (nrow(df) == 0) {
    message("No hay CpGs para gene_id=", gene_id,
            " con fdr<=", fdr_thr, " y pval<=", pval_thr)
    return(invisible(NULL))
  }
  
  # 4) Separación hipo/hiper por diff
  hipo  <- df %>% filter(diff < 0)
  hiper <- df %>% filter(diff > 0)
  
  # 5) Construcción BED (MISMO FORMATO que tu ejemplo):
  #    start = pos ; end = pos + 1 ; score = 0 ; strand fijo "+"
  to_bed <- function(d, label) {
    d %>%
      transmute(
        chr    = as.character(chr),
        start  = as.integer(pos),
        end    = as.integer(pos) + 1L,
        name   = paste0(name_tag, "_", label),
        score  = 0,
        strand = strand_out
      )
  }
  
  bed_hipo  <- to_bed(hipo,  "hipo")
  bed_hiper <- to_bed(hiper, "hiper")
  
  # 6) Guardar
  out_hipo  <- paste0(out_prefix, "_hipo.bed")
  out_hiper <- paste0(out_prefix, "_hiper.bed")
  
  write_tsv(bed_hipo,  out_hipo,  col_names = FALSE)
  write_tsv(bed_hiper, out_hiper, col_names = FALSE)
  
  message("OK: ", out_hipo,  " (n=", nrow(bed_hipo),  ")")
  message("OK: ", out_hiper, " (n=", nrow(bed_hiper), ")")
  
  invisible(list(hipo = bed_hipo, hiper = bed_hiper))
}

# ======= EJECUCIÓN PARA TU GEN =======
make_gene_beds(
  infile     = "DSS_diff_CpG_INF_F_vs_CTRL_F_with_overlap_genes copia.txt",
  gene_id    = "ENSDLAG00005014627",
  out_prefix = "zfp64_FxF",
  fdr_thr    = 0.001,
  pval_thr   = 0.001,
  name_tag   = "zfp64_FxF",
  strand_out = "-"
)


#####Archivo .bed enhancer miR1301#####
library(readr)
enh <- data.frame(
  chr = "CAJNNU010000015.1",
  start = 5081586L,
  end = 5081733L,
  name = "zfp64_openchromatin2",
  score = 0,
  strand = "-"
)

write_tsv(enh, "zfp64_openchromatin2.bed", col_names = FALSE)


#####Archivo .bed EMAR miR1301#####
library(readr)

emar <- data.frame(
  chr = "CAJNNU010000001.1",
  start = 12265446L,
  end = 12265738L,
  name = "Open chromatin_region",
  score = 0,
  strand = "+"
)

write_tsv(emar, "miR-181b_openchromatin_region.bed", col_names = FALSE)

#####Regulatory regions miR130a#####
library(readr)

reg <- data.frame(
  chr   = "CAJNNU010000002.1",
  start = c(18028532L, 18028633L, 18040586L, 18047410L, 18057414L),
  end   = c(18028632L, 18030444L, 18040959L, 18047675L, 18057936L),
  name  = c("promoter","enhancerandEMAR","enhancerandEMAR2","openchromatin","enhancerandEMAR3"),
  score = 0,
  strand= "+"
)

write_tsv(reg, "regulatory_region_mAP4K3.bed", col_names = FALSE)









