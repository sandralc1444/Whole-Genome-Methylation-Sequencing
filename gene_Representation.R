####Gene representation####

library(biomaRt)
packageVersion("biomaRt")
ensembl <- useEnsembl(biomart = "genes")

# Ver datasets disponibles
listDatasets(ensembl)



library(biomaRt)
library(Gviz)
library(GenomicRanges)

options(ucscChromosomeNames = FALSE)   # para que acepte CAJNNU010000024.1

ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "dlabrax_gene_ensembl"
)

##impdh
atts <- listAttributes(ensembl)
atts[grep("transcript", atts$name), ][1:20, ]  # solo para curiosear


tx_id <- "ENSDLAT00005065337"   # sin .2

bm <- getBM(
  attributes = c(
    "chromosome_name",
    "exon_chrom_start", "exon_chrom_end",
    "strand",
    "ensembl_gene_id",
    "external_gene_name",
    "ensembl_transcript_id",
    "rank"
  ),
  filters  = "ensembl_transcript_id",
  values   = tx_id,
  mart     = ensembl
)

head(bm)
unique(bm$chromosome_name)




# Ajusta la ruta al bed
hyper <- read.table("FxF_impdh_hiper.bed", sep = "\t", header = FALSE,
                    col.names = c("chr","start0","end","name","score","strand"))

# BED → 1-based
hyper$pos <- hyper$start0 + 1

# Asegurar que el nombre de cromosoma coincide con Ensembl
hyper$chr <- as.character(hyper$chr)
bm$chromosome_name <- as.character(bm$chromosome_name)

# Si por ejemplo bm tiene "CAJNN010000024.1" y hyper "CAJNN010000024.1", perfecto.
# Si fueran distintos, forzamos:
# hyper$chr <- unique(bm$chromosome_name)

chr_region <- unique(bm$chromosome_name)

## CpGs como DOTS (puntos redondos)
cpg_gr <- GRanges(
  seqnames = hyper$chr,
  ranges   = IRanges(start = hyper$pos, end = hyper$pos)
)

hyper_track <- DataTrack(
  range = cpg_gr,
  data  = rep(1, length(cpg_gr)),   # todos al mismo nivel
  genome = "dlabrax2021",
  chromosome = chr_region,
  name  = "Hyper CpGs",
  type  = "p",                      # <<< PUNTOS
  col   = "#548B54",                # verde
  cex   = 1.2,                      # tamaño del punto (ajusta si quieres)
  ylim  = c(0.5, 1.5),
  yAxis = FALSE,
  background.title = "white"
)




gene_gr <- GRanges(
  seqnames = bm$chromosome_name,
  ranges   = IRanges(start = bm$exon_chrom_start,
                     end   = bm$exon_chrom_end),
  strand   = ifelse(bm$strand[1] == 1, "+", "-"),
  gene     = bm$external_gene_name,
  exon     = bm$rank,
  transcript = bm$ensembl_transcript_id
)

grtrack <- GeneRegionTrack(
  gene_gr,
  genome = "dlabrax2021",
  chromosome = chr_region,
  name = "IMPDH-209",
  showId = TRUE,
  transcriptAnnotation = "transcript"
)

grtrack <- GeneRegionTrack(
  gene_gr,
  genome = "dlabrax2021",
  chromosome = chr_region,
  name = "IMPDH-209",
  showId = FALSE,
  transcriptAnnotation = "transcript",
  stacking = "squish",      # <<< compacta los exones
  exonHeight = 1,        # <<< hace los exones FINOS
  arrowHeadWidth = 8,       # <<< flechas más pequeñas
  geneSymbol = TRUE,
  fill = "#EDD5B3",
  col  = "#EDD5B3",
  background.title = "white",
  cex.title = 0.9
)



axisTrack <- GenomeAxisTrack()

from_region <- min(bm$exon_chrom_start) - 1000
to_region   <- max(bm$exon_chrom_end)   + 1000

plotTracks(
  list(axisTrack, grtrack, hyper_track),
  from = from_region,
  to   = to_region,
  chromosome = chr_region,
  main = "IMPDH-209 (canonical transcript)"
)


#####Añadiendo promotores y enhancers del gen
# Asegúrate de tener esto antes
# chr_region <- unique(bm$chromosome_name)

## PROMOTORES (2 regiones)
promoter_gr <- GRanges(
  seqnames = chr_region,
  ranges = IRanges(
    start = c(15200164, 15195311),   # 15,200,164  y 15,195,311
    end   = c(15200892, 15195442)    # 15,200,892  y 15,195,442
  )
)

## ENHANCERS (3 regiones)
enhancer_gr <- GRanges(
  seqnames = chr_region,
  ranges = IRanges(
    start = c(15195592, 15191731, 15189309),  # 15,195,592; 15,191,731; 15,189,309
    end   = c(15197046, 15192020, 15190254)   # 15,197,046; 15,192,020; 15,190,254
  )
)

promoterTrack <- AnnotationTrack(
  promoter_gr,
  genome = "dlabrax2021",
  chromosome = chr_region,
  name = "Promoters",
  fill = "darkred",
  col  = NA,
  alpha = 0.4      # transparente
)

enhancerTrack <- AnnotationTrack(
  enhancer_gr,
  genome = "dlabrax2021",
  chromosome = chr_region,
  name = "Enhancers",
  fill = "darkorange",
  col  = NA,
  alpha = 0.4
)

from_region <- min(
  c(bm$exon_chrom_start,
    start(promoter_gr),
    start(enhancer_gr))
) - 2000

to_region <- max(
  c(bm$exon_chrom_end,
    end(promoter_gr),
    end(enhancer_gr))
) + 2000



##axisTrack <- GenomeAxisTrack()

plotTracks(
  list(axisTrack, grtrack, promoterTrack, enhancerTrack, hyper_track),
  from = from_region,
  to   = to_region,
  chromosome = chr_region,
  main = "IMPDH-209 (canonical transcript)",
  sizes = c(1, 3, 0.8, 0.8, 1)   # ajusta proporción de cada pista si quieres
)



##Añadiendo la EMAR que se solapa 
chr_region <- unique(bm$chromosome_name)
chr_region
# debería devolver "CAJNNU010000024.1"

library(Gviz)
library(GenomicRanges)

## EMARs (todas las de Ensembl en la región del gen)
emar_starts <- c(15186537L, 15189309L, 15191731L, 15194983L, 15195592L, 15200039L)
emar_ends   <- c(15187265L, 15190254L, 15192020L, 15195419L, 15197046L, 15200163L)

emar_gr <- GRanges(
  seqnames = chr_region,
  ranges   = IRanges(start = emar_starts, end = emar_ends)
)

emarTrack <- AnnotationTrack(
  range = emar_gr,
  genome = "dlabrax2021",
  chromosome = chr_region,
  name = "EMAR",
  fill = "darkgreen",
  col  = NA,
  alpha = 0.4
)

from_region <- min(c(bm$exon_chrom_start,
                     start(promoter_gr),
                     start(enhancer_gr),
                     emar_starts)) - 2000

to_region   <- max(c(bm$exon_chrom_end,
                     end(promoter_gr),
                     end(enhancer_gr),
                     emar_ends)) + 2000

plotTracks(
  list(axisTrack, grtrack, promoterTrack, enhancerTrack, emarTrack, hyper_track),
  from = from_region,
  to   = to_region,
  chromosome = chr_region,
  main = "IMPDH-209 (reverse strand, displayed 5' → 3')",
  reverseStrand = TRUE,    # <<< ESTE ES EL BUENO
  sizes = c(0.6, 1.0, 0.3, 0.3, 0.3, 0.5)
)

plotTracks(
  list(axisTrack, grtrack, promoterTrack, enhancerTrack, hyper_track, emarTrack),
  from = from_region,
  to   = to_region,
  chromosome = chr_region,
  main = "IMPDH-209 (reverse strand, displayed 5' → 3')",
  reverseStrand = TRUE,
  sizes = c(1, 0.4, 0.4, 0.4, 0.4, 0.4)   # todos los tracks con la MISMA altura
)

###akna Ovaries
###Gene representation####

library(biomaRt)

ensembl <- useEnsembl(biomart = "genes")

# Ver datasets disponibles
listDatasets(ensembl)



library(biomaRt)
library(Gviz)
library(GenomicRanges)

options(ucscChromosomeNames = FALSE)   # para que acepte CAJNNU010000024.1

ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "dlabrax_gene_ensembl"
)

atts <- listAttributes(ensembl)
atts[grep("transcript", atts$name), ][1:20, ]  # solo para curiosear


tx_id <- "ENSDLAT00005037644"   # sin .2

bm <- getBM(
  attributes = c(
    "chromosome_name",
    "exon_chrom_start", "exon_chrom_end",
    "strand",
    "ensembl_gene_id",
    "external_gene_name",
    "ensembl_transcript_id",
    "rank"
  ),
  filters  = "ensembl_transcript_id",
  values   = tx_id,
  mart     = ensembl
)

head(bm)
unique(bm$chromosome_name)




hypo <- read.table("FxF_akna_hipo.bed", sep = "\t", header = FALSE,
                   col.names = c("chr","start0","end","name","score","strand"))

# BED → 1-based
hypo$pos <- hypo$start0 + 1


# Asegurar que el nombre de cromosoma coincide con Ensembl
hypo$chr <- as.character(hyper$chr)
bm$chromosome_name <- as.character(bm$chromosome_name)

# Si por ejemplo bm tiene "CAJNN010000024.1" y hyper "CAJNN010000024.1", perfecto.
# Si fueran distintos, forzamos:
# hyper$chr <- unique(bm$chromosome_name)

chr_region <- unique(bm$chromosome_name)

## CpGs como DOTS (puntos redondos)
cpg_gr <- GRanges(
  seqnames = hypo$chr,
  ranges   = IRanges(start = hypo$pos, end = hypo$pos)
)

hypo_track <- DataTrack(
  range = cpg_gr,
  data  = rep(1, length(cpg_gr)),   # todos al mismo nivel
  genome = "dlabrax2021",
  chromosome = chr_region,
  name  = "Hypo CpGs",
  type  = "p",                      # <<< PUNTOS
  col   = "#8B4789",                # verde
  cex   = 1.2,                      # tamaño del punto (ajusta si quieres)
  ylim  = c(0.5, 1.5),
  yAxis = FALSE,
  background.title = "white"
)




gene_gr <- GRanges(
  seqnames = bm$chromosome_name,
  ranges   = IRanges(start = bm$exon_chrom_start,
                     end   = bm$exon_chrom_end),
  strand   = ifelse(bm$strand[1] == 1, "+", "-"),
  gene     = bm$external_gene_name,
  exon     = bm$rank,
  transcript = bm$ensembl_transcript_id
)

grtrack <- GeneRegionTrack(
  gene_gr,
  genome = "dlabrax2021",
  chromosome = chr_region,
  name = "akna-202",
  showId = TRUE,
  transcriptAnnotation = "transcript"
)

grtrack <- GeneRegionTrack(
  gene_gr,
  genome = "dlabrax2021",
  chromosome = chr_region,
  name = "akna-202",
  showId = FALSE,
  transcriptAnnotation = "transcript",
  stacking = "squish",      # <<< compacta los exones
  exonHeight = 1,        # <<< hace los exones FINOS
  arrowHeadWidth = 8,       # <<< flechas más pequeñas
  geneSymbol = TRUE,
  fill = "#EDD5B3",
  col  = "#EDD5B3",
  background.title = "white",
  cex.title = 0.9
)



axisTrack <- GenomeAxisTrack()

from_region <- min(bm$exon_chrom_start) - 1000
to_region   <- max(bm$exon_chrom_end)   + 1000

plotTracks(
  list(axisTrack, grtrack, hypo_track),
  from = from_region,
  to   = to_region,
  chromosome = chr_region,
  main = "akna-202 (canonical transcript)"
)


#####Añadiendo promotores y enhancers del gen
# Asegúrate de tener esto antes
# chr_region <- unique(bm$chromosome_name)

## PROMOTORES (2 regiones)
promoter_gr <- GRanges(
  seqnames = chr_region,
  ranges = IRanges(
    start = c(19300117),   
    end   = c(19300217)    
  )
)

## ENHANCERS (3 regiones)
enhancer_gr <- GRanges(
  seqnames = chr_region,
  ranges = IRanges(
    start = c(19299331, 19296500),  
    end   = c(19300116, 19296764)   
  )
)

promoterTrack <- AnnotationTrack(
  promoter_gr,
  genome = "dlabrax2021",
  chromosome = chr_region,
  name = "Promoters",
  fill = "darkred",
  col  = NA,
  alpha = 0.4      # transparente
)

enhancerTrack <- AnnotationTrack(
  enhancer_gr,
  genome = "dlabrax2021",
  chromosome = chr_region,
  name = "Enhancers",
  fill = "darkorange",
  col  = NA,
  alpha = 0.4
)

from_region <- min(
  c(bm$exon_chrom_start,
    start(promoter_gr),
    start(enhancer_gr))
) - 2000

to_region <- max(
  c(bm$exon_chrom_end,
    end(promoter_gr),
    end(enhancer_gr))
) + 2000



##axisTrack <- GenomeAxisTrack()

plotTracks(
  list(axisTrack, grtrack, promoterTrack, enhancerTrack, hypo_track),
  from = from_region,
  to   = to_region,
  chromosome = chr_region,
  main = "akna-202 (canonical transcript)",
  sizes = c(1, 3, 0.8, 0.8, 1)   # ajusta proporción de cada pista si quieres
)



##Añadiendo la EMAR que se solapa 
chr_region <- unique(bm$chromosome_name)
chr_region
# debería devolver "CAJNNU010000024.1"

library(Gviz)
library(GenomicRanges)

## EMARs (todas las de Ensembl en la región del gen)
emar_starts <- c(19299331L, 19296500L)
emar_ends   <- c(19300139L, 19296764L)

emar_gr <- GRanges(
  seqnames = chr_region,
  ranges   = IRanges(start = emar_starts, end = emar_ends)
)

emarTrack <- AnnotationTrack(
  range = emar_gr,
  genome = "dlabrax2021",
  chromosome = chr_region,
  name = "EMAR",
  fill = "darkgreen",
  col  = NA,
  alpha = 0.4
)

from_region <- min(c(bm$exon_chrom_start,
                     start(promoter_gr),
                     start(enhancer_gr),
                     emar_starts)) - 2000

to_region   <- max(c(bm$exon_chrom_end,
                     end(promoter_gr),
                     end(enhancer_gr),
                     emar_ends)) + 2000

## Open chromatin (todas las de Ensembl en la región)
open_starts <- c(19288809L, 19284448L, 19284019L, 19272735L, 19270021L)
open_ends   <- c(19289214L, 19285001L, 19284249L, 19272950L, 19270236L)

open_gr <- GRanges(
  seqnames = chr_region,
  ranges   = IRanges(start = open_starts, end = open_ends)
)

openTrack <- AnnotationTrack(
  range  = open_gr,
  genome = "dlabrax2021",
  chromosome = chr_region,
  name  = "Open chromatin",
  fill  = "darkgray",   # gris
  col   = NA,
  alpha = 0.5
)




plotTracks(
  list(axisTrack, grtrack, promoterTrack, enhancerTrack, emarTrack, hypo_track, openTrack),
  from = from_region,
  to   = to_region,
  chromosome = chr_region,
  main = "akna-202 (reverse strand, displayed 5' → 3')",
  reverseStrand = TRUE,    # <<< ESTE ES EL BUENO
  sizes = c(0.6, 1.0, 0.3, 0.3, 0.3, 0.5, 0.3)
)

plotTracks(
  list(axisTrack, grtrack, promoterTrack, enhancerTrack, hypo_track, emarTrack, openTrack),
  from = from_region,
  to   = to_region,
  chromosome = chr_region,
  main = "akna-202 (reverse strand, displayed 5' → 3')",
  reverseStrand = TRUE,
  sizes = c(1, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4)   # todos los tracks con la MISMA altura
)


###akna Testes
###Gene representation####

library(biomaRt)

ensembl <- useEnsembl(biomart = "genes")

# Ver datasets disponibles
listDatasets(ensembl)



library(biomaRt)
library(Gviz)
library(GenomicRanges)

options(ucscChromosomeNames = FALSE)   # para que acepte CAJNNU010000024.1

ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "dlabrax_gene_ensembl"
)

atts <- listAttributes(ensembl)
atts[grep("transcript", atts$name), ][1:20, ]  # solo para curiosear


tx_id <- "ENSDLAT00005037644"   # sin .2

bm <- getBM(
  attributes = c(
    "chromosome_name",
    "exon_chrom_start", "exon_chrom_end",
    "strand",
    "ensembl_gene_id",
    "external_gene_name",
    "ensembl_transcript_id",
    "rank"
  ),
  filters  = "ensembl_transcript_id",
  values   = tx_id,
  mart     = ensembl
)

head(bm)
unique(bm$chromosome_name)




hypo <- read.table("MxM_akna_hiper.bed", sep = "\t", header = FALSE,
                   col.names = c("chr","start0","end","name","score","strand"))

# BED → 1-based
hypo$pos <- hypo$start0 + 1


# Asegurar que el nombre de cromosoma coincide con Ensembl
hypo$chr <- as.character(hyper$chr)
bm$chromosome_name <- as.character(bm$chromosome_name)

# Si por ejemplo bm tiene "CAJNN010000024.1" y hyper "CAJNN010000024.1", perfecto.
# Si fueran distintos, forzamos:
# hyper$chr <- unique(bm$chromosome_name)

chr_region <- unique(bm$chromosome_name)

## CpGs como DOTS (puntos redondos)
cpg_gr <- GRanges(
  seqnames = hypo$chr,
  ranges   = IRanges(start = hypo$pos, end = hypo$pos)
)

hypo_track <- DataTrack(
  range = cpg_gr,
  data  = rep(1, length(cpg_gr)),   # todos al mismo nivel
  genome = "dlabrax2021",
  chromosome = chr_region,
  name  = "Hypo CpGs",
  type  = "p",                      # <<< PUNTOS
  col   = "#548B54",                # verde
  cex   = 1.2,                      # tamaño del punto (ajusta si quieres)
  ylim  = c(0.5, 1.5),
  yAxis = FALSE,
  background.title = "white"
)




gene_gr <- GRanges(
  seqnames = bm$chromosome_name,
  ranges   = IRanges(start = bm$exon_chrom_start,
                     end   = bm$exon_chrom_end),
  strand   = ifelse(bm$strand[1] == 1, "+", "-"),
  gene     = bm$external_gene_name,
  exon     = bm$rank,
  transcript = bm$ensembl_transcript_id
)

grtrack <- GeneRegionTrack(
  gene_gr,
  genome = "dlabrax2021",
  chromosome = chr_region,
  name = "akna-202",
  showId = TRUE,
  transcriptAnnotation = "transcript"
)

grtrack <- GeneRegionTrack(
  gene_gr,
  genome = "dlabrax2021",
  chromosome = chr_region,
  name = "akna-202",
  showId = FALSE,
  transcriptAnnotation = "transcript",
  stacking = "squish",      # <<< compacta los exones
  exonHeight = 1,        # <<< hace los exones FINOS
  arrowHeadWidth = 8,       # <<< flechas más pequeñas
  geneSymbol = TRUE,
  fill = "#EDD5B3",
  col  = "#EDD5B3",
  background.title = "white",
  cex.title = 0.9
)



axisTrack <- GenomeAxisTrack()

from_region <- min(bm$exon_chrom_start) - 1000
to_region   <- max(bm$exon_chrom_end)   + 1000

plotTracks(
  list(axisTrack, grtrack, hypo_track),
  from = from_region,
  to   = to_region,
  chromosome = chr_region,
  main = "akna-202 (canonical transcript)"
)


#####Añadiendo promotores y enhancers del gen
# Asegúrate de tener esto antes
# chr_region <- unique(bm$chromosome_name)

## PROMOTORES (2 regiones)
promoter_gr <- GRanges(
  seqnames = chr_region,
  ranges = IRanges(
    start = c(19300117),   
    end   = c(19300217)    
  )
)

## ENHANCERS (3 regiones)
enhancer_gr <- GRanges(
  seqnames = chr_region,
  ranges = IRanges(
    start = c(19299331, 19296500),  
    end   = c(19300116, 19296764)   
  )
)

promoterTrack <- AnnotationTrack(
  promoter_gr,
  genome = "dlabrax2021",
  chromosome = chr_region,
  name = "Promoters",
  fill = "darkred",
  col  = NA,
  alpha = 0.4      # transparente
)

enhancerTrack <- AnnotationTrack(
  enhancer_gr,
  genome = "dlabrax2021",
  chromosome = chr_region,
  name = "Enhancers",
  fill = "darkorange",
  col  = NA,
  alpha = 0.4
)

from_region <- min(
  c(bm$exon_chrom_start,
    start(promoter_gr),
    start(enhancer_gr))
) - 2000

to_region <- max(
  c(bm$exon_chrom_end,
    end(promoter_gr),
    end(enhancer_gr))
) + 2000



##axisTrack <- GenomeAxisTrack()

plotTracks(
  list(axisTrack, grtrack, promoterTrack, enhancerTrack, hypo_track),
  from = from_region,
  to   = to_region,
  chromosome = chr_region,
  main = "akna-202 (canonical transcript)",
  sizes = c(1, 3, 0.8, 0.8, 1)   # ajusta proporción de cada pista si quieres
)



##Añadiendo la EMAR que se solapa 
chr_region <- unique(bm$chromosome_name)
chr_region
# debería devolver "CAJNNU010000024.1"

library(Gviz)
library(GenomicRanges)

## EMARs (todas las de Ensembl en la región del gen)
emar_starts <- c(19299331L, 19296500L)
emar_ends   <- c(19300139L, 19296764L)

emar_gr <- GRanges(
  seqnames = chr_region,
  ranges   = IRanges(start = emar_starts, end = emar_ends)
)

emarTrack <- AnnotationTrack(
  range = emar_gr,
  genome = "dlabrax2021",
  chromosome = chr_region,
  name = "EMAR",
  fill = "darkgreen",
  col  = NA,
  alpha = 0.4
)

from_region <- min(c(bm$exon_chrom_start,
                     start(promoter_gr),
                     start(enhancer_gr),
                     emar_starts)) - 2000

to_region   <- max(c(bm$exon_chrom_end,
                     end(promoter_gr),
                     end(enhancer_gr),
                     emar_ends)) + 2000

## Open chromatin (todas las de Ensembl en la región)
open_starts <- c(19288809L, 19284448L, 19284019L, 19272735L, 19270021L)
open_ends   <- c(19289214L, 19285001L, 19284249L, 19272950L, 19270236L)

open_gr <- GRanges(
  seqnames = chr_region,
  ranges   = IRanges(start = open_starts, end = open_ends)
)

openTrack <- AnnotationTrack(
  range  = open_gr,
  genome = "dlabrax2021",
  chromosome = chr_region,
  name  = "Open chromatin",
  fill  = "darkgray",   # gris
  col   = NA,
  alpha = 0.5
)




plotTracks(
  list(axisTrack, grtrack, promoterTrack, enhancerTrack, emarTrack, hypo_track, openTrack),
  from = from_region,
  to   = to_region,
  chromosome = chr_region,
  main = "akna-202 (reverse strand, displayed 5' → 3')",
  reverseStrand = TRUE,    # <<< ESTE ES EL BUENO
  sizes = c(0.6, 1.0, 0.3, 0.3, 0.3, 0.5, 0.3)
)

plotTracks(
  list(axisTrack, grtrack, promoterTrack, enhancerTrack, hypo_track, emarTrack, openTrack),
  from = from_region,
  to   = to_region,
  chromosome = chr_region,
  main = "akna-202 (reverse strand, displayed 5' → 3')",
  reverseStrand = TRUE,
  sizes = c(1, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4)   # todos los tracks con la MISMA altura
)

##vash2 ovaries
####### CONFIG BÁSICA #######
library(biomaRt)
library(Gviz)
library(GenomicRanges)

options(ucscChromosomeNames = FALSE)

ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "dlabrax_gene_ensembl"
)

## IDs de las DOS isoformas (canónica primero)
tx_ids <- c("ENSDLAT00005081048",   # canónica
            "ENSDLAT00005068485")   # no canónica
gene_label <- "GENE-vash2"              # pon aquí el nombre del gen (AKNA, etc.)

####### EXONES DE LAS DOS ISOFORMAS #######
bm <- getBM(
  attributes = c(
    "chromosome_name",
    "exon_chrom_start", "exon_chrom_end",
    "strand",
    "ensembl_gene_id",
    "external_gene_name",
    "ensembl_transcript_id",
    "rank"
  ),
  filters  = "ensembl_transcript_id",
  values   = tx_ids,
  mart     = ensembl
)

chr_region <- unique(bm$chromosome_name)

gene_gr <- GRanges(
  seqnames = bm$chromosome_name,
  ranges   = IRanges(start = bm$exon_chrom_start,
                     end   = bm$exon_chrom_end),
  strand   = ifelse(bm$strand[1] == 1, "+", "-"),
  gene     = bm$external_gene_name,
  exon     = bm$rank,
  transcript = bm$ensembl_transcript_id
)

## Forzar orden: canónica arriba, otra abajo
gene_gr$transcript <- factor(
  gene_gr$transcript,
  levels = tx_ids
)

grtrack <- GeneRegionTrack(
  gene_gr,
  genome = "dlabrax2021",
  chromosome = chr_region,
  name = gene_label,
  showId = FALSE,
  transcriptAnnotation = "transcript",
  stacking = "squish",
  exonHeight = 0.4,
  arrowHeadWidth = 6,
  geneSymbol = TRUE,
  fill = "#EDD5B3",
  col  = "#EDD5B3",
  background.title = "white",
  cex.title = 0.9
)

####### CARGA DE CpGs HIPER (archivo x.tmmb.bed) #######
hyper <- read.table("FxF_vash2_hiper.bed", sep = "\t", header = FALSE,
                    col.names = c("chr","start0","end","name","score","strand"))
hyper$pos <- hyper$start0 + 1

hyper_gr <- GRanges(
  seqnames = hyper$chr,
  ranges   = IRanges(start = hyper$pos, end = hyper$pos)
)

hyper_track <- DataTrack(
  range = hyper_gr,
  data  = rep(1, length(hyper_gr)),
  genome = "dlabrax2021",
  chromosome = chr_region,
  name  = "Hyper CpGs",
  type  = "p",
  col   = "#548B54",      # verde (como IMPDH)
  cex   = 1.2,
  ylim  = c(0.5, 1.5),
  yAxis = FALSE,
  background.title = "white"
)

####### CARGA DE CpGs HIPO (archivo x2.tmmb.bed) #######
hypo <- read.table("FxF_vash2_hipo.bed", sep = "\t", header = FALSE,
                   col.names = c("chr","start0","end","name","score","strand"))
hypo$pos <- hypo$start0 + 1

hypo_gr <- GRanges(
  seqnames = hypo$chr,
  ranges   = IRanges(start = hypo$pos, end = hypo$pos)
)

hypo_track <- DataTrack(
  range = hypo_gr,
  data  = rep(1, length(hypo_gr)),
  genome = "dlabrax2021",
  chromosome = chr_region,
  name  = "Hypo CpGs",
  type  = "p",
  col   = "#8B4789",      # azul (como en el otro gen)
  cex   = 1.2,
  ylim  = c(0.5, 1.5),
  yAxis = FALSE,
  background.title = "white"
)

####### RANGO GENÓMICO (de momento solo exones) #######
from_region <- min(bm$exon_chrom_start) - 2000
to_region   <- max(bm$exon_chrom_end)   + 2000

axisTrack <- GenomeAxisTrack()

####### PLOT BÁSICO (solo gen + CpGs, ya 5'→3') #######
plotTracks(
  list(axisTrack,
       grtrack,
       hyper_track,
       hypo_track),
  from = from_region,
  to   = to_region,
  chromosome = chr_region,
  main = paste0(gene_label, " (reverse strand, displayed 5' → 3')"),
  sizes = c(0.6, 1.2, 0.6, 0.6)   # eje, gen, hyper, hypo
)

## PROMOTOR
promoter_gr <- GRanges(
  seqnames = chr_region,   # debería ser "CAJNNU010000008.1"
  ranges = IRanges(
    start = 20128143L,     # 20,128,143
    end   = 20128243L      # 20,128,243
  )
)

promoterTrack <- AnnotationTrack(
  promoter_gr,
  genome = "dlabrax2021",
  chromosome = chr_region,
  name = "Promoter",
  fill = "darkred",
  col  = NA,
  alpha = 0.4
)

plotTracks(
  list(axisTrack,
       grtrack,
       hyper_track,
       hypo_track,
       promoterTrack),
  from = from_region,
  to   = to_region,
  chromosome = chr_region,
  main = paste0(gene_label, " (reverse strand, displayed 5' → 3')"),
  sizes = c(0.6, 1.2, 0.6, 0.6, 0.6)   # eje, gen, hyper, hypo
)

## ENHANCERS
enhancer_starts <- c(
  20128643L,  # 20,128,643–20,129,784
  20119872L,  # 20,119,872–20,120,182
  20108814L,  # 20,108,814–20,109,155
  20108048L,  # 20,108,048–20,108,544
  20103391L,  # 20,103,391–20,103,826
  20102504L,  # 20,102,504–20,102,851
  20097108L,  # 20,097,108–20,097,515
  20086738L   # 20,086,738–20,087,721
)

enhancer_ends <- c(
  20129784L,  # 20,129,784
  20120182L,  # 20,120,182
  20109155L,  # 20,109,155
  20108544L,  # 20,108,544
  20103826L,  # 20,103,826
  20102851L,  # 20,102,851
  20097515L,  # 20,097,515
  20087721L   # 20,087,721
)

enhancer_gr <- GRanges(
  seqnames = chr_region,
  ranges   = IRanges(start = enhancer_starts, end = enhancer_ends)
)

enhancerTrack <- AnnotationTrack(
  enhancer_gr,
  genome = "dlabrax2021",
  chromosome = chr_region,
  name = "Enhancers",
  fill = "darkorange",
  col  = NA,
  alpha = 0.4
)

plotTracks(
  list(axisTrack,
       grtrack,
       hyper_track,
       hypo_track,
       promoterTrack,
       enhancerTrack),
  from = from_region,
  to   = to_region,
  chromosome = chr_region,
  main = paste0(gene_label, " (reverse strand, displayed 5' → 3')"),
  sizes = c(0.6, 1.2, 0.6, 0.6, 0.6, 0.6)   # eje, gen, hyper, hypo
)

###from_region <- min(c(
  ##bm$exon_chrom_start,
 ## start(promoter_gr),
  ##enhancer_starts
###)) - 2000

###to_region <- max(c(
 ## bm$exon_chrom_end,
 ## end(promoter_gr),
 ## enhancer_ends
###)) + 2000

open_starts <- c(
  20182099L,
  20148830L,
  20141779L,
  20137516L,
  20126424L,
  20125216L,
  20100103L,
  20099313L
)

open_ends <- c(
  20182387L,
  20149118L,
  20142138L,
  20138325L,
  20126624L,
  20125833L,
  20101139L,
  20099514L
)

open_gr <- GRanges(
  seqnames = chr_region,
  ranges   = IRanges(start = open_starts, end = open_ends)
)

openTrack <- AnnotationTrack(
  range = open_gr,
  genome = "dlabrax2021",
  chromosome = chr_region,
  name = "Open chromatin",
  fill = "darkgray",
  col  = NA,
  alpha = 0.5
)

###gene_range <- GRanges(
  ###seqnames = chr_region,
  ##ranges = IRanges(
   ## start = min(bm$exon_chrom_start),
    ##end   = max(bm$exon_chrom_end)
  )
)

##open_gr <- subsetByOverlaps(open_gr, gene_range)
##openTrack <- openTrack[open_gr]

plotTracks(
  list(axisTrack,
       grtrack,
       hyper_track,
       hypo_track,
       promoterTrack,
       enhancerTrack,
       openTrack),
  from = from_region,
  to   = to_region,
  chromosome = chr_region,
  main = paste0(gene_label, " (reverse strand, displayed 5' → 3')"),
  sizes = c(0.6, 1.2, 0.6, 0.6, 0.6, 0.6, 0.6)   # eje, gen, hyper, hypo
)


## EMAR (todas las de Ensembl en la región del gen)
emar_starts <- c(
  20128643L,
  20127886L,
  20119872L,
  20108814L,
  20108048L,
  20103391L,
  20102504L,
  20097108L,
  20086738L
)

emar_ends <- c(
  20129784L,
  20128202L,
  20120182L,
  20109155L,
  20108544L,
  20103826L,
  20102851L,
  20097515L,
  20087721L
)

## --- EMAR (filtradas al rango del gen) ---
emar_gr <- GRanges(
  seqnames = chr_region,
  ranges   = IRanges(start = emar_starts, end = emar_ends)
)

##gene_range <- GRanges(
#  seqnames = chr_region,
 # ranges = IRanges(
  #  start = min(bm$exon_chrom_start),
  #  end   = max(bm$exon_chrom_end)
#  )
#)

#emar_gr <- subsetByOverlaps(emar_gr, gene_range)

emarTrack <- AnnotationTrack(
  range = emar_gr,
  genome = "dlabrax2021",
  chromosome = chr_region,
  name = "EMAR",
  fill = "darkgreen",
  col  = NA,
  alpha = 0.5
)

plotTracks(
  list(axisTrack,
       grtrack,
       hyper_track,
       hypo_track,
       promoterTrack,
       enhancerTrack,
       openTrack,
       emarTrack),
  from = from_region,
  to   = to_region,
  chromosome = chr_region,
  main = paste0(gene_label, " (reverse strand, displayed 5' → 3')"),
  sizes = c(1, 1, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4)   # eje, gen, hyper, hypo
)

## --- Open chromatin (filtradas igual) ---
##open_gr <- GRanges(
#  seqnames = chr_region,
 # ranges   = IRanges(start = open_starts, end = open_ends)
#)

#open_gr <- subsetByOverlaps(open_gr, gene_range)

#openTrack <- AnnotationTrack(
  #range  = open_gr,
  #genome = "dlabrax2021",
 # chromosome = chr_region,
 # name  = "Open chromatin",
 # fill  = "darkgray",
  #col   = NA,
 # alpha = 0.5
#)

#axisTrack <- GenomeAxisTrack(
  #chromosome  = chr_region,
 # littleTicks = TRUE,
#  labelPos    = "alternating",
#  scale       = 1e-6      # convierte bp a Mb en el eje
#)
#displayPars(axisTrack) <- list(
# showAxis       = TRUE,
 # showTickLabels = TRUE,
 # fontsize       = 10
#)



# Forzar orden: primero la canónica, luego la otra
tx_ids <- c(
  "ENSDLAT00005081048",  # <- CANÓNICA (arriba)
  "ENSDLAT00005068485"   # <- no canónica (abajo)
)

gene_gr <- GRanges(
  seqnames = bm$chromosome_name,
  ranges   = IRanges(start = bm$exon_chrom_start,
                     end   = bm$exon_chrom_end),
  strand   = ifelse(bm$strand[1] == 1, "+", "-"),
  gene     = bm$external_gene_name,
  exon     = bm$rank,
  transcript = bm$ensembl_transcript_id
)

## Asegurar que las isoformas están en el orden deseado
gene_gr$transcript <- factor(
  gene_gr$transcript,
  levels = tx_ids           # canónica primero
)

## Reordenar las filas del GRanges según ese factor
gene_gr <- gene_gr[order(gene_gr$transcript), ]


grtrack <- GeneRegionTrack(
  gene_gr,
  genome = "dlabrax2021",
  chromosome = chr_region,
  name = gene_label,
  showId = FALSE,
  transcriptAnnotation = "transcript",
  stacking = "squish",
  exonHeight = 0.4,
  arrowHeadWidth = 6,
  geneSymbol = TRUE,
  fill = "#EDD5B3",
  col  = "#EDD5B3",
  background.title = "white",
  cex.title = 0.9
)



levels(gene_gr$transcript) ##comprobar cual isoforma nos da primero


plotTracks(
  list(axisTrack,
       grtrack,
       hyper_track,
       hypo_track,
       promoterTrack,
       enhancerTrack,
       openTrack,
       emarTrack),
  from = from_region,
  to   = to_region,
  chromosome = chr_region,
  main = paste0(gene_label, " (reverse strand, displayed 5' → 3')"),
  sizes = c(1, 1, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4)   # eje, gen, hyper, hypo
)




##vash2 testes
####### CONFIG BÁSICA #######
library(biomaRt)
library(Gviz)
library(GenomicRanges)

options(ucscChromosomeNames = FALSE)

ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "dlabrax_gene_ensembl"
)

## IDs de las DOS isoformas (canónica primero)
tx_ids <- c("ENSDLAT00005081048",   # canónica
            "ENSDLAT00005068485")   # no canónica
gene_label <- "GENE-vash2"              # pon aquí el nombre del gen (AKNA, etc.)

####### EXONES DE LAS DOS ISOFORMAS #######
bm <- getBM(
  attributes = c(
    "chromosome_name",
    "exon_chrom_start", "exon_chrom_end",
    "strand",
    "ensembl_gene_id",
    "external_gene_name",
    "ensembl_transcript_id",
    "rank"
  ),
  filters  = "ensembl_transcript_id",
  values   = tx_ids,
  mart     = ensembl
)

chr_region <- unique(bm$chromosome_name)

gene_gr <- GRanges(
  seqnames = bm$chromosome_name,
  ranges   = IRanges(start = bm$exon_chrom_start,
                     end   = bm$exon_chrom_end),
  strand   = ifelse(bm$strand[1] == 1, "+", "-"),
  gene     = bm$external_gene_name,
  exon     = bm$rank,
  transcript = bm$ensembl_transcript_id
)

## Forzar orden: canónica arriba, otra abajo
gene_gr$transcript <- factor(
  gene_gr$transcript,
  levels = tx_ids
)

grtrack <- GeneRegionTrack(
  gene_gr,
  genome = "dlabrax2021",
  chromosome = chr_region,
  name = gene_label,
  showId = FALSE,
  transcriptAnnotation = "transcript",
  stacking = "squish",
  exonHeight = 0.4,
  arrowHeadWidth = 6,
  geneSymbol = TRUE,
  fill = "#EDD5B3",
  col  = "#EDD5B3",
  background.title = "white",
  cex.title = 0.9
)

####### CARGA DE CpGs HIPER (archivo x.tmmb.bed) #######
hyper <- read.table("MxM_vash2_hiper.bed", sep = "\t", header = FALSE,
                    col.names = c("chr","start0","end","name","score","strand"))
hyper$pos <- hyper$start0 + 1

hyper_gr <- GRanges(
  seqnames = hyper$chr,
  ranges   = IRanges(start = hyper$pos, end = hyper$pos)
)

hyper_track <- DataTrack(
  range = hyper_gr,
  data  = rep(1, length(hyper_gr)),
  genome = "dlabrax2021",
  chromosome = chr_region,
  name  = "Hyper CpGs",
  type  = "p",
  col   = "#548B54",      # verde (como IMPDH)
  cex   = 1.2,
  ylim  = c(0.5, 1.5),
  yAxis = FALSE,
  background.title = "white"
)

####### CARGA DE CpGs HIPO (archivo x2.tmmb.bed) #######
hypo <- read.table("MxM_vash2_hipo.bed", sep = "\t", header = FALSE,
                   col.names = c("chr","start0","end","name","score","strand"))
hypo$pos <- hypo$start0 + 1

hypo_gr <- GRanges(
  seqnames = hypo$chr,
  ranges   = IRanges(start = hypo$pos, end = hypo$pos)
)

hypo_track <- DataTrack(
  range = hypo_gr,
  data  = rep(1, length(hypo_gr)),
  genome = "dlabrax2021",
  chromosome = chr_region,
  name  = "Hypo CpGs",
  type  = "p",
  col   = "#8B4789",      # azul (como en el otro gen)
  cex   = 1.2,
  ylim  = c(0.5, 1.5),
  yAxis = FALSE,
  background.title = "white"
)

####### RANGO GENÓMICO (de momento solo exones) #######
from_region <- min(bm$exon_chrom_start) - 2000
to_region   <- max(bm$exon_chrom_end)   + 2000

axisTrack <- GenomeAxisTrack()

####### PLOT BÁSICO (solo gen + CpGs, ya 5'→3') #######
plotTracks(
  list(axisTrack,
       grtrack,
       hyper_track,
       hypo_track),
  from = from_region,
  to   = to_region,
  chromosome = chr_region,
  main = paste0(gene_label, " (reverse strand, displayed 5' → 3')"),
  sizes = c(0.6, 1.2, 0.6, 0.6)   # eje, gen, hyper, hypo
)

## PROMOTOR
promoter_gr <- GRanges(
  seqnames = chr_region,   # debería ser "CAJNNU010000008.1"
  ranges = IRanges(
    start = 20128143L,     # 20,128,143
    end   = 20128243L      # 20,128,243
  )
)

promoterTrack <- AnnotationTrack(
  promoter_gr,
  genome = "dlabrax2021",
  chromosome = chr_region,
  name = "Promoter",
  fill = "darkred",
  col  = NA,
  alpha = 0.4
)

plotTracks(
  list(axisTrack,
       grtrack,
       hyper_track,
       hypo_track,
       promoterTrack),
  from = from_region,
  to   = to_region,
  chromosome = chr_region,
  main = paste0(gene_label, " (reverse strand, displayed 5' → 3')"),
  sizes = c(0.6, 1.2, 0.6, 0.6, 0.6)   # eje, gen, hyper, hypo
)

## ENHANCERS
enhancer_starts <- c(
  20128643L,  # 20,128,643–20,129,784
  20119872L,  # 20,119,872–20,120,182
  20108814L,  # 20,108,814–20,109,155
  20108048L,  # 20,108,048–20,108,544
  20103391L,  # 20,103,391–20,103,826
  20102504L,  # 20,102,504–20,102,851
  20097108L,  # 20,097,108–20,097,515
  20086738L   # 20,086,738–20,087,721
)

enhancer_ends <- c(
  20129784L,  # 20,129,784
  20120182L,  # 20,120,182
  20109155L,  # 20,109,155
  20108544L,  # 20,108,544
  20103826L,  # 20,103,826
  20102851L,  # 20,102,851
  20097515L,  # 20,097,515
  20087721L   # 20,087,721
)

enhancer_gr <- GRanges(
  seqnames = chr_region,
  ranges   = IRanges(start = enhancer_starts, end = enhancer_ends)
)

enhancerTrack <- AnnotationTrack(
  enhancer_gr,
  genome = "dlabrax2021",
  chromosome = chr_region,
  name = "Enhancers",
  fill = "darkorange",
  col  = NA,
  alpha = 0.4
)

plotTracks(
  list(axisTrack,
       grtrack,
       hyper_track,
       hypo_track,
       promoterTrack,
       enhancerTrack),
  from = from_region,
  to   = to_region,
  chromosome = chr_region,
  main = paste0(gene_label, " (reverse strand, displayed 5' → 3')"),
  sizes = c(0.6, 1.2, 0.6, 0.6, 0.6, 0.6)   # eje, gen, hyper, hypo
)

open_starts <- c(
  20182099L,
  20148830L,
  20141779L,
  20137516L,
  20126424L,
  20125216L,
  20100103L,
  20099313L
)

open_ends <- c(
  20182387L,
  20149118L,
  20142138L,
  20138325L,
  20126624L,
  20125833L,
  20101139L,
  20099514L
)

open_gr <- GRanges(
  seqnames = chr_region,
  ranges   = IRanges(start = open_starts, end = open_ends)
)

openTrack <- AnnotationTrack(
  range = open_gr,
  genome = "dlabrax2021",
  chromosome = chr_region,
  name = "Open chromatin",
  fill = "darkgray",
  col  = NA,
  alpha = 0.5
)

plotTracks(
  list(axisTrack,
       grtrack,
       hyper_track,
       hypo_track,
       promoterTrack,
       enhancerTrack,
       openTrack),
  from = from_region,
  to   = to_region,
  chromosome = chr_region,
  main = paste0(gene_label, " (reverse strand, displayed 5' → 3')"),
  sizes = c(0.6, 1.2, 0.6, 0.6, 0.6, 0.6, 0.6)   # eje, gen, hyper, hypo
)


## EMAR (todas las de Ensembl en la región del gen)
emar_starts <- c(
  20128643L,
  20127886L,
  20119872L,
  20108814L,
  20108048L,
  20103391L,
  20102504L,
  20097108L,
  20086738L
)

emar_ends <- c(
  20129784L,
  20128202L,
  20120182L,
  20109155L,
  20108544L,
  20103826L,
  20102851L,
  20097515L,
  20087721L
)

## --- EMAR (filtradas al rango del gen) ---
emar_gr <- GRanges(
  seqnames = chr_region,
  ranges   = IRanges(start = emar_starts, end = emar_ends)
)

emarTrack <- AnnotationTrack(
  range = emar_gr,
  genome = "dlabrax2021",
  chromosome = chr_region,
  name = "EMAR",
  fill = "darkgreen",
  col  = NA,
  alpha = 0.5
)

plotTracks(
  list(axisTrack,
       grtrack,
       hyper_track,
       hypo_track,
       promoterTrack,
       enhancerTrack,
       openTrack,
       emarTrack),
  from = from_region,
  to   = to_region,
  chromosome = chr_region,
  main = paste0(gene_label, " (reverse strand, displayed 5' → 3')"),
  sizes = c(1, 1, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4)   # eje, gen, hyper, hypo
)

# Forzar orden: primero la canónica, luego la otra
tx_ids <- c(
  "ENSDLAT00005081048",  # <- CANÓNICA (arriba)
  "ENSDLAT00005068485"   # <- no canónica (abajo)
)

gene_gr <- GRanges(
  seqnames = bm$chromosome_name,
  ranges   = IRanges(start = bm$exon_chrom_start,
                     end   = bm$exon_chrom_end),
  strand   = ifelse(bm$strand[1] == 1, "+", "-"),
  gene     = bm$external_gene_name,
  exon     = bm$rank,
  transcript = bm$ensembl_transcript_id
)

## Asegurar que las isoformas están en el orden deseado
gene_gr$transcript <- factor(
  gene_gr$transcript,
  levels = tx_ids           # canónica primero
)

## Reordenar las filas del GRanges según ese factor
gene_gr <- gene_gr[order(gene_gr$transcript), ]


grtrack <- GeneRegionTrack(
  gene_gr,
  genome = "dlabrax2021",
  chromosome = chr_region,
  name = gene_label,
  showId = FALSE,
  transcriptAnnotation = "transcript",
  stacking = "squish",
  exonHeight = 0.4,
  arrowHeadWidth = 6,
  geneSymbol = TRUE,
  fill = "#EDD5B3",
  col  = "#EDD5B3",
  background.title = "white",
  cex.title = 0.9
)



levels(gene_gr$transcript) ##comprobar cual isoforma nos da primero


plotTracks(
  list(axisTrack,
       grtrack,
       hyper_track,
       hypo_track,
       promoterTrack,
       enhancerTrack,
       openTrack,
       emarTrack),
  from = from_region,
  to   = to_region,
  chromosome = chr_region,
  main = paste0(gene_label, " (reverse strand, displayed 5' → 3')"),
  sizes = c(1, 1, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4)   # eje, gen, hyper, hypo
)



