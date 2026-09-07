#!/usr/bin/env Rscript
# =============================================================================
#  prefetch_drug_signatures.R
#
#  Builds a SIGNED drug-response resource for the decoupler app's
#  pathway / signature panel.
#
#  Source: Enrichr "Drug_Perturbations_from_GEO_up" / "_down" - curated GEO
#  drug-treatment differential-expression signatures. Genes up on treatment
#  -> mor +1, genes down -> mor -1. All GEO series / samples for one drug are
#  collapsed to ONE consensus signature (sign of the vote across experiments;
#  genes that disagree net to 0 and are dropped).
#
#  Output:  decoupler_cache/drugpert_geo_human.csv   (source, target, mor)
#  The resource "Drug Perturbations from GEO (signed)" then appears in the
#  pathway / signature dropdown. Re-run any time to refresh.
#
#  Interpreting the score (ULM, like CytoSig):
#    + : your contrast resembles this drug's effect
#    - : your contrast resembles the REVERSE of this drug (candidate reversal)
#
#  Run:  Rscript prefetch_drug_signatures.R
#  Needs: data.table  (+ internet to maayanlab.cloud)
# =============================================================================
suppressPackageStartupMessages({ library(data.table) })
options(timeout = 1200)

CACHE_DIR <- "decoupler_cache"
dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)
URL <- "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=%s"
msg <- function(...) cat(sprintf(...), "\n")

# set name looks like "<drug words> <DBxxxx|PubChemCID> <organism> GSExxxx sample n"
drug_of <- function(nm) {
  d <- sub("\\s+(DB[0-9]+|[0-9]{3,}|mouse|human|rat)\\b.*$", "", nm,
           ignore.case = TRUE, perl = TRUE)
  d <- gsub("&#?[a-z0-9]+;", "", d, ignore.case = TRUE)   # stray HTML entities
  d <- gsub("\\s+", " ", d)
  tolower(trimws(d))
}

read_lib <- function(lib, mor) {
  f <- tempfile(fileext = ".gmt")
  utils::download.file(sprintf(URL, lib), f, quiet = TRUE, mode = "wb")
  ln <- strsplit(readLines(f, warn = FALSE), "\t", fixed = TRUE)
  rbindlist(lapply(ln, function(x) {
    if (length(x) < 3L) return(NULL)
    tg <- toupper(trimws(x[-(1:2)])); tg <- tg[nzchar(tg)]
    dr <- drug_of(x[1])
    if (!nzchar(dr) || !length(tg)) return(NULL)
    data.table(source = dr, target = tg, mor = mor)
  }))
}

msg("downloading Drug_Perturbations_from_GEO_up / _down ...")
d <- rbind(read_lib("Drug_Perturbations_from_GEO_up",    1L),
           read_lib("Drug_Perturbations_from_GEO_down", -1L))
msg("  raw: %d edges, %d drugs", nrow(d), uniqueN(d$source))

# drop clone IDs / loci that map to nothing in a human symbol matrix
d <- d[!grepl("RIK$|^GM[0-9]+$|^LOC[0-9]+$|^[0-9]", target)]
# consensus per drug: net up/down vote across all GEO series for that drug
d <- d[, .(v = sum(mor)), by = .(source, target)][v != 0]
# cap each drug to its 400 strongest / most cross-series-consistent genes so a
# signature stays a signature (some drugs span many series -> huge raw unions)
d <- d[order(source, -abs(v))][, utils::head(.SD, 400L), by = source]
d[, mor := sign(v)]
d <- d[, if (.N >= 15L) .SD, by = source]      # keep drugs with a usable set

out <- file.path(CACHE_DIR, "drugpert_geo_human.csv")
fwrite(d[, .(source, target, mor)], out)
gpd <- d[, .N, by = source]
msg("wrote %s", out)
msg("  %d drugs, %d edges, %d..%d genes/drug (median %d)",
    nrow(gpd), nrow(d), min(gpd$N), max(gpd$N), as.integer(median(gpd$N)))
msg("  e.g. %s", paste(head(sort(unique(d$source)), 12), collapse = ", "))
