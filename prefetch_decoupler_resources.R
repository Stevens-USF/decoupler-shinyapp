#!/usr/bin/env Rscript
# =============================================================================
#  prefetch_decoupler_resources.R
#
#  Downloads every prior-knowledge resource the decoupler app can use and
#  writes it to decoupler_cache/ as a plain <id>_<organism>.csv so the app
#  works offline (OmniPath's server is frequently down).
#
#  Run once (repeat when you want to refresh):
#    Rscript prefetch_decoupler_resources.R              # human
#    Rscript prefetch_decoupler_resources.R human mouse  # both
#
#  Needs: decoupleR, OmnipathR  (BiocManager::install(c("decoupleR","OmnipathR")))
#  The app's VIPER / GSEA regulator-scoring options need two more packages that
#  decoupleR only Suggests (not installed automatically):
#    BiocManager::install(c("viper", "fgsea"))
#  Without them, selecting VIPER or GSEA silently falls back to ULM with a
#  warning ("there is no package called 'viper'/'fgsea'").
#
#  Resources fetched:
#    collectri        CollecTRI            TF regulon, curated, signed
#    dorothea         DoRothEA A-C         TF regulon, signed, confidence-tiered
#    progeny          PROGENy (top 500)    14 signalling-pathway responsive genes
#    cytosig          CytoSig (top 500)    cytokine-response signatures, signed
#    msigdb_hallmark  MSigDB Hallmark      50 curated gene sets
#    msigdb_reactome  MSigDB Reactome      pathway gene sets
#    msigdb_kegg      MSigDB KEGG          pathway gene sets
#    msigdb_gobp      MSigDB GO:BP         (large) biological-process gene sets
#
#  ChEA / ENCODE / custom .gmt files are read straight from geneset_cache/ by
#  the app - nothing to prefetch for those.
# =============================================================================

suppressPackageStartupMessages({ library(decoupleR); library(data.table) })
options(timeout = 1200)

args      <- commandArgs(trailingOnly = TRUE)
organisms <- if (length(args)) args else "human"
CACHE_DIR <- "decoupler_cache"
dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)

msg  <- function(...) cat(sprintf(...), "\n")
save_csv <- function(df, id, org) {
  df <- as.data.frame(df)
  f  <- file.path(CACHE_DIR, sprintf("%s_%s.csv", id, org))
  data.table::fwrite(df, f)
  msg("  [ok] %-18s %7d edges  %5d sources  -> %s",
      id, nrow(df), length(unique(df$source)), basename(f))
}
try_fetch <- function(id, org, expr) {
  msg("- %s (%s)", id, org)
  df <- tryCatch(eval.parent(substitute(expr)), error = function(e) {
    msg("  [FAIL] %s", conditionMessage(e)); NULL
  })
  if (!is.null(df) && nrow(df)) save_csv(df, id, org)
}

msigdb <- function(collection, org) {
  r <- as.data.frame(get_resource("MSigDB", organism = org))
  r <- unique(r[r$collection == collection, c("geneset", "genesymbol")])
  data.frame(source = r$geneset, target = toupper(r$genesymbol), mor = 1)
}

for (org in organisms) {
  msg("\n================  %s  ================", org)

  try_fetch("collectri", org,
    get_collectri(organism = org, split_complexes = FALSE))

  try_fetch("dorothea", org, {
    d <- as.data.frame(get_dorothea(organism = org, levels = c("A", "B", "C")))
    d[, c("source", "target", "mor", "confidence")]
  })

  try_fetch("progeny", org, {
    d <- as.data.frame(get_progeny(organism = org, top = 500))
    data.frame(source = d$source, target = d$target, mor = d$weight)
  })

  try_fetch("cytosig", org, {
    r <- as.data.frame(get_resource("CytoSig", organism = org))
    r <- data.frame(source = r$cytokine_genesymbol,
                    target = r$target_genesymbol, mor = r$score)
    r <- r[is.finite(r$mor) & r$mor != 0 & r$source != "" & r$target != "", ]
    setDT(r)
    r <- r[order(-abs(mor))][!duplicated(paste(source, target))]   # unique edges
    as.data.frame(r[, head(.SD, 500L), by = source])
  })

  try_fetch("msigdb_hallmark", org, msigdb("hallmark", org))
  try_fetch("msigdb_reactome", org, msigdb("reactome_pathways", org))
  try_fetch("msigdb_kegg",     org, msigdb("kegg_pathways", org))
  try_fetch("msigdb_gobp",     org, msigdb("go_biological_process", org))
}

msg("\nDone. Files in %s/", CACHE_DIR)
