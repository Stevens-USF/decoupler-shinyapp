# =============================================================================
#  decoupler_module.R  --  Shiny module: upstream TF + pathway activity
#  ("IPA Upstream Regulator / Canonical Pathways" alternative, free stack)
#
#  Regulator (TF panel) sets, user-selectable:
#    CollecTRI            curated TF->target, signed          (get_collectri)
#    DoRothEA A-C         TF->target, signed, confidence      (get_dorothea)
#    ChEA 2022 .gmt       ChIP-based TF->target, unsigned     (geneset_cache/*.gmt)
#  Pathway / signature sets, user-selectable:
#    PROGENy              14 signalling-pathway responsive genes
#    CytoSig             ~40 cytokine-response signatures, signed
#    Drug Perturbations from GEO   signed drug up/down signatures, connectivity
#                                  (run prefetch_drug_signatures.R to build it)
#    MSigDB Hallmark / Reactome / KEGG / GO:BP   curated gene sets
#
#  Resources load from decoupler_cache/<id>_<organism>.csv when present
#  (run prefetch_decoupler_resources.R once), else live OmniPath, else fail
#  with a message. Add more .gmt files to geneset_cache/ and they appear
#  automatically in the regulator dropdown.
#
#  Input: MS-DAP differential-abundance results in LONG format - a data.frame
#  with columns  gene, contrast,  and at least one of
#  foldchange.log2 / pvalue / qvalue / effectsize.
#  (dea_app.R turns differential_abundance_analysis.xlsx into that shape;
#   in MS-DAP_visualizer pass store$dea after renaming its columns.)
#
#  NOTE these regulons were trained on transcriptomics. They work on
#  proteomics (target-protein levels still track regulator activity) - but
#  ALWAYS sanity-check on a positive control (e.g. a TNFa / LPS contrast
#  should call NFkB / TNF / interferon pathways active) before trusting the
#  novel calls.
#
#  Tabs: regulator activity (heatmap / table / top movers), pathway activity,
#  regulon inspector, "Contrast comparison" (signature concordance = score each
#  contrast against one contrast's signed gene signature; a differ-most regulator
#  table; and descriptive correlation visuals), and "Help / interpretation".
#
#  Integration recipe (MS-DAP_visualizer_FINAL.R):
#    library(decoupleR); source("decoupler_module.R")
#    menuItem("Upstream activity", tabName="decoupler", icon=icon("diagram-project"))
#    tabItem("decoupler", decouplerTabUI("dc"))
#    decouplerServer("dc", dea = reactive(store$dea))
#  store$dea already has columns gene / contrast / log2fc / pvalue / qvalue /
#  effectsize / algorithm - the module aliases log2fc and shows a DAA-algorithm
#  picker itself. Needs decoupler_cache/ and geneset_cache/ in the working dir
#  (run prefetch_decoupler_resources.R once).
# =============================================================================

suppressPackageStartupMessages({
  library(shiny); library(shinydashboard); library(shinyWidgets)
  library(DT); library(plotly); library(ggplot2); library(data.table)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# where the bundled regulon files live (ship decoupler_cache/ next to this file)
CACHE_DIR   <- getOption("decoupler.cache",    "decoupler_cache")
# where plain Enrichr-style .gmt gene sets live (ChEA, ENCODE, custom signatures)
GENESET_DIR <- getOption("decoupler.genesets", "geneset_cache")
# regulons are always scored in human symbol space. MS-DAP emits upper-cased
# gene symbols; for a 1:1 mouse/human ortholog the symbol is identical once
# upper-cased, so human regulons match the large majority of a mouse dataset
# directly (on the PBIT mouse in-vitro set: 94% of genes hit ChEA 2022, only
# 5.7% hit no human regulon at all - and that residual is mostly RIKEN clone
# IDs and rodent-specific paralogs like Abcb1a/Abcb1b that have no 1:1 human
# ortholog, i.e. formal ortholog mapping would not recover them either).
# Override only if you have pre-mapped: options(decoupler.organism = "mouse")
DC_ORGANISM <- getOption("decoupler.organism", "human")

# ---- resource registry ------------------------------------------------------
# Every entry resolves to a data.frame(source, target, mor):
#   mor > 0  activating edge, mor < 0  repressing edge, mor == 1  unsigned set.
# panel = "tf"      -> shown in the regulator (TF) panels
# panel = "pathway" -> shown in the pathway / signature panel
# `fetch(organism)` is the live OmniPath/decoupleR call; used only when neither
# a bundled CSV nor an .rds cache is present (OmniPath is frequently down).
.DC_RESOURCES <- list(
  collectri = list(label = "CollecTRI  (TF, curated, signed)", panel = "tf",
    fetch = function(org)
      decoupleR::get_collectri(organism = org, split_complexes = FALSE)),
  # DoRothEA is disabled in the deployed build - it has no bundled CSV and
  # OmniPath is not reachable at run time. To re-enable: run
  # prefetch_decoupler_resources.R to create decoupler_cache/dorothea_human.csv,
  # then un-comment this entry.
  # dorothea = list(label = "DoRothEA A-C  (TF, signed)", panel = "tf",
  #   fetch = function(org) {
  #     d <- as.data.frame(decoupleR::get_dorothea(organism = org,
  #                                                levels = c("A", "B", "C")))
  #     d[, c("source", "target", "mor", "confidence")]
  #   }),
  progeny = list(label = "PROGENy  (14 signalling pathways)", panel = "pathway",
    fetch = function(org) {
      d <- as.data.frame(decoupleR::get_progeny(organism = org, top = 500))
      data.frame(source = d$source, target = d$target, mor = d$weight)
    }),
  drugpert_geo = list(
    label = "Drug Perturbations from GEO  (signed, up/down)", panel = "pathway",
    fetch = function(org) stop(
      "Build decoupler_cache/drugpert_geo_human.csv first: Rscript prefetch_drug_signatures.R",
      call. = FALSE)),
  cytosig = list(label = "CytoSig  (cytokine-response signatures)", panel = "pathway",
    fetch = function(org) {
      r <- as.data.frame(decoupleR::get_resource("CytoSig", organism = org))
      r <- data.frame(source = r$cytokine_genesymbol,
                      target = r$target_genesymbol, mor = r$score)
      r <- r[is.finite(r$mor) & r$mor != 0 & r$source != "" & r$target != "", ]
      data.table::setDT(r)
      as.data.frame(r[order(-abs(mor))][, utils::head(.SD, 500L), by = source])
    }),
  msigdb_hallmark = list(label = "MSigDB Hallmark", panel = "pathway",
    fetch = function(org) .dc_msigdb("hallmark", org)),
  msigdb_reactome = list(label = "MSigDB Reactome", panel = "pathway",
    fetch = function(org) .dc_msigdb("reactome_pathways", org)),
  msigdb_kegg = list(label = "MSigDB KEGG", panel = "pathway",
    fetch = function(org) .dc_msigdb("kegg_pathways", org)),
  msigdb_gobp = list(label = "MSigDB GO:BP  (large, slower)", panel = "pathway",
    fetch = function(org) .dc_msigdb("go_biological_process", org))
)

# discovered .gmt files in GENESET_DIR become resources "gmt:<filename>"
.dc_gmt_resources <- function() {
  if (!dir.exists(GENESET_DIR)) return(list())
  fs <- list.files(GENESET_DIR, pattern = "\\.gmt$", full.names = FALSE)
  stats::setNames(lapply(fs, function(f)
    list(label = sprintf("%s  (.gmt, unsigned)", tools::file_path_sans_ext(f)),
         panel = "tf", gmt = f)), paste0("gmt:", fs))
}

# all resources, keyed by id
.dc_all_resources <- function() c(.DC_RESOURCES, .dc_gmt_resources())

.dc_msigdb <- function(collection, organism = "human") {
  r <- as.data.frame(decoupleR::get_resource("MSigDB", organism = organism))
  r <- unique(r[r$collection == collection, c("geneset", "genesymbol")])
  data.frame(source = r$geneset, target = toupper(r$genesymbol), mor = 1)
}

# Enrichr-style .gmt -> data.frame(source, target, mor = 1). The set name's
# first whitespace token is taken as the regulator symbol (ChEA/ENCODE names
# look like "STAT1 HeLa-S3 hg19"); rows are de-duplicated so repeated
# experiments for one TF widen its target set rather than double-count.
.dc_read_gmt <- function(path) {
  ln <- strsplit(readLines(path, warn = FALSE), "\t", fixed = TRUE)
  ln <- ln[lengths(ln) >= 3L]
  out <- data.table::rbindlist(lapply(ln, function(x) {
    tg <- toupper(trimws(x[-(1:2)]))
    tg <- tg[tg != ""]
    if (!length(tg)) return(NULL)
    data.table::data.table(source = sub("\\s.*$", "", trimws(x[1])),
                           target = tg, mor = 1)
  }))
  unique(out)
}

# decoupleR rejects a network with repeated source-target edges; collapse them
# to the signed weight of largest magnitude.
.dc_dedup_edges <- function(net) {
  if (is.null(net) || !nrow(net)) return(net)
  d <- as.data.table(net)
  d <- d[!is.na(source) & !is.na(target) & source != "" & target != ""]
  if (!"mor" %in% names(d)) return(as.data.frame(unique(d)))
  d <- d[order(-abs(mor))][!duplicated(paste(source, target))]
  as.data.frame(d)
}

# Resolve a resource id to data.frame(source, target, mor).
# Order: bundled CSV -> .rds cache -> live fetch (then cached as .rds) -> NULL.
.dc_net <- function(kind, organism = "human") {
  res <- .dc_all_resources()[[kind]]
  if (is.null(res)) return(NULL)

  if (!is.null(res$gmt)) {
    f <- file.path(GENESET_DIR, res$gmt)
    return(if (file.exists(f)) .dc_dedup_edges(.dc_read_gmt(f)) else NULL)
  }

  dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)
  local <- file.path(CACHE_DIR, sprintf("%s_%s.csv", kind, organism))
  if (file.exists(local)) {
    net <- data.table::fread(local)
    if (!"mor" %in% names(net) && "weight" %in% names(net))
      data.table::setnames(net, "weight", "mor")
    return(.dc_dedup_edges(net))
  }
  rds <- file.path(CACHE_DIR, sprintf("%s_%s.rds", kind, organism))
  if (file.exists(rds)) return(.dc_dedup_edges(readRDS(rds)))

  net <- tryCatch(as.data.frame(res$fetch(organism)), error = function(e) NULL)
  if (!is.null(net) && nrow(net)) saveRDS(net, rds) else net <- NULL
  .dc_dedup_edges(net)
}

# LONG dea df -> gene x contrast matrix of the chosen signed statistic
.dc_matrix <- function(dea, stat, contrasts) {
  d <- as.data.table(dea)
  d <- d[contrast %in% contrasts & !is.na(gene) & gene != ""]
  d[, gene := toupper(sub(";.*$", "", gene))]           # first symbol, upper
  d[, value := switch(stat,
      "signed_-log10(p)"  = sign(foldchange.log2) * -log10(pmax(pvalue, 1e-300)),
      "signed_-log10(q)"  = sign(foldchange.log2) * -log10(pmax(qvalue, 1e-300)),
      "effectsize"        = effectsize,
      "foldchange.log2"   = foldchange.log2)]
  d <- d[is.finite(value)]
  # collapse duplicate gene symbols per contrast: keep the largest |value|
  d <- d[d[, .I[which.max(abs(value))], by = .(gene, contrast)]$V1]
  m <- data.table::dcast(d, gene ~ contrast, value.var = "value")
  mat <- as.matrix(m[, -1]); rownames(mat) <- m$gene
  mat[is.na(mat)] <- 0
  mat
}

# per-gene log2FC / p-value for ONE contrast, gene symbols upper-cased and
# de-duplicated the same way .dc_matrix does (keep the largest |log2FC|). Used
# by the inspector volcano; columns absent in the DAA table come back all-NA.
.dc_gene_stats <- function(dea, cc) {
  d <- as.data.frame(dea)
  d <- d[!is.na(d$contrast) & d$contrast == cc & !is.na(d$gene) & d$gene != "", ,
         drop = FALSE]
  g  <- toupper(sub(";.*$", "", d$gene))
  fc <- if ("foldchange.log2" %in% names(d)) d$foldchange.log2 else rep(NA_real_, nrow(d))
  pv <- if ("pvalue" %in% names(d)) d$pvalue else rep(NA_real_, nrow(d))
  o  <- order(-abs(ifelse(is.na(fc), 0, fc)))
  out <- data.table::data.table(gene = g[o], logfc = fc[o], pval = pv[o])
  out[!duplicated(out$gene)]
}

.dc_available_stats <- function(dea) {
  cn <- names(dea)
  c(if (all(c("foldchange.log2","pvalue") %in% cn)) "signed_-log10(p)",
    if (all(c("foldchange.log2","qvalue") %in% cn)) "signed_-log10(q)",
    if ("effectsize" %in% cn) "effectsize",
    if ("foldchange.log2" %in% cn) "foldchange.log2")
}

# ---- signed vs unsigned resources -----------------------------------------
# Signed resources carry a real mode of regulation (mor = +/-1 or a signed
# weight), so the score sign is an activation call. Everything else (ChEA,
# MSigDB) is mor = 1: the score sign is only "did this gene set move up or
# down as a group". Labels and the pathway scoring method switch on this.
.DC_SIGNED    <- c("collectri", "dorothea", "progeny", "cytosig", "drugpert_geo")
.dc_is_signed <- function(id) isTRUE(id %in% .DC_SIGNED)
.dc_score_lab <- function(signed) if (isTRUE(signed)) "activity" else "enrichment"
.dc_score_note <- function(signed) if (isTRUE(signed)) NULL else paste(
  "Unsigned set: the sign shows whether the set's genes moved up or down",
  "together, not whether the regulator/pathway is activated. See the Help tab.")

# tidy MSigDB set names for display: drop the collection prefix, "_" -> space,
# lower-case.  REACTOME_INTERFERON_SIGNALING -> "interferon signaling"
.dc_pretty <- function(x) {
  x <- sub("^(HALLMARK|KEGG_MEDICUS|KEGG|REACTOME|GOBP|GOCC|GOMF|WP|PID|BIOCARTA|NABA|MODULE)_",
           "", x, perl = TRUE)
  tolower(gsub("_", " ", x))
}

# Score `mat` against `net`. ULM runs as asked. MLM / consensus are tried, but
# unsigned collections (MSigDB Reactome/KEGG/GO especially) have heavily
# collinear, overlapping sources that a linear model cannot fit - on that
# error fall back to ULM and warn instead of failing the whole run.
.dc_score <- function(mat, net, method, minsize, warn = function(m) NULL) {
  run1 <- function(meth) switch(meth,
    ulm = decoupleR::run_ulm(mat, net, .source = "source", .target = "target",
                             .mor = "mor", minsize = minsize),
    mlm = decoupleR::run_mlm(mat, net, .source = "source", .target = "target",
                             .mor = "mor", minsize = minsize),
    consensus = {
      r <- decoupleR::decouple(mat, net, .source = "source", .target = "target",
             statistics = c("ulm", "mlm", "wsum"), consensus = TRUE,
             args = list(ulm = list(minsize = minsize),
                         mlm = list(minsize = minsize),
                         wsum = list(minsize = minsize)))
      r[r$statistic == "consensus", ]
    })
  if (identical(method, "ulm")) return(run1("ulm"))
  tryCatch(run1(method), error = function(e) {
    warn(sprintf("%s could not fit this set (collinear sources); scored with ULM instead.",
                 toupper(method)))
    run1("ulm")
  })
}

# ---- figure export -------------------------------------------------------------
# register PNG (300 dpi raster) + PDF (vector) download handlers for a reactive
# that returns a ggplot. ids become "<id>_png" / "<id>_pdf".
.dc_reg_dl <- function(output, id, plot_fun, base, w = 10, h = 7) {
  output[[paste0(id, "_png")]] <- shiny::downloadHandler(
    filename = function() sprintf("%s_%s.png", base, Sys.Date()),
    content  = function(f) ggplot2::ggsave(f, plot_fun(), width = w, height = h,
                                           dpi = 300, bg = "white"))
  output[[paste0(id, "_pdf")]] <- shiny::downloadHandler(
    filename = function() sprintf("%s_%s.pdf", base, Sys.Date()),
    content  = function(f) ggplot2::ggsave(f, plot_fun(), width = w, height = h,
                                           device = "pdf"))
}
# the paired PNG / PDF buttons for the UI
.dc_dl_ui <- function(ns, id) div(style = "margin:6px 0",
  downloadButton(ns(paste0(id, "_png")), "PNG (300 dpi)", class = "btn-xs"),
  downloadButton(ns(paste0(id, "_pdf")), "PDF (vector)",  class = "btn-xs"))

# heatmap helpers: shared matrix prep + a ggplot tile version for export
.dc_heat_mat <- function(dt, topn) {
  w <- data.table::dcast(dt, source ~ condition, value.var = "score")
  m <- as.matrix(w[, -1, drop = FALSE]); rownames(m) <- w$source
  if (nrow(m) > topn) {
    # rank rows for the "top N" cut: across-contrast SD when there are >=2
    # contrasts, else |score| (SD of one value is NA -> would fall back to the
    # dcast row order, i.e. alphabetical, and drop the real top hits)
    v <- if (ncol(m) >= 2) apply(m, 1, function(x) stats::sd(x, na.rm = TRUE))
         else abs(m[, 1])
    m <- m[order(-v)[seq_len(topn)], , drop = FALSE]
  }
  m[order(rowMeans(m, na.rm = TRUE)), , drop = FALSE]
}
.dc_heat_gg <- function(m, fill_lab) {
  d <- as.data.frame(as.table(as.matrix(m)))
  names(d) <- c("row", "col", "z")
  ggplot2::ggplot(d, ggplot2::aes(col, factor(row, levels = rownames(m)), fill = z)) +
    ggplot2::geom_tile(color = "grey92") +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                  midpoint = 0, name = fill_lab) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1))
}

# ---------------------------------------------------------------- help text --
.DC_HELP_HTML <- '
<div style="max-width:820px;padding:4px 8px;line-height:1.5">
<h3>Upstream regulator &amp; pathway activity &mdash; help</h3>

<p><b>What this does.</b> For each regulator (TF) or pathway it fits a linear
model of your per-gene differential statistic against that target set and reports
the model t-value as an <i>activity score</i>. This is activity <i>inference</i>
(it uses every measured gene), not over-representation of a cut list.</p>

<h4>Workflow</h4>
<ol>
<li>Pick the DAA algorithm (if several) and the contrasts to score.</li>
<li>Pick the per-gene statistic (see below).</li>
<li>Pick a regulator set, a pathway set, and the scoring method; Run.</li>
<li><b>Validate on a positive control first.</b> A TNF-alpha / LPS contrast must
call NF-kB / TNF / interferon active. If it does not, do not trust the rest.</li>
</ol>

<h4>Choosing the regulator / pathway set</h4>
<table border="1" cellpadding="4" style="border-collapse:collapse">
<tr><th>Set</th><th>Signed?</th><th>Use for</th></tr>
<tr><td>CollecTRI, DoRothEA A-C</td><td>yes (mor = +/-1)</td>
<td>TF activation / inhibition calls</td></tr>
<tr><td>PROGENy, CytoSig</td><td>yes (weights)</td>
<td>signalling / cytokine response</td></tr>
<tr><td>Drug Perturbations from GEO</td><td>yes (mor = +/-1)</td>
<td>connectivity: + your contrast resembles the drug, - resembles its reverse
(needs <code>prefetch_drug_signatures.R</code>)</td></tr>
<tr><td>ChEA 2022 (.gmt), MSigDB collections</td><td><b>no</b> (mor = 1)</td>
<td>direction of the target set only &mdash; NOT an activation call</td></tr>
</table>

<p><b>Unsigned sets.</b> Every edge weight is 1, so a positive score only means
"this regulator target set moved up as a group". It says nothing about whether
the regulator activates or represses those genes. Example: KDM5B is a repressor,
so if its targets go up its activity has most likely gone <i>down</i> &mdash; the
opposite of the sign shown. For activation calls use CollecTRI / DoRothEA; for
unsigned sets prefer GOAT / GSEA and supply the direction from biology yourself.</p>

<h4>Per-gene statistic</h4>
<ul>
<li><b>signed -log10(p)</b> &mdash; magnitude and precision. Standard default;
can favour well-measured (abundant) proteins.</li>
<li><b>effectsize</b> &mdash; magnitude scaled by variability; least
abundance-biased; most stable input for MLM.</li>
<li><b>foldchange.log2</b> &mdash; raw magnitude; sensitive to noisy large
changes; avoid with MLM.</li>
<li><b>signed -log10(q)</b> &mdash; thresholded and conservative; result depends
on how many proteins pass FDR in that contrast.</li>
</ul>
<p>Run with two of these and trust only regulators that survive both. Match the
choice to any GOAT / GSEA comparison (GOAT ranked by signed effect size is close
to <i>effectsize</i> here).</p>

<h4>Method</h4>
<ul>
<li><b>ULM</b> &mdash; one regulator at a time; fast; default.</li>
<li><b>MLM</b> &mdash; all regulators jointly; deconvolves overlapping target
sets (useful for redundant sets such as ChEA); more outlier-sensitive.</li>
<li><b>consensus</b> &mdash; mean of ULM / MLM / wsum; slow, robust.</li>
</ul>

<h4>Reading the score and the p-value</h4>
<p>The score is a t-value: sign = direction, magnitude = strength and
consistency. The <code>p_value</code> column is the p-value <i>of that t</i>, so
within one run score and p move together &mdash; a large |score| always has a
small p. It is <b>not</b> a hypergeometric / Fisher overlap p-value (IPA reports
those separately; this module does not), and it is <b>not FDR-corrected</b>.
Apply your own correction across regulators, and use the <b>Regulon inspector</b>
to confirm a call is not driven by two or three targets.</p>

<h4>Inspector tabs</h4>
<ul>
<li><b>Regulon inspector</b> &mdash; the chosen regulator&rsquo;s targets on this
contrast&rsquo;s volcano (log2FC vs -log10 p), as in the decoupleR TF vignette.
For a <i>signed</i> set (CollecTRI, DoRothEA) red = the gene moved the way its
edge predicts (activating target up, repressing target down) and so
<i>supports</i> an active call; blue opposes it. For an <i>unsigned</i> set
(ChEA, .gmt) every edge weight is 1, so colour is only the target&rsquo;s
direction &mdash; there is no activation call to support. With no p-value column
it falls back to a target bar.</li>
<li><b>Pathway inspector</b> &mdash; the same idea for the pathway / signature
set. PROGENy and CytoSig (continuous weights) get the decoupleR pathway-vignette
scatter of set weight vs your per-gene statistic (the &ldquo;MAPK&rdquo; view) &mdash;
top-right and bottom-left genes drive a positive score. Drug Perturbations from
GEO (&plusmn;1 up/down membership) gets the volcano instead, red where a gene
moved the way the signature predicts (your contrast resembles the drug). Unsigned
sets (MSigDB, ChEA) show a ranked target bar.</li>
</ul>

<h4>Contrast comparison tab</h4>
<ol>
<li><b>Signature concordance</b> &mdash; the quantitative comparison. Pick a
reference contrast; every contrast is scored (ULM) against the signed gene
signature of the reference (its top-N genes by |signed stat|). <b>+</b> =
resembles the reference, <b>-</b> = reverses it; the reference scores highest
against itself. Directional, and <i>not</i> inflated when contrasts share a
control group (which a raw correlation is). p-values are approximate &mdash;
rank by the score.</li>
<li><b>Regulators that differ most</b> &mdash; table of every regulator ranked by
|score_X - score_Y|; the condition-specific activity.</li>
<li><b>Exploratory visuals</b> &mdash; the contrast-by-contrast Spearman matrix
and the X-vs-Y scatter. Descriptive only: a positive correlation between two
contrasts that share a group is partly structural, so read a strong
<i>negative</i> correlation or an off-diagonal regulator, not the coefficient
itself.</li>
</ol>

<h4>Common pitfalls</h4>
<ul>
<li>Reading an unsigned-set sign as an activation call.</li>
<li>Trusting raw p-values without FDR.</li>
<li>These regulons were trained on transcriptomics &mdash; always sanity-check a
positive control on your proteomics data.</li>
<li>A regulon carried by a few targets &mdash; check the inspector.</li>
<li>Comparing decoupleR and GOAT run on different per-gene statistics.</li>
</ul>
</div>'

# ------------------------------------------------------------------- UI ------
decouplerTabUI <- function(id) {
  ns <- NS(id)
  fluidRow(
    box(width = 3, title = "Settings", status = "primary", solidHeader = TRUE,
      uiOutput(ns("ui_algo")),
      uiOutput(ns("ui_contrasts")),
      uiOutput(ns("ui_stat")),
      uiOutput(ns("ui_tf_resource")),
      uiOutput(ns("ui_pw_resource")),
      radioButtons(ns("method"), "Regulator scoring method",
                   c("ULM (fast, default)" = "ulm", "MLM" = "mlm",
                     "consensus (slow)" = "consensus"), selected = "ulm"),
      numericInput(ns("minsize"), "Min. targets per regulator", 5, 3, 50),
      numericInput(ns("topn"), "Rows shown in heatmaps (most variable)", 30, 10, 100),
      actionButton(ns("run"), "Run activity inference", icon = icon("play"),
                   class = "btn-primary", width = "100%"),
      uiOutput(ns("ui_status"))),
    tabBox(width = 9, id = ns("tabs"),
      tabPanel("Regulator activity - heatmap",
        uiOutput(ns("tf_scorenote")),
        plotlyOutput(ns("tf_heat"), height = "600px"),
        .dc_dl_ui(ns, "tf_heat")),
      tabPanel("Regulator activity - table",
        uiOutput(ns("ui_tf_contrast")),
        DTOutput(ns("tf_tbl")),
        downloadButton(ns("dl_tf"), "Download all regulator activities (.csv)")),
      tabPanel("Regulator activity - top movers",
        uiOutput(ns("ui_tf_contrast2")),
        plotOutput(ns("tf_bar"), height = "560px"),
        .dc_dl_ui(ns, "tf_bar")),
      tabPanel("Pathway / signature activity",
        uiOutput(ns("ui_pw_title")),
        uiOutput(ns("pw_scorenote")),
        plotlyOutput(ns("pw_heat"), height = "420px"),
        .dc_dl_ui(ns, "pw_heat"),
        DTOutput(ns("pw_tbl")),
        downloadButton(ns("dl_pw"), "Download pathway / signature activities (.csv)")),
      tabPanel("Regulon inspector",
        fluidRow(column(6, uiOutput(ns("ui_insp_tf"))),
                 column(6, uiOutput(ns("ui_insp_contrast")))),
        uiOutput(ns("insp_note")),
        plotOutput(ns("insp_plot"), height = "480px"),
        .dc_dl_ui(ns, "insp_plot"),
        DTOutput(ns("insp_tbl"))),
      tabPanel("Pathway inspector",
        fluidRow(column(6, uiOutput(ns("ui_pw_insp_src"))),
                 column(6, uiOutput(ns("ui_pw_insp_contrast")))),
        uiOutput(ns("pw_insp_note")),
        plotOutput(ns("pw_insp_plot"), height = "500px"),
        .dc_dl_ui(ns, "pw_insp_plot"),
        DTOutput(ns("pw_insp_tbl"))),
      tabPanel("Contrast comparison",
        fluidRow(
          column(4, uiOutput(ns("ui_corr_x"))),
          column(4, uiOutput(ns("ui_corr_y"))),
          column(4, radioButtons(ns("corr_src"), "Score source",
                     c("regulators" = "tf", "pathways" = "pw"), inline = TRUE))),
        checkboxInput(ns("corr_sig"),
          "Only regulators significant (p < 0.05) in >=1 selected contrast", TRUE),

        tags$h4("1. Signature concordance  (does one contrast resemble / reverse another?)"),
        fluidRow(
          column(6, uiOutput(ns("ui_sig_ref"))),
          column(6, numericInput(ns("sig_n"),
                    "Signature size (top genes of the reference, by |signed stat|)",
                    200, 50, 1000, step = 50))),
        plotOutput(ns("sig_bar"), height = "360px"),
        .dc_dl_ui(ns, "sig_bar"),
        helpText(HTML(paste(
          "Scores every contrast against the reference contrast's signed gene",
          "signature (ULM). <b>+</b> = resembles the reference, <b>-</b> = reverses it;",
          "the reference scores highest against itself. This is directional and,",
          "unlike a raw correlation, is not inflated when contrasts share a control",
          "group. p-values are approximate (genes are not fully independent) -",
          "rank by the score."))),
        DTOutput(ns("sig_tbl")),

        tags$hr(),
        tags$h4("2. Regulators that differ most between X and Y"),
        DTOutput(ns("corr_tbl")),
        downloadButton(ns("dl_corr"), "Download paired scores (.csv)"),

        tags$hr(),
        tags$h4("3. Exploratory visuals"),
        helpText(HTML(paste(
          "Descriptive only. A positive correlation between two contrasts is",
          "partly structural when they share a group (most \"control vs X\"",
          "designs) - do not read significance into it. Useful signals: a strong",
          "<i>negative</i> correlation, or individual regulators far off the",
          "diagonal (see table 2)."))),
        uiOutput(ns("ui_corr_matrix")),
        plotlyOutput(ns("corr_scatter"), height = "460px"),
        .dc_dl_ui(ns, "corr")),
      tabPanel("Help / interpretation", uiOutput(ns("help_html")))
    )
  )
}

# ---------------------------------------------------------------- server -----
#  dea(): a long data.frame with columns gene, contrast and >=1 of
#  foldchange.log2 (or log2fc) / pvalue / qvalue / effectsize. An optional
#  `algorithm` column (>1 DAA method, as in MS-DAP_visualizer's store$dea) is
#  exposed as a picker and filtered to one method before scoring.
decouplerServer <- function(id, dea) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    RES <- reactiveVal(NULL)

    # normalise column names + collapse to the selected DAA algorithm
    dea_n <- reactive({
      d <- dea(); req(d); d <- as.data.frame(d)
      if (!"foldchange.log2" %in% names(d) && "log2fc" %in% names(d))
        d$foldchange.log2 <- d$log2fc
      if ("algorithm" %in% names(d) && !is.null(input$algo) && nzchar(input$algo))
        d <- d[d$algorithm == input$algo, , drop = FALSE]
      d
    })

    output$ui_algo <- renderUI({
      d <- dea(); req(d)
      if (!"algorithm" %in% names(d)) return(NULL)
      a <- sort(unique(as.character(d$algorithm)))
      if (length(a) < 2) return(NULL)
      selectInput(ns("algo"), "DAA algorithm", choices = a, selected = a[1])
    })

    contrasts_all <- reactive({
      d <- dea_n(); req(d); sort(unique(as.character(d$contrast)))
    })

    output$ui_contrasts <- renderUI({
      cc <- contrasts_all(); req(length(cc))
      pickerInput(ns("contrasts"), "Contrasts", choices = cc, selected = cc,
                  multiple = TRUE, options = list(`actions-box` = TRUE,
                  `selected-text-format` = "count > 3"))
    })
    output$ui_stat <- renderUI({
      d <- dea_n(); req(d)
      st <- .dc_available_stats(d)
      selectInput(ns("stat"), "Per-gene statistic", choices = st, selected = st[1])
    })

    res_choices <- function(panel) {
      all <- .dc_all_resources()
      all <- all[vapply(all, function(x) identical(x$panel, panel), logical(1))]
      stats::setNames(names(all), vapply(all, `[[`, character(1), "label"))
    }
    output$ui_tf_resource <- renderUI({
      ch <- res_choices("tf")
      selectInput(ns("tf_resource"), "Regulator set (TF panels)",
                  choices = ch,
                  selected = if ("collectri" %in% ch) "collectri" else ch[1])
    })
    output$ui_pw_resource <- renderUI({
      ch <- res_choices("pathway")
      selectInput(ns("pw_resource"), "Pathway / signature set",
                  choices = ch,
                  selected = if ("progeny" %in% ch) "progeny" else ch[1])
    })

    observeEvent(input$run, {
      d <- dea_n(); req(d, input$contrasts, input$stat, input$tf_resource, input$pw_resource)
      if (!requireNamespace("decoupleR", quietly = TRUE)) {
        showNotification("Install decoupleR: BiocManager::install('decoupleR')",
                         type = "error", duration = NULL); return(invisible())
      }
      reg   <- .dc_all_resources()
      tf_id <- input$tf_resource; pw_id <- input$pw_resource
      # any failure below is shown to the user instead of failing silently
      ok <- tryCatch(withProgress(message = "Activity inference", value = 0.1, {
        mat <- .dc_matrix(d, input$stat, input$contrasts)
        if (nrow(mat) < 50)
          stop(sprintf("Only %d usable gene x contrast values - check that the DAA table has a %s column and real gene symbols.",
                       nrow(mat), input$stat), call. = FALSE)

        incProgress(0.2, detail = paste("loading", reg[[tf_id]]$label))
        col <- .dc_net(tf_id, DC_ORGANISM)
        if (is.null(col) || !nrow(col))
          stop(sprintf("Could not load '%s' (offline and not cached). Run prefetch_decoupler_resources.R.",
                       reg[[tf_id]]$label), call. = FALSE)
        ov <- sum(unique(rownames(mat)) %in% unique(col$target))
        if (ov < 3 * input$minsize)
          stop(sprintf("Only %d of your %d genes are in '%s'. Symbols look like: %s. Expected human symbols (MS-DAP upper-cases them); if these are non-human native symbols, ortholog-map before scoring.",
                       ov, nrow(mat), reg[[tf_id]]$label,
                       paste(utils::head(rownames(mat), 6), collapse = ", ")), call. = FALSE)
        incProgress(0.15, detail = paste("loading", reg[[pw_id]]$label))
        pw  <- .dc_net(pw_id, DC_ORGANISM)

        incProgress(0.25, detail = "scoring regulator activity")
        warn <- function(m) showNotification(m, type = "warning", duration = 10)
        tf <- .dc_score(mat, col, input$method, input$minsize, warn)
        setDT(tf)

        incProgress(0.2, detail = "scoring pathway / signature activity")
        pwr <- if (!is.null(pw) && nrow(pw) > 0)
          as.data.table(.dc_score(mat, pw,
            if (pw_id %in% c("progeny", "cytosig")) "mlm" else "ulm",
            input$minsize, warn))
          else NULL

        RES(list(tf = tf, pw = pwr, mat = mat, net = as.data.table(col),
                 pw_net = if (!is.null(pw) && nrow(pw) > 0) as.data.table(pw) else NULL,
                 contrasts = input$contrasts, stat = input$stat,
                 tf_id = tf_id, pw_id = pw_id,
                 tf_label = reg[[tf_id]]$label, pw_label = reg[[pw_id]]$label))
        TRUE
      }), error = function(e) {
        showNotification(conditionMessage(e), type = "error", duration = NULL)
        FALSE
      })
      if (isTRUE(ok))
        showNotification("Activity inference complete.", type = "message")
    })

    tf_signed <- reactive(.dc_is_signed(RES()$tf_id))
    pw_signed <- reactive(.dc_is_signed(RES()$pw_id))
    # MSigDB set names are tidied for display (REACTOME_X_Y -> "x y"); other
    # pathway resources (PROGENy/CytoSig/drug) keep their names as-is
    pw_pretty <- reactive(
      if (isTRUE(startsWith(RES()$pw_id %||% "", "msigdb"))) .dc_pretty else identity)

    output$ui_status <- renderUI({
      r <- RES(); if (is.null(r)) return(NULL)
      if (tf_signed() && pw_signed()) return(NULL)
      helpText(strong("Unsigned set in use: "), .dc_score_note(FALSE))
    })

    output$tf_scorenote <- renderUI({
      req(RES()); n <- .dc_score_note(tf_signed()); if (is.null(n)) return(NULL)
      helpText(strong("Note: "), n)
    })
    output$pw_scorenote <- renderUI({
      req(RES()); n <- .dc_score_note(pw_signed()); if (is.null(n)) return(NULL)
      helpText(strong("Note: "), n)
    })
    output$ui_pw_title <- renderUI({
      r <- RES(); if (is.null(r)) return(NULL)
      tags$p(tags$b(r$pw_label))
    })

    # ---- TF heatmap ----
    tf_heat_m <- reactive({ r <- RES(); req(r); .dc_heat_mat(r$tf, input$topn) })
    output$tf_heat <- renderPlotly({
      m <- tf_heat_m()
      plot_ly(x = colnames(m), y = rownames(m), z = m, type = "heatmap",
              colors = colorRamp(c("#2166AC", "white", "#B2182B")),
              zmid = 0, colorbar = list(title = .dc_score_lab(tf_signed()))) %>%
        layout(xaxis = list(tickangle = -40), margin = list(l = 90, b = 120))
    })
    .dc_reg_dl(output, "tf_heat",
      reactive(.dc_heat_gg(tf_heat_m(), .dc_score_lab(tf_signed()))),
      "regulator_activity_heatmap", w = 9, h = 8)

    # ---- TF table ----
    output$ui_tf_contrast <- renderUI(selectInput(ns("tfc"), "Contrast",
      choices = RES()$contrasts %||% NULL))
    output$tf_tbl <- renderDT({
      r <- RES(); req(r, input$tfc)
      t <- r$tf[condition == input$tfc][order(-abs(score))]
      t[, `:=`(score = round(score, 3), p_value = signif(p_value, 3))]
      tt <- t[, .(TF = source, score, p_value)]
      data.table::setnames(tt, "score", paste0(.dc_score_lab(tf_signed()), "_score"))
      datatable(tt, rownames = FALSE, filter = "top", options = list(pageLength = 20))
    })

    # ---- TF bar ----
    output$ui_tf_contrast2 <- renderUI(selectInput(ns("tfc2"), "Contrast",
      choices = RES()$contrasts %||% NULL))
    tf_bar_gg <- reactive({
      r <- RES(); req(r, input$tfc2)
      t <- r$tf[condition == input$tfc2]
      t <- rbind(head(t[order(-score)], 15), head(t[order(score)], 15))
      t[, source := factor(source, levels = source[order(score)])]
      ggplot(t, aes(score, source, fill = score > 0)) +
        geom_col() +
        scale_fill_manual(values = c("#2166AC", "#B2182B"), guide = "none") +
        geom_vline(xintercept = 0, color = "grey40") +
        labs(x = paste(.dc_score_lab(tf_signed()), "score"), y = NULL,
             title = paste(if (tf_signed()) "Top activated / inhibited regulators -"
                           else "Sets with targets most coordinately up / down -",
                           input$tfc2)) +
        theme_minimal(base_size = 13)
    })
    output$tf_bar <- renderPlot(tf_bar_gg())
    .dc_reg_dl(output, "tf_bar", tf_bar_gg, "regulator_top_movers", w = 9, h = 8)

    # ---- pathway ----
    pw_heat_m <- reactive({
      r <- RES(); req(r, !is.null(r$pw))
      m <- .dc_heat_mat(r$pw, input$topn)
      rownames(m) <- pw_pretty()(rownames(m))
      m
    })
    output$pw_heat <- renderPlotly({
      m <- pw_heat_m()
      plot_ly(x = colnames(m), y = rownames(m), z = m, type = "heatmap",
              colors = colorRamp(c("#2166AC", "white", "#B2182B")), zmid = 0,
              colorbar = list(title = .dc_score_lab(pw_signed()))) %>%
        layout(xaxis = list(tickangle = -40), margin = list(l = 90, b = 120))
    })
    .dc_reg_dl(output, "pw_heat",
      reactive(.dc_heat_gg(pw_heat_m(), .dc_score_lab(pw_signed()))),
      "pathway_activity_heatmap", w = 9, h = 7)
    output$pw_tbl <- renderDT({
      r <- RES(); req(r, !is.null(r$pw))
      p <- copy(r$pw)[order(condition, -abs(score))]
      p[, `:=`(score = round(score, 3), p_value = signif(p_value, 3))]
      pt <- p[, .(pathway = pw_pretty()(source), contrast = condition, score, p_value)]
      data.table::setnames(pt, "score", paste0(.dc_score_lab(pw_signed()), "_score"))
      datatable(pt, rownames = FALSE, filter = "top", options = list(pageLength = 15))
    })

    # ---- regulon inspector ----
    output$ui_insp_tf <- renderUI({
      r <- RES(); req(r)
      selectizeInput(ns("insp_tf"), "Transcription factor",
        choices = sort(unique(r$tf$source)),
        selected = r$tf[order(-abs(score))][1]$source)
    })
    output$ui_insp_contrast <- renderUI(selectInput(ns("insp_c"), "Contrast",
      choices = RES()$contrasts %||% NULL))

    # regulon target frame: mode of regulation + this contrast's log2FC / p /
    # signed score, one row per target gene present in the data. `effect` is an
    # agreement call for SIGNED sets (does the change match the edge sign) and
    # only a direction for UNSIGNED sets (ChEA, .gmt - every mor is 1, so there
    # is no activation direction to agree with).
    insp_data <- reactive({
      r <- RES(); req(r, input$insp_tf, input$insp_c)
      tg <- as.data.table(r$net)[source == input$insp_tf, .(gene = target, mor)]
      validate(need(nrow(tg) > 0, "No targets for that regulator."))
      gs <- .dc_gene_stats(dea_n(), input$insp_c)
      sm <- data.table(gene = rownames(r$mat), stat = r$mat[, input$insp_c])
      m  <- merge(merge(tg, gs, by = "gene", all.x = TRUE), sm, by = "gene")
      validate(need(nrow(m) > 0, "No overlap between this regulon and your data."))
      if (isTRUE(tf_signed())) {
        m[, agree := sign(mor) * sign(stat)]
        m[, effect := ifelse(agree > 0, "supports activation",
                      ifelse(agree < 0, "opposes activation", "neutral"))]
      } else {
        m[, effect := ifelse(stat > 0, "target up",
                      ifelse(stat < 0, "target down", "flat"))]
      }
      m[order(-abs(stat))]
    })

    # volcano of the regulon's targets (decoupleR TF vignette style): log2FC vs
    # -log10 p. Signed sets colour by agreement with the edge sign; unsigned
    # sets colour by direction only (no activation call). Falls back to a target
    # bar when the DAA table has no p-value column.
    insp_gg <- reactive({
      dd <- insp_data(); req(nrow(dd) > 0); tf <- input$insp_tf; cc <- input$insp_c
      signed <- isTRUE(tf_signed())
      cols   <- if (signed) c("supports activation" = "#B2182B",
                              "opposes activation" = "#2166AC", "neutral" = "grey70")
                else        c("target up" = "#B2182B",
                              "target down" = "#2166AC", "flat" = "grey70")
      volc <- mean(is.finite(dd$logfc)) > 0.5 && mean(is.finite(dd$pval)) > 0.5
      if (volc) {
        dv  <- dd[is.finite(logfc) & is.finite(pval)]
        lab <- head(dv[order(-abs(stat))], 15)
        ggplot(dv, aes(logfc, -log10(pmax(pval, 1e-300)), color = effect)) +
          geom_vline(xintercept = 0, linetype = 3, color = "grey60") +
          geom_point(size = 2.3, alpha = 0.85) +
          geom_text(data = lab, aes(label = gene), size = 3, vjust = -0.7,
                    check_overlap = TRUE, show.legend = FALSE) +
          scale_color_manual(values = cols, name = NULL) +
          labs(x = "log2 fold change", y = "-log10 p-value",
               title = if (signed)
                 sprintf("%s targets on the %s volcano", tf, cc)
               else
                 sprintf("%s targets on the %s volcano  (unsigned set - direction only, not an activation call)",
                         tf, cc)) +
          theme_minimal(base_size = 12)
      } else {
        db <- data.table::copy(head(dd[order(-abs(stat))], 40))
        if (signed) {
          db[, grp := ifelse(mor > 0, "activating target", "repressing target")]
          ggplot(db, aes(stat, reorder(gene, stat), fill = grp)) +
            geom_col() +
            scale_fill_manual(values = c("activating target" = "#B2182B",
                                         "repressing target" = "#2166AC"), name = NULL) +
            geom_vline(xintercept = 0, color = "grey40") +
            labs(x = sprintf("%s (%s)", RES()$stat, cc), y = NULL,
                 title = sprintf("%s targets  (no p-value column; activating up / repressing down = active)",
                                 tf)) +
            theme_minimal(base_size = 12)
        } else {
          ggplot(db, aes(stat, reorder(gene, stat), fill = stat > 0)) +
            geom_col() +
            scale_fill_manual(values = c("#2166AC", "#B2182B"), guide = "none") +
            geom_vline(xintercept = 0, color = "grey40") +
            labs(x = sprintf("%s (%s)", RES()$stat, cc), y = NULL,
                 title = sprintf("%s targets  (unsigned set - direction of the target set only)", tf)) +
            theme_minimal(base_size = 12)
        }
      }
    })
    output$insp_plot <- renderPlot(insp_gg())
    .dc_reg_dl(output, "insp_plot", insp_gg, "regulon_inspector", w = 9, h = 7)
    output$insp_note <- renderUI({
      req(RES())
      if (isTRUE(tf_signed()))
        helpText("Each point is one target of the regulon. ",
                 strong("Red"), " = the gene's change agrees with its edge sign ",
                 "(activating target up, or repressing target down), i.e. it ",
                 "supports an 'active' call; ", strong("blue"), " opposes it. ",
                 "A call carried by only two or three red points is weak.")
      else
        helpText(strong("Unsigned set (ChEA / .gmt): "), "every edge weight is ",
                 "1, so there is no activation direction to agree with. Colour ",
                 "is just whether each target went up or down in your data - a ",
                 "positive score means the target set moved up ",
                 strong("as a group"), ", not that the TF is activated.")
    })
    output$insp_tbl <- renderDT({
      dd <- insp_data(); req(nrow(dd) > 0); signed <- isTRUE(tf_signed())
      cols <- if (signed)
        dd[, .(target = gene, mode = ifelse(mor > 0, "+", "-"),
               log2fc = round(logfc, 3), p_value = signif(pval, 3),
               stat = round(stat, 3), effect)]
      else
        dd[, .(target = gene, log2fc = round(logfc, 3), p_value = signif(pval, 3),
               stat = round(stat, 3), direction = effect)]
      datatable(cols, rownames = FALSE, options = list(pageLength = 15))
    })

    # ---- pathway inspector --------------------------------------------------
    # Three views, picked from the set's weights:
    #   continuous-weight signed (PROGENy, CytoSig)  -> weight x stat scatter
    #        (decoupleR pathway vignette / "MAPK" view)
    #   +/-1 signed (Drug Perturbations from GEO)     -> log2FC vs -log10 p
    #        volcano coloured by agreement (decoupleR TF-vignette view)
    #   unsigned (MSigDB, ChEA)                       -> ranked target bar
    output$ui_pw_insp_src <- renderUI({
      r <- RES(); req(r, !is.null(r$pw))
      selectizeInput(ns("pw_insp_src"), "Pathway / signature",
        choices = sort(unique(r$pw$source)),
        selected = r$pw[order(-abs(score))][1]$source)
    })
    output$ui_pw_insp_contrast <- renderUI(selectInput(ns("pw_insp_c"), "Contrast",
      choices = RES()$contrasts %||% NULL))

    pw_insp_data <- reactive({
      r <- RES(); req(r, !is.null(r$pw_net), input$pw_insp_src, input$pw_insp_c)
      tg <- as.data.table(r$pw_net)[source == input$pw_insp_src,
                                    .(gene = target, weight = mor)]
      validate(need(nrow(tg) > 0, "No footprint genes for that set."))
      gs <- .dc_gene_stats(dea_n(), input$pw_insp_c)
      sm <- data.table(gene = rownames(r$mat), stat = r$mat[, input$pw_insp_c])
      m  <- merge(merge(tg, gs, by = "gene", all.x = TRUE), sm, by = "gene")
      validate(need(nrow(m) > 0, "No overlap between this set and your data."))
      if (isTRUE(pw_signed()))
        m[, effect := ifelse(weight == 0 | stat == 0, "neutral",
                      ifelse(sign(weight) == sign(stat), "supports activity",
                             "opposes activity"))]
      else
        m[, effect := ifelse(stat > 0, "gene up",
                      ifelse(stat < 0, "gene down", "flat"))]
      m[order(-abs(weight * stat))]
    })

    pw_insp_gg <- reactive({
      dd <- pw_insp_data(); req(nrow(dd) > 0)
      src <- input$pw_insp_src; cc <- input$pw_insp_c
      signed   <- isTRUE(pw_signed())
      weighted <- signed && data.table::uniqueN(round(abs(dd$weight), 8)) > 1
      cols     <- c("supports activity" = "#B2182B",
                    "opposes activity"  = "#2166AC", "neutral" = "grey70")
      if (weighted) {
        lab <- head(dd, 18)
        ggplot(dd, aes(weight, stat, color = effect)) +
          geom_hline(yintercept = 0, linetype = 3, color = "grey60") +
          geom_vline(xintercept = 0, linetype = 3, color = "grey60") +
          geom_point(size = 2.4, alpha = 0.85) +
          geom_text(data = lab, aes(label = gene), size = 3, vjust = -0.7,
                    check_overlap = TRUE, show.legend = FALSE) +
          scale_color_manual(values = cols, name = NULL) +
          labs(x = sprintf("%s footprint weight", src),
               y = sprintf("%s  (%s)", RES()$stat, cc),
               title = sprintf("%s  -  top-right / bottom-left genes drive a positive score",
                               src)) +
          theme_minimal(base_size = 12)
      } else if (signed && mean(is.finite(dd$logfc)) > 0.5 &&
                 mean(is.finite(dd$pval)) > 0.5) {
        dv  <- dd[is.finite(logfc) & is.finite(pval)]
        lab <- head(dv[order(-abs(stat))], 15)
        ggplot(dv, aes(logfc, -log10(pmax(pval, 1e-300)), color = effect)) +
          geom_vline(xintercept = 0, linetype = 3, color = "grey60") +
          geom_point(size = 2.3, alpha = 0.85) +
          geom_text(data = lab, aes(label = gene), size = 3, vjust = -0.7,
                    check_overlap = TRUE, show.legend = FALSE) +
          scale_color_manual(values = cols, name = NULL) +
          labs(x = "log2 fold change", y = "-log10 p-value",
               title = sprintf("%s signature on the %s volcano", src, cc)) +
          theme_minimal(base_size = 12)
      } else if (signed) {
        db <- data.table::copy(head(dd[order(-abs(stat))], 40))
        db[, dir := ifelse(weight > 0, "up in signature", "down in signature")]
        ggplot(db, aes(stat, reorder(gene, stat), fill = dir)) +
          geom_col() +
          scale_fill_manual(values = c("up in signature" = "#B2182B",
                                       "down in signature" = "#2166AC"), name = NULL) +
          geom_vline(xintercept = 0, color = "grey40") +
          labs(x = sprintf("%s  (%s)", RES()$stat, cc), y = NULL,
               title = sprintf("%s signature genes  (up-genes up + down-genes down = resembles it)",
                               src)) +
          theme_minimal(base_size = 12)
      } else {
        db <- head(dd[order(-abs(stat))], 30)
        ggplot(db, aes(stat, reorder(gene, stat), fill = stat > 0)) +
          geom_col() +
          scale_fill_manual(values = c("#2166AC", "#B2182B"), guide = "none") +
          geom_vline(xintercept = 0, color = "grey40") +
          labs(x = sprintf("%s  (%s)", RES()$stat, cc), y = NULL,
               title = sprintf("%s targets (unsigned set - group direction only)", src)) +
          theme_minimal(base_size = 12)
      }
    })
    output$pw_insp_plot <- renderPlot(pw_insp_gg())
    .dc_reg_dl(output, "pw_insp_plot", pw_insp_gg, "pathway_inspector", w = 9, h = 7)
    output$pw_insp_note <- renderUI({
      r <- RES(); req(r); id <- r$pw_id %||% ""
      if (id %in% c("progeny", "cytosig"))
        helpText("Each point is a footprint gene: its set weight (x) vs your ",
                 "per-gene statistic (y), as in the decoupleR pathway vignette. ",
                 "Genes in the top-right and bottom-left quadrants push the ",
                 "activity score up; off-quadrant genes pull it down.")
      else if (isTRUE(pw_signed()))
        helpText("This signature has up / down membership only (weight +/-1), ",
                 "so it gets the volcano view: ", strong("red"), " = the gene ",
                 "moved the way the signature predicts (a signature 'up' gene ",
                 "up, or a 'down' gene down) - your contrast resembles it; ",
                 strong("blue"), " = it reverses the signature.")
      else
        helpText(strong("Unsigned set: "), "every weight is 1, so only the ",
                 "target statistic is shown - the direction the set's genes ",
                 "moved as a group, not an activation call.")
    })
    output$pw_insp_tbl <- renderDT({
      dd <- pw_insp_data(); req(nrow(dd) > 0)
      cols <- if (isTRUE(pw_signed()))
        dd[, .(gene, weight = round(weight, 3),
               log2fc = round(logfc, 3), p_value = signif(pval, 3),
               stat = round(stat, 3), effect)]
      else
        dd[, .(gene, log2fc = round(logfc, 3), p_value = signif(pval, 3),
               stat = round(stat, 3), direction = effect)]
      datatable(cols, rownames = FALSE, options = list(pageLength = 15))
    })

    # ---- contrast correlation ----
    output$ui_corr_x <- renderUI({
      cc <- RES()$contrasts %||% NULL; req(length(cc))
      selectInput(ns("corr_x"), "Contrast X", choices = cc, selected = cc[1])
    })
    output$ui_corr_y <- renderUI({
      cc <- RES()$contrasts %||% NULL; req(length(cc))
      selectInput(ns("corr_y"), "Contrast Y", choices = cc,
                  selected = cc[min(2L, length(cc))])
    })

    corr_src_dt <- reactive({
      r <- RES(); req(r)
      d <- if (identical(input$corr_src, "pw")) r$pw else r$tf
      validate(need(!is.null(d) && nrow(d) > 0,
        "No scores for that source - run inference (pathway table may be empty)."))
      d <- data.table::copy(as.data.table(d))
      if (identical(input$corr_src, "pw")) d[, source := pw_pretty()(source)]
      d
    })

    corr_pts <- reactive({
      d <- corr_src_dt(); x <- input$corr_x; y <- input$corr_y
      req(x, y)
      validate(need(x != y, "Pick two different contrasts."))
      sc <- data.table::dcast(d, source ~ condition, value.var = "score")
      pv <- data.table::dcast(d, source ~ condition, value.var = "p_value")
      req(x %in% names(sc), y %in% names(sc))
      p <- data.table(source = sc$source, sx = sc[[x]], sy = sc[[y]],
                      px = pv[[x]], py = pv[[y]])
      p <- p[is.finite(sx) & is.finite(sy)]
      if (isTRUE(input$corr_sig))
        p <- p[(!is.na(px) & px < 0.05) | (!is.na(py) & py < 0.05)]
      p[, concordance := ifelse(sign(sx) == sign(sy),
          ifelse(sx > 0, "concordant up", "concordant down"), "discordant")]
      p[]
    })

    # score-source label for the correlation axes ("activity" / "enrichment")
    corr_slab <- reactive(.dc_score_lab(.dc_is_signed(
      if (identical(input$corr_src, "pw")) RES()$pw_id else RES()$tf_id)))
    # descriptive only - no CI / p, and every arg forced to length 1 so it can
    # never be sprintf(character(0)) -> a plotly/ggplot "is.character(txt)" crash
    corr_caption <- function(d) {
      rho <- suppressWarnings(stats::cor(d$sx, d$sy, method = "spearman"))
      pr  <- suppressWarnings(stats::cor(d$sx, d$sy))
      sprintf("Spearman rho = %.2f   Pearson r = %.2f   (descriptive; %d regulators, not independent)",
              if (length(rho)) rho else NA_real_,
              if (length(pr)) pr else NA_real_, nrow(d))
    }

    # ---- signature concordance: score every contrast against one contrast's
    # signed gene signature (ULM). Directional, and robust to a shared control.
    output$ui_sig_ref <- renderUI({
      cc <- RES()$contrasts %||% NULL; req(length(cc))
      selectInput(ns("sig_ref"), "Reference contrast", choices = cc, selected = cc[1])
    })
    sig_conc <- reactive({
      r <- RES(); req(r, input$sig_ref, length(r$contrasts) >= 2)
      m <- r$mat; ref <- input$sig_ref
      req(ref %in% colnames(m))
      v <- m[, ref]; v <- v[is.finite(v) & v != 0]
      n <- min(as.integer(input$sig_n %||% 200L), length(v))
      req(n >= 10)
      top <- names(sort(abs(v), decreasing = TRUE))[seq_len(n)]
      net <- data.frame(source = "signature", target = top,
                        mor = sign(v[top]), stringsAsFactors = FALSE)
      res <- tryCatch(decoupleR::run_ulm(m, net, .source = "source",
               .target = "target", .mor = "mor", minsize = 10),
               error = function(e) NULL)
      validate(need(!is.null(res) && nrow(res) > 0, "Could not score the signature."))
      data.table::as.data.table(res)[order(-score)]
    })
    sig_bar_gg <- reactive({
      r <- sig_conc(); req(nrow(r) > 0)
      ref <- input$sig_ref
      r[, is_ref := condition == ref]
      r[, condition := factor(condition, levels = condition[order(score)])]
      ggplot(r, aes(score, condition, fill = is_ref)) +
        geom_col() +
        scale_fill_manual(values = c(`FALSE` = "#4472A8", `TRUE` = "grey55"),
                          guide = "none") +
        geom_vline(xintercept = 0, color = "grey40") +
        labs(x = "concordance with reference signature  (ULM t;  + resembles, - reverses)",
             y = NULL, title = paste("Contrasts vs signature of:", ref,
                                     " (grey = reference itself)")) +
        theme_minimal(base_size = 12)
    })
    output$sig_bar <- renderPlot(sig_bar_gg())
    .dc_reg_dl(output, "sig_bar", sig_bar_gg, "signature_concordance", w = 8, h = 6)
    output$sig_tbl <- renderDT({
      r <- sig_conc(); req(nrow(r) > 0)
      out <- r[, .(contrast = condition, concordance = round(score, 2),
                   p_value = signif(p_value, 3))]
      datatable(out, rownames = FALSE, options = list(pageLength = 15),
        caption = sprintf("signature = top %s genes of %s",
                          min(as.integer(input$sig_n %||% 200L), 99999), input$sig_ref))
    })

    output$ui_corr_matrix <- renderUI({
      if (length(RES()$contrasts %||% NULL) < 3) return(NULL)
      plotlyOutput(ns("corr_matrix"), height = "320px")
    })
    output$corr_matrix <- renderPlotly({
      d <- corr_src_dt(); cc <- RES()$contrasts
      w <- data.table::dcast(d, source ~ condition, value.var = "score")
      m <- as.matrix(w[, -1]); rownames(m) <- w$source
      m <- m[, intersect(cc, colnames(m)), drop = FALSE]
      req(ncol(m) >= 2)
      cm <- suppressWarnings(stats::cor(m, use = "pairwise.complete.obs",
                                        method = "spearman"))
      plot_ly(x = colnames(cm), y = rownames(cm), z = cm, type = "heatmap",
              zmin = -1, zmax = 1, zmid = 0,
              colors = colorRamp(c("#2166AC", "white", "#B2182B")),
              colorbar = list(title = "Spearman")) %>%
        layout(title = "Contrast x contrast correlation (all regulators)",
               xaxis = list(tickangle = -40), margin = list(l = 110, b = 110))
    })

    output$corr_scatter <- renderPlotly({
      d <- corr_pts(); req(nrow(d) > 3)
      req(stats::sd(d$sx, na.rm = TRUE) > 0, stats::sd(d$sy, na.rm = TRUE) > 0)
      x <- as.character(input$corr_x)[1]; y <- as.character(input$corr_y)[1]
      slab <- corr_slab()
      fit <- stats::lm(sy ~ sx, d)
      rng <- range(c(d$sx, d$sy), na.rm = TRUE)
      lab <- d[order(-(sx^2 + sy^2))][seq_len(min(10L, .N))]
      plot_ly(d, x = ~sx, y = ~sy, type = "scatter", mode = "markers",
              color = ~concordance, text = ~source, hoverinfo = "text+x+y",
              marker = list(size = 7, opacity = 0.8)) %>%
        add_lines(x = rng, y = rng, inherit = FALSE, name = "y = x",
                  line = list(dash = "dot", color = "grey"),
                  hoverinfo = "none") %>%
        add_lines(x = rng,
                  y = as.numeric(predict(fit, data.frame(sx = rng))),
                  inherit = FALSE, name = "fit",
                  line = list(color = "black"), hoverinfo = "none") %>%
        add_annotations(x = lab$sx, y = lab$sy, text = as.character(lab$source),
                        showarrow = FALSE, xshift = 4, yshift = 9,
                        font = list(size = 10)) %>%
        layout(title = corr_caption(d),
          xaxis = list(title = paste0(slab, " score - ", x)),
          yaxis = list(title = paste0(slab, " score - ", y)))
    })

    corr_gg <- reactive({
      d <- corr_pts(); req(nrow(d) > 3)
      x <- as.character(input$corr_x)[1]; y <- as.character(input$corr_y)[1]
      slab <- corr_slab()
      lab <- d[order(-(sx^2 + sy^2))][seq_len(min(12L, .N))]
      ggplot(d, aes(sx, sy)) +
        geom_abline(slope = 1, intercept = 0, linetype = 3, color = "grey60") +
        geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                    color = "black", linewidth = 0.5) +
        geom_point(aes(color = concordance), size = 1.9, alpha = 0.85) +
        geom_text(data = lab, aes(label = source), size = 3, vjust = -0.7) +
        labs(x = paste0(slab, " score - ", x), y = paste0(slab, " score - ", y),
             color = NULL, title = corr_caption(d)) +
        theme_minimal(base_size = 12)
    })
    .dc_reg_dl(output, "corr", corr_gg, "contrast_correlation", w = 8, h = 7)

    output$corr_tbl <- renderDT({
      d <- corr_pts(); req(nrow(d) > 0)
      d <- d[order(-abs(sx - sy))]
      out <- d[, .(regulator = source,
                   score_X = round(sx, 3), p_X = signif(px, 3),
                   score_Y = round(sy, 3), p_Y = signif(py, 3),
                   delta = round(sx - sy, 3), concordance)]
      datatable(out, rownames = FALSE, filter = "top",
                options = list(pageLength = 20),
                caption = sprintf("X = %s   Y = %s", input$corr_x, input$corr_y))
    })

    output$dl_corr <- downloadHandler(
      filename = function() sprintf("contrast_correlation_%s.csv", Sys.Date()),
      content  = function(f) utils::write.csv(corr_pts(), f, row.names = FALSE))

    # ---- help ----
    output$help_html <- renderUI(HTML(.DC_HELP_HTML))

    output$dl_tf <- downloadHandler(
      filename = function() sprintf("TF_activity_%s.csv", Sys.Date()),
      content = function(f) utils::write.csv(RES()$tf, f, row.names = FALSE))
    output$dl_pw <- downloadHandler(
      filename = function() sprintf("pathway_activity_%s.csv", Sys.Date()),
      content = function(f) utils::write.csv(RES()$pw, f, row.names = FALSE))
  })
}
