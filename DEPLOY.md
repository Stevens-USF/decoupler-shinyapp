# decoupler_shinyapp — deployment

Standalone Shiny app for upstream regulator + pathway activity inference on
MS-DAP differential-abundance results (a free "IPA Upstream Regulator /
Canonical Pathways" alternative). All prior-knowledge resources are bundled, so
the app needs no internet at run time.

## Contents

```
app.R                          entry point (shinyapps.io looks for this name)
decoupler_module.R             all UI + server logic, sourced by app.R
decoupler_cache/               signed/curated resources as <id>_human.csv (source,target,mor)
  collectri_human.csv          CollecTRI  TF regulon (signed)
  progeny_human.csv            PROGENy    14 signalling pathways (signed)
  cytosig_human.csv            CytoSig    cytokine-response signatures (signed)
  drugpert_geo_human.csv       Drug Perturbations from GEO (signed, up/down)
  msigdb_hallmark_human.csv    MSigDB Hallmark
  msigdb_kegg_human.csv        MSigDB KEGG
  msigdb_reactome_human.csv    MSigDB Reactome
  msigdb_gobp_human.csv        MSigDB GO:BP   (55 MB — see "Memory" below)
geneset_cache/
  ChEA_2022.gmt                ChIP-based TF targets (unsigned)
prefetch_decoupler_resources.R  rebuild decoupler_cache/ from OmniPath/decoupleR (run locally)
prefetch_drug_signatures.R      rebuild drugpert_geo_human.csv from Enrichr (run locally)
```

Bundle size ≈ 78 MB (well under the 1 GB limit).

## Run locally

```r
setwd("path/to/decoupler_shinyapp")
shiny::runApp()
```

## Deploy to shinyapps.io

```r
install.packages("rsconnect")
rsconnect::setAccountInfo(name="<acct>", token="<token>", secret="<secret>")   # from shinyapps.io > Account > Tokens

# decoupleR is a Bioconductor package — make sure rsconnect can see it:
options(repos = BiocManager::repositories())          # needs install.packages("BiocManager")
# (decoupleR, and the CRAN packages below, must be installed in this R library
#  so rsconnect can snapshot their versions)

rsconnect::deployApp(
  appDir     = "path/to/decoupler_shinyapp",
  appName    = "decoupler",
  appFiles   = c("app.R", "decoupler_module.R",
                 list.files("decoupler_cache", full.names = TRUE),
                 list.files("geneset_cache",  full.names = TRUE))
)
```

`appFiles` excludes the two `prefetch_*.R` scripts (maintenance only). Drop the
`appFiles` argument to deploy everything in the folder.

### Package dependencies

CRAN: `shiny shinydashboard shinyWidgets DT plotly ggplot2 data.table readxl`
Bioconductor: `decoupleR`

`OmnipathR` is **not** required at run time (all resources are cached). It is a
dependency of `decoupleR` so it will be installed regardless; that is fine.

## Memory (shinyapps.io free tier = 1 GB)

`msigdb_gobp_human.csv` is 55 MB and expands to thousands of gene sets; scoring
it (ULM over ~10k genes) is the heaviest operation in the app and can exhaust a
1 GB instance. If you hit out-of-memory restarts:

* delete `decoupler_cache/msigdb_gobp_human.csv`, **and**
* comment out the `msigdb_gobp` entry in `decoupler_module.R` `.DC_RESOURCES`
  (otherwise selecting it shows a load error)

or move to a paid instance with more RAM.

## Refreshing / adding resources

* `.gmt` files in `geneset_cache/` are auto-discovered — add one, redeploy, it
  appears in the regulator dropdown (unsigned).
* Signed CSVs and the MSigDB collections are fixed entries in
  `decoupler_module.R` `.DC_RESOURCES`. Rebuild their CSVs with
  `prefetch_decoupler_resources.R` / `prefetch_drug_signatures.R` (run locally,
  needs internet), then redeploy.
* **DoRothEA is disabled** in this build (no bundled CSV, OmniPath unreachable at
  run time). To enable: run `prefetch_decoupler_resources.R` to create
  `decoupler_cache/dorothea_human.csv`, then un-comment the `dorothea` entry in
  `decoupler_module.R`.

## Notes

* Regulons are scored in **human** symbol space. MS-DAP upper-cases gene
  symbols, so mouse datasets map to human orthologs directly for ~94% of genes.
* Always validate on a positive control (a TNFα / LPS contrast should call
  NF-κB / STAT1 / interferon active) before trusting novel calls.
