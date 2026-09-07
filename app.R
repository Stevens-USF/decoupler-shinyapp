# =============================================================================
#  app.R  --  standalone "upstream regulator activity" app  (shinyapps.io entry)
#  Wraps decoupler_module.R. Accepts an MS-DAP differential_abundance_analysis.xlsx
#  (or the *_dea long-format .csv/.tsv it also writes).
#
#  Local:   shiny::runApp()          # from this folder
#  Deploy:  rsconnect::deployApp()   # see DEPLOY.md
#  Needs: decoupleR (Bioconductor) + shiny, shinydashboard, shinyWidgets, DT,
#         plotly, ggplot2, data.table, readxl.  All prior-knowledge resources
#         are bundled under decoupler_cache/ and geneset_cache/, so OmnipathR
#         is NOT needed at run time.
# =============================================================================

options(shiny.maxRequestSize = 512 * 1024^2)
suppressPackageStartupMessages({
  library(shiny); library(shinydashboard); library(shinyWidgets)
  library(DT); library(plotly); library(ggplot2); library(data.table); library(readxl)
})
source("decoupler_module.R")

# ---- MS-DAP statistics sheet (wide)  ->  long (gene, contrast, stats) -------
parse_msdap_dea <- function(path, algo = NULL) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("csv", "tsv", "txt")) {
    d <- as.data.frame(data.table::fread(path))
    need <- c("gene", "contrast")
    if (all(need %in% names(d))) return(list(long = d, algos = NA))
  }
  dea <- as.data.frame(read_excel(path, sheet = "statistics"))
  cn  <- names(dea)
  gcol <- intersect(c("gene_symbols_or_id", "gene_symbols", "gene_symbol"), cn)[1]
  m <- regmatches(cn, regexec(
    "^(foldchange\\.log2|pvalue|qvalue|effectsize)_(\\w+)_contrast: (.+?) # ", cn))
  keep <- lengths(m) == 4
  info <- data.table(col = cn[keep],
                     stat = sapply(m[keep], `[`, 2),
                     algo = sapply(m[keep], `[`, 3),
                     contrast = sapply(m[keep], `[`, 4))
  algos <- sort(unique(info$algo))
  use   <- if (!is.null(algo) && algo %in% algos) algo else algos[1]
  info  <- info[algo == use]

  long <- rbindlist(lapply(unique(info$contrast), function(cc) {
    sub <- info[contrast == cc]
    dt <- data.table(gene = as.character(dea[[gcol]]), contrast = cc)
    for (st in unique(sub$stat)) dt[[st]] <- dea[[sub[stat == st, col][1]]]
    dt
  }), fill = TRUE)
  setnames(long, "foldchange.log2", "foldchange.log2", skip_absent = TRUE)
  list(long = long, algos = algos, used = use)
}

ui <- dashboardPage(
  dashboardHeader(title = "decoupleR (MS-DAP)"),
  dashboardSidebar(sidebarMenu(
    menuItem("Data", tabName = "data", icon = icon("upload")),
    menuItem("decoupleR", tabName = "dc", icon = icon("diagram-project")))),
  dashboardBody(tabItems(
    tabItem("data",
      box(width = 6, title = "Load MS-DAP DAA", status = "primary", solidHeader = TRUE,
        fileInput("f", "differential_abundance_analysis.xlsx",
                  accept = c(".xlsx", ".csv", ".tsv", ".txt")),
        uiOutput("ui_algo"),
        actionButton("load", "Load", class = "btn-primary", icon = icon("check"))),
      box(width = 6, title = "Status", status = "info", solidHeader = TRUE,
        verbatimTextOutput("msg"))),
    tabItem("dc", decouplerTabUI("dc"))
  ))
)

server <- function(input, output, session) {
  raw <- reactiveVal(NULL)
  observeEvent(input$f, raw(parse_msdap_dea(input$f$datapath)))
  output$ui_algo <- renderUI({
    r <- raw(); if (is.null(r) || length(r$algos) < 2 || any(is.na(r$algos))) return(NULL)
    selectInput("algo", "DAA algorithm", choices = r$algos, selected = r$used)
  })
  store <- reactiveValues(long = NULL)
  observeEvent(input$load, {
    req(input$f)
    p <- parse_msdap_dea(input$f$datapath, algo = input$algo %||% NULL)
    store$long <- as.data.frame(p$long)
  })
  output$msg <- renderText({
    if (is.null(store$long)) return("Load a file.")
    d <- store$long
    sprintf("%d rows | %d contrasts | %d unique genes\ncolumns: %s\ncontrasts: %s",
            nrow(d), length(unique(d$contrast)), length(unique(d$gene)),
            paste(names(d), collapse = ", "),
            paste(unique(d$contrast), collapse = " | "))
  })
  decouplerServer("dc", dea = reactive(store$long))
}

shinyApp(ui, server)
