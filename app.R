# dssd Survey Designer — a Shiny front end for the dssd package
# (Distance Sampling Survey Design; Marshall & Rexstad, CRAN / GitHub
#  DistanceDevelopment/dssd).
#
# Workflow: region -> design -> generate transects -> coverage simulation
#           -> export (shapefile / GPX / CSV) + reproducible R script.
#
# Single-file app so it can be exported with shinylive::export().
# Packages: shiny, bslib, DT, sf, dssd, ggplot2, zip

library(shiny)
library(bslib)
library(DT)
library(sf)
library(dssd)
library(ggplot2)

# zip is only called as zip::zipr(); this never-run line makes sure
# Shinylive bundles it into the WebAssembly build.
if (FALSE) library(zip)
# sf and dssd use tibble only as an optional ("Suggests") package, which
# Shinylive does not bundle automatically; this line makes sure it is included.
if (FALSE) library(tibble)

options(shiny.maxRequestSize = 100 * 1024^2)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

line_designs  <- c("Systematic parallel"                 = "systematic",
                   "Random parallel"                     = "random",
                   "Equal-spaced zigzag"                 = "eszigzag",
                   "Equal-spaced zigzag (complementary)" = "eszigzagcom",
                   "Segmented grid"                      = "segmentedgrid")
point_designs <- c("Systematic grid" = "systematic",
                   "Random"          = "random")

example_regions <- c(
  "St Andrews bay (1 stratum, no CRS)" = "StAndrew",
  "Area R, 6 strata (projected)"       = "AreaRProjStrata",
  "Tentsmuir, 2 strata (lat/lon)"      = "TentsmuirUnproj"
)

# Parse "1" or "1, 2, 3" into numeric; NULL if blank. Errors are user-facing.
parse_num <- function(txt, label, n_strata, allow_blank = TRUE) {
  txt <- trimws(txt %||% "")
  if (!nzchar(txt)) {
    if (allow_blank) return(numeric(0))
    stop(label, " is required.", call. = FALSE)
  }
  vals <- suppressWarnings(as.numeric(strsplit(txt, "[,;[:space:]]+")[[1]]))
  if (anyNA(vals)) stop(label, ": could not read '", txt, "' as numbers.", call. = FALSE)
  if (!length(vals) %in% c(1L, n_strata))
    stop(label, ": supply 1 value or ", n_strata, " (one per stratum); got ",
         length(vals), ".", call. = FALSE)
  vals
}

# Read an uploaded set of files (zip, gpkg, geojson, or shapefile parts)
read_upload <- function(files, layer = NULL) {
  dir <- tempfile("upload_"); dir.create(dir)
  # restore original names so shapefile parts sit together
  file.copy(files$datapath, file.path(dir, files$name))
  zips <- list.files(dir, "\\.zip$", full.names = TRUE, ignore.case = TRUE)
  for (z in zips) utils::unzip(z, exdir = dir)
  cand <- list.files(dir, "\\.(shp|gpkg|geojson|json|kml)$",
                     full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  cand <- cand[!grepl("__MACOSX", cand)]
  if (!length(cand))
    stop("No .shp, .gpkg, .geojson or .kml found. For shapefiles, upload a .zip ",
         "or select all parts (.shp, .shx, .dbf, .prj) together.", call. = FALSE)
  if (length(cand) > 1)
    stop("More than one spatial file found (", paste(basename(cand), collapse = ", "),
         "). Upload one region at a time.", call. = FALSE)
  layers <- sf::st_layers(cand)$name
  lyr <- if (!is.null(layer) && layer %in% layers) layer else layers[1]
  list(sf = sf::st_zm(sf::st_read(cand, layer = lyr, quiet = TRUE)), path = cand,
       layers = layers, layer = lyr)
}

# UTM EPSG code from a geographic layer's centroid
utm_epsg <- function(x) {
  ctr <- suppressWarnings(sf::st_coordinates(
    sf::st_centroid(sf::st_union(sf::st_transform(x, 4326)))))
  zone <- floor((ctr[1] + 180) / 6) + 1
  if (ctr[2] < 0) 32700 + zone else 32600 + zone
}

fmt <- function(x, d = 0) formatC(x, format = "f", digits = d, big.mark = ",")

# Run an expression, turning errors/warnings into notifications.
# Returns NULL on error; warnings are shown but the value is kept.
safely <- function(expr, what) {
  withCallingHandlers(
    tryCatch(expr, error = function(e) {
      showNotification(paste0(what, " failed: ", conditionMessage(e)),
                       type = "error", duration = 12)
      NULL
    }),
    warning = function(w) {
      showNotification(paste0(what, ": ", conditionMessage(w)),
                       type = "warning", duration = 10)
      invokeRestart("muffleWarning")
    })
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# Value of a reactive, or NULL if it errors / hasn't fired (req() etc.)
maybe <- function(r) tryCatch(r(), error = function(e) NULL)

# Dissolve features into one MULTIPOLYGON per value of grp, keeping first-seen order
dissolve_by <- function(x, grp) {
  grp <- as.character(grp)
  grp[is.na(grp)] <- "NA"
  g <- split(sf::st_geometry(x), factor(grp, levels = unique(grp)))
  geoms <- do.call(c, lapply(g, sf::st_union))
  sf::st_sf(stratum = names(g), geometry = geoms)
}

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------

ui <- page_navbar(
  title = "dssd Survey Designer",
  theme = bs_theme(version = 5, preset = "flatly"),
  fillable = FALSE,

  # ---- 1. Region ----------------------------------------------------------
  nav_panel(
    "1. Region",
    layout_sidebar(
      sidebar = sidebar(
        width = 340,
        radioButtons("region_src", "Region source",
                     c("Upload file" = "upload", "dssd example" = "example")),
        conditionalPanel(
          "input.region_src == 'upload'",
          fileInput("region_file", "Shapefile (.zip or all parts), GeoPackage, GeoJSON",
                    multiple = TRUE,
                    accept = c(".zip", ".shp", ".shx", ".dbf", ".prj", ".cpg",
                               ".gpkg", ".geojson", ".json", ".kml")),
          uiOutput("layer_ui")
        ),
        conditionalPanel(
          "input.region_src == 'example'",
          selectInput("example", NULL, example_regions)
        ),
        textInput("region_name", "Region name", "region"),
        uiOutput("strata_ui"),
        hr(),
        uiOutput("crs_ui"),
        actionButton("build_region", "Build region", class = "btn-primary w-100")
      ),
      card(card_header("Survey region"), plotOutput("region_plot", height = 520)),
      card(card_header("Strata"), DTOutput("strata_tbl"),
           verbatimTextOutput("region_info"))
    )
  ),

  # ---- 2. Design ----------------------------------------------------------
  nav_panel(
    "2. Design",
    layout_sidebar(
      sidebar = sidebar(
        width = 360,
        radioButtons("transect_type", "Transect type",
                     c("Line" = "line", "Point" = "point"), inline = TRUE),
        selectInput("design", "Design", line_designs),
        radioButtons("effort_type", "Specify effort by",
                     c("Number of samplers" = "samplers",
                       "Total line length"  = "line.length",
                       "Spacing"            = "spacing")),
        textInput("effort_val", "Effort value(s)", "20"),
        uiOutput("strata_order_hint"),
        conditionalPanel(
          "input.design == 'segmentedgrid'",
          textInput("seg_length", "Segment length(s)", ""),
          numericInput("seg_threshold", "Segment threshold (%)", 50, 0, 100)
        ),
        conditionalPanel(
          "input.design == 'eszigzag' || input.design == 'eszigzagcom'",
          radioButtons("bounding_shape", "Zigzag bounding shape",
                       c("Rectangle" = "rectangle", "Convex hull" = "convex.hull"),
                       inline = TRUE)
        ),
        textInput("design_angle", "Design angle(s) (degrees)", "0"),
        radioButtons("edge_protocol", "Edge protocol",
                     c("Minus" = "minus", "Plus" = "plus"), inline = TRUE),
        numericInput("truncation", "Truncation distance (region units)", 50, min = 0),
        textInput("effort_alloc", "Effort allocation (proportions, optional)", ""),
        helpText("Used only when a single total effort is given for a multi-strata",
                 "region. Blank = allocate by stratum area."),
        hr(),
        numericInput("seed", "Random seed (blank = none)", 123),
        actionButton("generate", "Generate transects", class = "btn-primary w-100")
      ),
      navset_card_tab(
        nav_panel("Map", plotOutput("transect_plot", height = 560)),
        nav_panel("Summary", DTOutput("transect_summary"),
                  verbatimTextOutput("transect_print")),
        nav_panel("Coordinates", DTOutput("transect_coords"))
      )
    )
  ),

  # ---- 3. Coverage --------------------------------------------------------
  nav_panel(
    "3. Coverage",
    layout_sidebar(
      sidebar = sidebar(
        width = 320,
        helpText("Simulates the design many times to estimate the coverage",
                 "probability at a grid of points and the distribution of effort."),
        numericInput("n_grid", "Approx. number of grid points", 1000, 50, 20000, 50),
        numericInput("reps", "Replicates", 100, 2, 5000, 10),
        helpText("Few replicates give noisy coverage scores; use several hundred",
                 "or more for a final assessment. Run time grows with replicates",
                 "and grid size."),
        actionButton("run_cov", "Run coverage simulation", class = "btn-primary w-100")
      ),
      layout_columns(
        col_widths = c(7, 5),
        card(card_header("Coverage probability"), plotOutput("cov_map", height = 480)),
        card(card_header("Coverage score distribution"),
             plotOutput("cov_hist", height = 300), uiOutput("cov_flag"))
      ),
      card(
        card_header(
          div(class = "d-flex justify-content-between align-items-center",
              "Design statistics",
              div(style = "width:260px", selectInput("stat", NULL, choices = NULL)))
        ),
        DTOutput("cov_stats")
      )
    )
  ),

  # ---- 4. Export ----------------------------------------------------------
  nav_panel(
    "4. Export",
    layout_columns(
      col_widths = c(5, 7),
      card(
        card_header("Download transects"),
        downloadButton("dl_shp", "Shapefile (.zip)", class = "mb-2"),
        downloadButton("dl_gpx", "GPX (WGS84)", class = "mb-2"),
        downloadButton("dl_csv", "CSV of coordinates", class = "mb-2"),
        helpText("GPX: lines are written as tracks, points as waypoints, reprojected",
                 "to WGS84. Requires the region to have a CRS. CSV includes",
                 "projected coordinates and, where a CRS exists, longitude/latitude."),
        hr(),
        downloadButton("dl_cov", "Coverage statistics (.csv)"),
        downloadButton("dl_script", "R script (.R)")
      ),
      card(card_header("Reproducible R code"), verbatimTextOutput("r_code"))
    )
  ),

  nav_spacer(),
  nav_item(tags$a("dssd on GitHub", href = "https://github.com/DistanceDevelopment/dssd",
                  target = "_blank"))
)

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------

server <- function(input, output, session) {

  # ---- Raw input layer ----------------------------------------------------
  raw <- reactive({
    if (input$region_src == "example") {
      p <- system.file("extdata", paste0(input$example, ".shp"), package = "dssd")
      list(sf = sf::st_read(p, quiet = TRUE), path = p, layers = NULL, layer = NULL)
    } else {
      req(input$region_file)
      safely(read_upload(input$region_file, input$layer), "Reading file")
    }
  })

  output$layer_ui <- renderUI({
    r <- raw(); req(r, length(r$layers) > 1)
    selectInput("layer", "Layer", r$layers, r$layer)
  })

  output$strata_ui <- renderUI({
    r <- raw(); req(r)
    x <- sf::st_drop_geometry(r$sf)
    cols <- names(x)
    n <- nrow(x)
    # default: a field that is unique per feature, else one stratum per feature
    uniq <- cols[vapply(x, function(v) length(unique(v)) == n, logical(1))]
    sel <- if (n == 1) "(none)" else if (length(uniq)) uniq[1] else "(each)"
    tagList(
      selectInput("strata_col", "Strata field",
                  c("Single stratum (dissolve all)" = "(none)",
                    "Each feature is a stratum" = "(each)", cols), sel),
      helpText(n, "feature(s). Choosing a field dissolves features by its",
               "values, one stratum per value.")
    )
  })

  output$crs_ui <- renderUI({
    r <- raw(); req(r)
    crs <- sf::st_crs(r$sf)
    if (is.na(crs)) {
      tagList(
        tags$b("No CRS found."),
        radioButtons("na_crs", NULL,
                     c("Assign EPSG code" = "assign",
                       "Leave unset (no GPX / lon-lat export)" = "none")),
        conditionalPanel("input.na_crs == 'assign'",
                         numericInput("assign_epsg", "EPSG of the coordinates", NA)),
        conditionalPanel("input.na_crs == 'none'",
                         textInput("na_units", "Coordinate units", "m"))
      )
    } else if (isTRUE(crs$IsGeographic)) {
      tagList(
        tags$b("Geographic CRS (", crs$Name, ")."),
        helpText("dssd needs projected coordinates; the region will be reprojected."),
        radioButtons("proj_mode", "Project to",
                     c("Auto UTM zone" = "utm", "EPSG code" = "epsg")),
        conditionalPanel("input.proj_mode == 'epsg'",
                         numericInput("proj_epsg", "Target EPSG", 7855))
      )
    } else {
      tagList(
        tags$b("Projected CRS: ", crs$Name, " (units: ", crs$units_gdal, ")"),
        checkboxInput("reproject", "Reproject anyway", FALSE),
        conditionalPanel("input.reproject",
                         numericInput("proj_epsg2", "Target EPSG", 7855))
      )
    }
  })

  # ---- Build region -------------------------------------------------------
  region_res <- eventReactive(input$build_region, {
    r <- raw(); req(r)
    safely({
      x <- r$sf
      crs <- sf::st_crs(x)
      units <- character(0)
      proj_code <- NULL                       # for the R script

      if (is.na(crs)) {
        if (identical(input$na_crs, "assign")) {
          if (is.na(input$assign_epsg)) stop("Enter an EPSG code.", call. = FALSE)
          x <- sf::st_set_crs(x, input$assign_epsg)
          proj_code <- list(assign = input$assign_epsg)
          if (isTRUE(sf::st_crs(x)$IsGeographic)) {
            epsg <- utm_epsg(x); x <- sf::st_transform(x, epsg)
            proj_code$to <- epsg
          }
        } else {
          units <- input$na_units %||% "m"
        }
      } else if (isTRUE(crs$IsGeographic)) {
        epsg <- if (identical(input$proj_mode, "epsg")) input$proj_epsg else utm_epsg(x)
        x <- sf::st_transform(x, epsg); proj_code <- list(to = epsg)
      } else if (isTRUE(input$reproject)) {
        x <- sf::st_transform(x, input$proj_epsg2); proj_code <- list(to = input$proj_epsg2)
      }

      if (isTRUE(sf::st_is_longlat(x)))
        stop("Region is still in geographic coordinates.", call. = FALSE)

      x <- sf::st_make_valid(x)
      x <- x[sf::st_dimension(x) == 2, ]
      if (!nrow(x)) stop("No polygon features found.", call. = FALSE)

      scol <- input$strata_col %||% "(none)"
      if (scol == "(none)") {
        x <- sf::st_sf(stratum = input$region_name,
                       geometry = sf::st_union(x))
        strata <- input$region_name
      } else if (scol == "(each)") {
        x <- sf::st_sf(stratum = paste0("S", seq_len(nrow(x))),
                       geometry = sf::st_geometry(x))
        strata <- x$stratum
      } else {
        x <- dissolve_by(x, x[[scol]])
        strata <- x$stratum
      }
      x <- sf::st_cast(x, "MULTIPOLYGON")

      reg <- dssd::make.region(region.name = input$region_name,
                               strata.name = strata, shape = x, units = units)
      list(region = reg, shape = x, strata = strata, proj_code = proj_code,
           strata_col = scol, units = units, src = r)
    }, "Building region")
  })

  region <- reactive({ req(region_res()); region_res()$region })

  output$region_plot <- renderPlot({
    if (!is.null(maybe(region_res))) return(plot(region()))
    r <- raw(); req(r)
    ggplot(r$sf) + geom_sf(fill = "grey85") + theme_minimal() +
      labs(subtitle = "Preview of input — click 'Build region'")
  })

  output$strata_tbl <- renderDT({
    reg <- region()
    datatable(data.frame(Stratum = reg@strata.name,
                         Area = round(reg@area, 1),
                         Units = paste0(reg@units, "^2")),
              rownames = FALSE, options = list(dom = "t"))
  })

  output$region_info <- renderText({
    res <- region_res(); req(res)
    crs <- sf::st_crs(res$shape)
    paste0("Region: ", res$region@region.name,
           "\nCRS: ", if (is.na(crs)) "none" else crs$Name,
           "\nUnits: ", res$region@units,
           "\nTotal area: ", fmt(sum(res$region@area)), " ", res$region@units, "^2")
  })

  # ---- Design controls ----------------------------------------------------
  observeEvent(input$transect_type, {
    if (input$transect_type == "point") {
      updateSelectInput(session, "design", choices = point_designs)
      updateRadioButtons(session, "effort_type",
                         choices = c("Number of samplers" = "samplers",
                                     "Spacing" = "spacing"))
    } else {
      updateSelectInput(session, "design", choices = line_designs)
      updateRadioButtons(session, "effort_type",
                         choices = c("Number of samplers" = "samplers",
                                     "Total line length" = "line.length",
                                     "Spacing" = "spacing"))
    }
  })

  output$strata_order_hint <- renderUI({
    req(region_res())
    s <- region()@strata.name
    if (length(s) > 1)
      helpText("Strata order for per-stratum values:", paste(s, collapse = ", "))
  })

  # Build the list of make.design() arguments from the inputs
  design_args <- reactive({
    reg <- region()
    n <- length(reg@strata.name)
    args <- list(transect.type = input$transect_type,
                 design = input$design,
                 edge.protocol = input$edge_protocol,
                 truncation = input$truncation)
    eff <- parse_num(input$effort_val, "Effort", n, allow_blank = FALSE)
    args[[input$effort_type]] <- eff
    args$design.angle <- parse_num(input$design_angle, "Design angle", n)
    if (!length(args$design.angle)) args$design.angle <- 0
    alloc <- parse_num(input$effort_alloc, "Effort allocation", n)
    if (length(alloc)) {
      if (length(alloc) != n) stop("Effort allocation needs one value per stratum.", call. = FALSE)
      if (abs(sum(alloc) - 1) > 1e-6) stop("Effort allocation must sum to 1.", call. = FALSE)
      args$effort.allocation <- alloc
    }
    if (input$transect_type == "line" && input$design %in% c("eszigzag", "eszigzagcom"))
      args$bounding.shape <- input$bounding_shape
    if (input$transect_type == "line" && input$design == "segmentedgrid") {
      args$seg.length <- parse_num(input$seg_length, "Segment length", n, allow_blank = FALSE)
      args$seg.threshold <- input$seg_threshold
    }
    if (!is.numeric(args$truncation) || is.na(args$truncation) || args$truncation <= 0)
      stop("Truncation must be a positive number.", call. = FALSE)
    args
  })

  make_design <- function(coverage.grid = NULL) {
    args <- design_args()
    args$region <- region()
    if (!is.null(coverage.grid)) args$coverage.grid <- coverage.grid
    do.call(dssd::make.design, args)
  }

  set_seed <- function() if (isTruthy(input$seed)) set.seed(input$seed)

  # ---- Generate transects -------------------------------------------------
  transects <- eventReactive(input$generate, {
    req(region_res())
    withProgress(message = "Generating transects", value = 0.5, {
      safely({
        d <- make_design()
        set_seed()
        list(design = d, transects = dssd::generate.transects(d, quiet = TRUE),
             args = design_args())
      }, "Generating transects")
    })
  })

  tr <- reactive({ req(transects()); transects()$transects })

  output$transect_plot <- renderPlot(plot(region(), tr()))

  output$transect_summary <- renderDT({
    t <- tr()
    df <- data.frame(Stratum = t@strata.names,
                     Samplers = t@samp.count,
                     `Covered area` = round(t@cov.area, 1),
                     `% covered` = round(100 * t@cov.area / t@strata.area, 2),
                     check.names = FALSE)
    if (inherits(t, "Line.Transect") || inherits(t, "Segment.Transect")) {
      df$`Line length` <- round(t@line.length, 1)
      if (.hasSlot(t, "trackline"))       df$`Trackline` <- round(t@trackline, 1)
      if (.hasSlot(t, "cyclictrackline")) df$`Cyclic trackline` <- round(t@cyclictrackline, 1)
    }
    if (nrow(df) > 1) {
      tot <- lapply(df, function(v) if (is.numeric(v)) sum(v) else "Total")
      tot$`% covered` <- round(100 * sum(t@cov.area) / sum(t@strata.area), 2)
      df <- rbind(df, as.data.frame(tot, check.names = FALSE))
    }
    datatable(df, rownames = FALSE, options = list(dom = "t", pageLength = 50))
  })

  output$transect_print <- renderPrint(print(tr()))

  # coordinates via dssd's own CSV writer, plus lon/lat where possible
  coords <- reactive({
    t <- tr()
    f <- tempfile(fileext = ".csv")
    dssd::write.transects(t, dsn = f)
    df <- utils::read.csv(f)
    crs <- sf::st_crs(t@samplers)
    if (!is.na(crs)) {
      ll <- function(x, y) {
        p <- sf::st_transform(sf::st_as_sf(data.frame(x = x, y = y),
                                           coords = c("x", "y"), crs = crs), 4326)
        sf::st_coordinates(p)
      }
      if (all(c("x", "y") %in% names(df))) {
        c1 <- ll(df$x, df$y); df$lon <- c1[, 1]; df$lat <- c1[, 2]
      } else if (all(c("x.start", "y.start", "x.end", "y.end") %in% names(df))) {
        c1 <- ll(df$x.start, df$y.start); c2 <- ll(df$x.end, df$y.end)
        df$lon.start <- c1[, 1]; df$lat.start <- c1[, 2]
        df$lon.end <- c2[, 1];   df$lat.end <- c2[, 2]
      }
    }
    df
  })

  output$transect_coords <- renderDT({
    df <- coords()
    num <- vapply(df, is.numeric, logical(1))
    df[num] <- lapply(names(df)[num], function(nm)
      round(df[[nm]], if (grepl("^lon|^lat", nm)) 6 else 2))
    datatable(df, rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE))
  })

  # ---- Coverage simulation -----------------------------------------------
  coverage <- eventReactive(input$run_cov, {
    req(region_res())
    withProgress(message = "Running coverage simulation",
                 detail = paste(input$reps, "replicates"), value = 0.3, {
      safely({
        cg <- dssd::make.coverage(region(), n.grid.points = input$n_grid)
        d <- make_design(coverage.grid = cg)
        set_seed()
        dssd::run.coverage(d, reps = input$reps, quiet = TRUE)
      }, "Coverage simulation")
    })
  })

  observeEvent(coverage(), {
    st <- names(coverage()@design.statistics)
    lab <- c(sampler.count = "Number of samplers", cov.area = "Covered area",
             p.cov.area = "% of region covered", line.length = "Line length",
             trackline = "Trackline length", cyclictrackline = "Cyclic trackline length")
    ch <- stats::setNames(st, ifelse(st %in% names(lab), lab[st], st))
    updateSelectInput(session, "stat", choices = ch, selected = st[1])
  })

  output$cov_map <- renderPlot({ d <- coverage(); req(d); plot(d) })

  cov_scores <- reactive({
    d <- coverage(); req(d)
    s <- d@coverage.scores
    s[!is.na(s)]
  })

  output$cov_hist <- renderPlot({
    s <- cov_scores()
    ggplot(data.frame(s = s), aes(s)) +
      geom_histogram(bins = 30, fill = "#2c3e50", colour = "white") +
      labs(x = "Coverage probability", y = "Grid points") + theme_minimal()
  })

  output$cov_flag <- renderUI({
    s <- cov_scores(); req(length(s))
    tags$p(class = "small text-muted",
           sprintf("Grid points: %d. Coverage probability min %.3f, mean %.3f, max %.3f. ",
                   length(s), min(s), mean(s), max(s)),
           "A wide spread or points at zero indicate non-uniform coverage",
           "(e.g. edge effects or zigzag apex areas), which can bias abundance",
           "estimates if density varies across the region.")
  })

  cov_table <- reactive({
    d <- coverage(); req(d, input$stat)
    m <- d@design.statistics[[input$stat]]; req(m)
    data.frame(Statistic = rownames(m), round(m, 2), check.names = FALSE)
  })

  output$cov_stats <- renderDT(
    datatable(cov_table(), rownames = FALSE, options = list(dom = "t")))

  # ---- R script -----------------------------------------------------------
  r_script <- reactive({
    res <- region_res(); req(res)
    dep <- function(x) paste(deparse(x, width.cutoff = 500L), collapse = "")
    src <- res$src
    read_line <- if (input$region_src == "example") {
      sprintf('shape <- sf::st_read(system.file("extdata", "%s.shp", package = "dssd"), quiet = TRUE)',
              input$example)
    } else {
      sprintf('shape <- sf::st_read("%s"%s, quiet = TRUE)  # edit path',
              basename(src$path),
              if (!is.null(src$layer) && length(src$layers) > 1)
                sprintf(', layer = "%s"', src$layer) else "")
    }
    pc <- res$proj_code
    proj <- c(
      if (!is.null(pc$assign)) sprintf("shape <- sf::st_set_crs(shape, %s)", pc$assign),
      if (!is.null(pc$to)) sprintf("shape <- sf::st_transform(shape, %s)", pc$to)
    )
    strata_lines <- if (res$strata_col == "(none)") {
      sprintf('shape <- sf::st_sf(stratum = "%s", geometry = sf::st_union(sf::st_make_valid(shape)))',
              input$region_name)
    } else if (res$strata_col == "(each)") {
      c("shape <- sf::st_make_valid(shape)",
        'shape <- sf::st_sf(stratum = paste0("S", seq_len(nrow(shape))), geometry = sf::st_geometry(shape))')
    } else {
      c("shape <- sf::st_make_valid(shape)",
        sprintf('grp <- as.character(shape[["%s"]])', res$strata_col),
        "g <- split(sf::st_geometry(shape), factor(grp, levels = unique(grp)))",
        "shape <- sf::st_sf(stratum = names(g), geometry = do.call(c, lapply(g, sf::st_union)))")
    }
    region_line <- sprintf(
      'region <- make.region(region.name = "%s",\n                      strata.name = shape$stratum,\n                      shape = sf::st_cast(shape, "MULTIPOLYGON")%s)',
      input$region_name,
      if (length(res$units)) sprintf(', units = "%s"', res$units) else "")

    out <- c("library(dssd)", "library(sf)", "", "# Region", read_line, proj,
             strata_lines, region_line, "plot(region)")

    if (!is.null(maybe(transects))) {
      a <- transects()$args
      arg_txt <- paste0("  ", names(a), " = ", vapply(a, dep, ""), collapse = ",\n")
      cov_n <- if (!is.null(maybe(coverage))) input$n_grid else 1000
      out <- c(out, "",
               "# Coverage grid (needed for run.coverage)",
               sprintf("cover <- make.coverage(region, n.grid.points = %s)", cov_n),
               "", "# Design",
               paste0("design <- make.design(\n  region = region,\n", arg_txt,
                      ",\n  coverage.grid = cover)"),
               "", "# One realisation of the design",
               if (isTruthy(input$seed)) sprintf("set.seed(%s)", input$seed),
               "transects <- generate.transects(design)",
               "plot(region, transects)",
               "transects",
               "", "# Coverage simulation",
               if (isTruthy(input$seed)) sprintf("set.seed(%s)", input$seed),
               sprintf("design <- run.coverage(design, reps = %s)",
                       if (!is.null(maybe(coverage))) input$reps else 100),
               "plot(design)",
               "design",
               "", "# Export (write.transects() will not overwrite an existing .shp)",
               'write.transects(transects, dsn = "transects.shp")',
               'write.transects(transects, dsn = "transects.csv")',
               if (!is.na(sf::st_crs(res$shape)))
                 sprintf('write.transects(transects, dsn = "transects.gpx", layer = "%s", proj4string = sf::st_crs(4326))',
                         if (input$transect_type == "line") "lines" else "points"))
    }
    paste(out, collapse = "\n")
  })

  output$r_code <- renderText(r_script())

  # ---- Downloads ----------------------------------------------------------
  output$dl_shp <- downloadHandler(
    filename = function() "transects_shp.zip",
    content = function(file) {
      d <- tempfile("shp_"); dir.create(d)
      dssd::write.transects(tr(), dsn = file.path(d, "transects.shp"))
      zip::zipr(file, list.files(d, full.names = TRUE))
    })

  output$dl_gpx <- downloadHandler(
    filename = function() "transects.gpx",
    content = function(file) {
      t <- tr()
      if (is.na(sf::st_crs(t@samplers))) {
        showNotification("No CRS on the region, cannot write GPX.", type = "error")
        stop("No CRS")
      }
      f <- tempfile(fileext = ".gpx")
      dssd::write.transects(t, dsn = f,
                            layer = if (inherits(t, "Point.Transect")) "points" else "lines",
                            proj4string = sf::st_crs(4326))
      file.copy(f, file)
    })

  output$dl_csv <- downloadHandler(
    filename = function() "transects.csv",
    content = function(file) utils::write.csv(coords(), file, row.names = FALSE))

  output$dl_cov <- downloadHandler(
    filename = function() "coverage_statistics.csv",
    content = function(file) {
      d <- coverage(); req(d)
      out <- do.call(rbind, lapply(names(d@design.statistics), function(nm) {
        m <- d@design.statistics[[nm]]
        data.frame(measure = nm, statistic = rownames(m), m, check.names = FALSE)
      }))
      utils::write.csv(out, file, row.names = FALSE)
    })

  output$dl_script <- downloadHandler(
    filename = function() "dssd_design.R",
    content = function(file) writeLines(r_script(), file))
}

shinyApp(ui, server)
