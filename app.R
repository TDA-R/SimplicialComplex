library(Matrix)
library(shiny)
library(gtools)
library(igraph)

source("./R/Faces.R")
source("./R/Betti.R")
source("./R/Boundary.R")
source("./R/EulerCharacteristic.R")
source("./R/AbstractSimplicialComplex.R")
source("./R/VRComplex.R")
source("./R/FloodComplex.R")   # flood_complex, generate_landmarks, flood_persistence, as_filtration

parse_points <- function(txt) {
  lines <- strsplit(txt, "\n")[[1]]
  lines <- trimws(lines)
  lines <- lines[nchar(lines) > 0]
  mat <- do.call(rbind, lapply(lines, function(ln) as.numeric(strsplit(ln, "[,\\s]+")[[1]][1:2])))
  matrix(mat, ncol = 2, byrow = FALSE)
}

# 2D demo point clouds for the Flood mode
generate_cloud <- function(shape, n, noise = 0.05) {
  th <- runif(n, 0, 2 * pi)
  pts <- switch(shape,
                circle = cbind(cos(th), sin(th)),
                two_circles = {half <- n %/% 2 rbind(cbind(cos(th[1:half]) - 1.25, sin(th[1:half])), cbind(cos(th[(half + 1):n]) + 1.25, sin(th[(half + 1):n])))},
                annulus = {r <- sqrt(runif(n, 0.5^2, 1)) cbind(r * cos(th), r * sin(th))}
  )
  pts + matrix(rnorm(2 * length(pts[, 1]), sd = noise), ncol = 2)
}

# static explainer figure for the Flood Complex guide tab
flood_guide_svg <- '
<svg viewBox="0 0 800 330" xmlns="http://www.w3.org/2000/svg" role="img"
     font-family="system-ui, sans-serif" style="max-width:820px;width:100%;height:auto;">
  <title>Flood complex: from Delaunay to grid points</title>
  <text x="130" y="24" text-anchor="middle" font-size="13" font-weight="600" fill="#333">1. Point cloud + FPS landmarks</text>
  <g fill="#999" opacity="0.55">
    <circle cx="55" cy="70" r="2"/><circle cx="72" cy="55" r="2"/><circle cx="95" cy="48" r="2"/><circle cx="120" cy="45" r="2"/><circle cx="148" cy="52" r="2"/><circle cx="170" cy="65" r="2"/><circle cx="185" cy="88" r="2"/><circle cx="190" cy="115" r="2"/><circle cx="183" cy="142" r="2"/><circle cx="165" cy="163" r="2"/><circle cx="140" cy="175" r="2"/><circle cx="112" cy="178" r="2"/><circle cx="85" cy="170" r="2"/><circle cx="63" cy="152" r="2"/><circle cx="50" cy="128" r="2"/><circle cx="48" cy="98" r="2"/>
    <circle cx="66" cy="82" r="2"/><circle cx="88" cy="60" r="2"/><circle cx="135" cy="48" r="2"/><circle cx="160" cy="57" r="2"/><circle cx="180" cy="76" r="2"/><circle cx="189" cy="102" r="2"/><circle cx="188" cy="130" r="2"/><circle cx="174" cy="154" r="2"/><circle cx="152" cy="170" r="2"/><circle cx="126" cy="178" r="2"/><circle cx="98" cy="175" r="2"/><circle cx="73" cy="162" r="2"/><circle cx="55" cy="140" r="2"/><circle cx="47" cy="112" r="2"/><circle cx="57" cy="86" r="2"/><circle cx="79" cy="66" r="2"/>
  </g>
  <g fill="#c76b1d">
    <circle cx="120" cy="45" r="5"/><circle cx="190" cy="115" r="5"/><circle cx="112" cy="178" r="5"/><circle cx="48" cy="98" r="5"/><circle cx="165" cy="163" r="5"/><circle cx="72" cy="55" r="5"/>
  </g>
  <text x="130" y="215" text-anchor="middle" font-size="11" fill="#777">grey = data points (flood sources)</text>
  <text x="130" y="231" text-anchor="middle" font-size="11" fill="#777">orange = landmarks</text>
  <path d="M 225 115 L 255 115" stroke="#999" stroke-width="1.5" fill="none" marker-end="url(#arr)"/>
  <text x="390" y="24" text-anchor="middle" font-size="13" font-weight="600" fill="#333">2. Delaunay on landmarks</text>
  <g fill="#c76b1d" fill-opacity="0.12" stroke="#c76b1d" stroke-width="1.5">
    <polygon points="380,45 332,55 308,98"/>
    <polygon points="380,45 308,98 450,115"/>
    <polygon points="380,45 450,115 332,55" fill-opacity="0"/>
    <polygon points="308,98 372,178 450,115"/>
    <polygon points="450,115 372,178 425,163"/>
  </g>
  <g fill="#c76b1d">
    <circle cx="380" cy="45" r="5"/><circle cx="450" cy="115" r="5"/><circle cx="372" cy="178" r="5"/><circle cx="308" cy="98" r="5"/><circle cx="425" cy="163" r="5"/><circle cx="332" cy="55" r="5"/>
  </g>
  <text x="390" y="215" text-anchor="middle" font-size="11" fill="#777">vertex = 0-simplex, edge = 1-simplex</text>
  <text x="390" y="231" text-anchor="middle" font-size="11" fill="#777">triangle = 2-simplex (3D adds tetrahedra)</text>
  <path d="M 480 115 L 510 115" stroke="#999" stroke-width="1.5" fill="none" marker-end="url(#arr)"/>
  <text x="655" y="24" text-anchor="middle" font-size="13" font-weight="600" fill="#333">3. Grid on one simplex</text>
  <polygon points="560,60 750,120 620,240" fill="#c76b1d" fill-opacity="0.08" stroke="#c76b1d" stroke-width="1.5"/>
  <g fill="#444">
    <circle cx="620" cy="240" r="3"/><circle cx="652.5" cy="210" r="3"/><circle cx="685" cy="180" r="3"/><circle cx="717.5" cy="150" r="3"/><circle cx="750" cy="120" r="3"/><circle cx="605" cy="195" r="3"/><circle cx="637.5" cy="165" r="3"/><circle cx="670" cy="135" r="3"/><circle cx="702.5" cy="105" r="3"/><circle cx="590" cy="150" r="3"/><circle cx="622.5" cy="120" r="3"/><circle cx="655" cy="90" r="3"/><circle cx="575" cy="105" r="3"/><circle cx="607.5" cy="75" r="3"/><circle cx="560" cy="60" r="3"/>
  </g>
  <g fill="#c76b1d">
    <circle cx="560" cy="60" r="5"/><circle cx="750" cy="120" r="5"/><circle cx="620" cy="240" r="5"/>
  </g>
  <circle cx="700" cy="205" r="2.5" fill="#999"/>
  <circle cx="580" cy="130" r="2.5" fill="#999"/>
  <path d="M 670 135 L 698 202" stroke="#999" stroke-width="1" stroke-dasharray="3 3" fill="none"/>
  <path d="M 590 150 L 582 133" stroke="#999" stroke-width="1" stroke-dasharray="3 3" fill="none"/>
  <text x="655" y="275" text-anchor="middle" font-size="11" fill="#777">each grid point queries its nearest data point (dashed)</text>
  <text x="655" y="291" text-anchor="middle" font-size="11" fill="#777">filtration value = the maximum of those distances</text>
  <defs>
    <marker id="arr" markerWidth="8" markerHeight="8" refX="7" refY="4" orient="auto">
      <path d="M0,0 L8,4 L0,8 z" fill="#999"/>
    </marker>
  </defs>
</svg>'

default_simplices <- list(c(1, 2), c(3, 4), c(2, 1, 3), c(4, 2))
default_input <- paste0(
  "list(\n",
  paste(sapply(default_simplices, function(x) paste0("  c(", paste(x, collapse = ", "), ")")), collapse = ",\n"),
  "\n)"
)
default_points_txt <- "0,0
1,0
1,1
0,1"

ui <- fluidPage(
  titlePanel("Simplicial Complex Explorer"),

  sidebarLayout(
    sidebarPanel(
      radioButtons(
        "mode", "Input Mode",
        choices = c("Simplices (manual)" = "simp",
                    "Data points → VR simplices" = "pts",
                    "Data points → Flood complex" = "flood"),
        selected = "simp"
      ),

      conditionalPanel(
        condition = "input.mode == 'simp'",
        textAreaInput("simplices_input", "Input Simplices", value = default_input, rows = 10, width = "100%")
      ),

      conditionalPanel(
        condition = "input.mode == 'pts'",
        textAreaInput("points_input", "Points (each line: x,y)", value = default_points_txt, rows = 6, width = "100%"),
        numericInput("epsilon_pts", "ε (VR edge threshold)", value = 1.5, step = 0.1, min = 0)
      ),

      conditionalPanel(
        condition = "input.mode == 'flood'",
        selectInput("flood_shape", "Point cloud",
                    choices = c("Circle" = "circle", "Two circles" = "two_circles", "Annulus" = "annulus", "Custom (textarea)" = "custom")),
        conditionalPanel(
          condition = "input.flood_shape == 'custom'",
          textAreaInput("flood_points_input", "Points (each line: x,y)", rows = 6, width = "100%")
        ),
        numericInput("flood_n", "Number of witness points", value = 600, min = 50, step = 50),
        numericInput("flood_lms", "Number of landmarks (FPS)", value = 30, min = 4, step = 1),
        numericInput("flood_ppe", "Grid points per edge", value = 10, min = 2, max = 30),
        sliderInput("flood_t", "Flood threshold t", min = 0, max = 1, value = 0.15, step = 0.005),
        checkboxInput("flood_balls", "Show flood balls (radius t)", value = FALSE)
      ),

      tags$hr(),
      numericInput("dim", "Target Dimension k", value = 0, min = 0, max = 10),
      numericInput("tol", "tolerance (for Betti & Euler)", value = 0.1, step = 0.05),
      actionButton("compute", "Compute", class = "btn-primary", width = "100%"),
      width = 3
    ),

    mainPanel(
      conditionalPanel(
        condition = "input.mode == 'pts'",
        h4("Plot"),
        plotOutput("vr_plot"),
        tags$hr()
      ),
      conditionalPanel(
        condition = "input.mode == 'flood'",
        h4("Flood complex at threshold t"),
        plotOutput("flood_plot", height = "420px"),
        tags$hr()
      ),
      tabsetPanel(
        tabPanel("Faces",
                 h4("Faces Interpretation"),
                 withMathJax(HTML("simplicial complex outputs the faces of dimension \\( k \\)(points, edge,...)")),
                 verbatimTextOutput("faces_out")
        ),

        tabPanel("Boundary Matrix",
                 h4("Boundary Matrix \\( \\partial_k \\)"),
                 verbatimTextOutput("boundary_out")
        ),

        tabPanel("Betti Numbers",
                 h4("Betti Number \\( \\beta_k \\)"),
                 verbatimTextOutput("betti_out")
        ),

        tabPanel("Euler Characteristic",
                 h4("Euler Characteristic \\( \\chi \\)"),
                 verbatimTextOutput("euler_out")
        ),

        tabPanel("Abstract Simplicial Complex",
                 h4("Abstract Simplicial Complex"),
                 verbatimTextOutput("abstract_out")
        ),

        tabPanel("Flood Complex Guide",
                 h4("How the Flood Complex works"),
                 withMathJax(
                   HTML("The Flood complex (Graf et al., NeurIPS 2025) makes
                        persistent homology feasible on large point clouds by
                        separating two roles: a small set of <b>landmarks</b>
                        determines <i>which simplices can exist</i>, while the
                        <b>full point cloud</b> determines <i>when each simplex
                        appears</i>.")
                 ),
                 HTML(flood_guide_svg),
                 withMathJax(
                   HTML("<ol>
                        <li><b>Landmarks (FPS).</b> Pick a small number of
                        well-spread representatives via Farthest-Point Sampling:
                        repeatedly add the point farthest from those already
                        chosen.</li>
                        <li><b>Delaunay = the candidate simplices.</b> The
                        Delaunay triangulation of the landmarks <i>is</i> a
                        simplicial complex: each triangle (tetrahedron in 3D)
                        is a top simplex, and its vertex subsets are the lower
                        faces. No extra geometry is needed &mdash; only subset
                        enumeration. The complex size depends only on the
                        number of landmarks.</li>
                        <li><b>Flood time via a barycentric grid.</b> Imagine
                        water rising from every data point: at time \\( t \\)
                        each point carries a ball of radius \\( t \\). A simplex
                        \\( \\sigma \\) enters the complex once it is fully
                        submerged. To estimate that time, sample \\( \\sigma \\)
                        with a barycentric grid (weights \\( w_i \\ge 0, \\sum w_i = 1 \\),
                        discretized with <i>points_per_edge</i> steps) and set
                        $$ f(\\sigma) \\;=\\; \\max_{g \\in \\text{grid}(\\sigma)} \\; \\min_{x \\in X} \\; \\lVert g - x \\rVert $$
                        &mdash; the hardest-to-reach spot on the simplex decides
                        its flood time. Grid points are probes only; they are
                        discarded afterwards.</li>
                        <li><b>Persistence.</b> Sort simplices by flood time and
                        run the usual boundary-matrix reduction. Simplices that
                        stretch across empty regions (holes) flood late, so the
                        holes survive long and appear far from the diagonal in
                        the persistence diagram.</li>
                        </ol>
                        <b>Try it:</b> switch the input mode to
                        <i>Data points &rarr; Flood complex</i>, press Compute,
                        then drag the threshold slider to watch the complex
                        flood in. The tabs above (Faces, Boundary, Betti, Euler)
                        always describe the sub-complex at the current
                        threshold, and the Persistence tab shows which classes
                        are alive.")
                 )
        ),

        tabPanel("Persistence (Flood)",
                 h4("Persistence Diagram"),
                 withMathJax(
                   HTML("Each point is a homology class: born at x, dead at y.
                        Classes far above the diagonal are true topological
                        features; classes near it are noise. The dashed cross
                        marks the current threshold \\( t \\): classes with
                        birth \\( \\le t < \\) death are alive in the complex
                        shown above.")
                 ),
                 plotOutput("pd_plot", height = "420px"),
                 verbatimTextOutput("pd_summary")
        )
      ),
      width = 9
    )
  )
)

server <- function(input, output, session) {

  # Flood complex: computed once per Compute click
  flood_state <- eventReactive(input$compute, {
    req(input$mode == "flood")
    pts <- if (input$flood_shape == "custom") {
      parse_points(input$flood_points_input)
    } else {
      generate_cloud(input$flood_shape, input$flood_n)
    }
    validate(need(nrow(pts) >= 4, "Need at least 4 points"))
    n_lms <- min(input$flood_lms, nrow(pts))
    fc <- flood_complex(pts, landmarks = n_lms,
                        points_per_edge = input$flood_ppe, backend = "cpu")
    pairs <- flood_persistence(as_filtration(fc))
    t_max <- max(fc$filtration)
    updateSliderInput(session, "flood_t",
                      max = round(t_max * 1.05, 3),
                      value = round(t_max * 0.3, 3))
    list(points = pts, fc = fc, pairs = pairs)
  }, ignoreInit = TRUE)

  # simplices fed into the shared tabs
  # VR / manual: fixed on Compute. Flood: also reacts to the t slider,
  # so the tabs always describe the sub-complex currently plotted.
  simplices_current <- reactive({
    if (input$mode == "flood") {
      st <- flood_state()
      keep <- st$fc$filtration <= input$flood_t
      validate(need(any(keep), "No simplex alive at this threshold yet"))
      st$fc$simplices[keep]
    } else {
      input$compute
      isolate({
        if (input$mode == "simp") {
          eval(parse(text = input$simplices_input))
        } else {
          pts <- parse_points(input$points_input)
          VietorisRipsComplex(pts, input$epsilon_pts)$simplices
        }
      })
    }
  })

  # Flood plot
  output$flood_plot <- renderPlot({
    st <- flood_state()
    t <- input$flood_t
    pts <- st$points
    lms <- st$fc$landmarks
    filt <- st$fc$filtration
    simp <- st$fc$simplices
    dims <- lengths(simp) - 1L
    alive <- filt <= t

    plot(pts, pch = 16, cex = 0.35, col = "grey65", asp = 1,
         xlab = "", ylab = "",
         main = sprintf("t = %.3f | vertices %d, edges %d, triangles %d",
                        t, sum(alive & dims == 0), sum(alive & dims == 1),
                        sum(alive & dims == 2)))

    if (isTRUE(input$flood_balls)) {
      symbols(pts[, 1], pts[, 2], circles = rep(t, nrow(pts)),
              inches = FALSE, add = TRUE, fg = NA,
              bg = rgb(0.2, 0.55, 0.9, alpha = 0.04))
    }
    for (i in which(alive & dims == 2)) {
      polygon(lms[simp[[i]], 1], lms[simp[[i]], 2],
              col = rgb(0.2, 0.55, 0.9, alpha = 0.35), border = NA)
    }
    for (i in which(alive & dims == 1)) {
      lines(lms[simp[[i]], 1], lms[simp[[i]], 2], col = "steelblue4", lwd = 1.6)
    }
    points(lms, pch = 21, bg = "orange", cex = 1.4)
    text(lms, labels = seq_len(nrow(lms)), pos = 3, cex = 0.7)
  })

  # Persistence diagram
  output$pd_plot <- renderPlot({
    st <- flood_state()
    pairs <- st$pairs
    t <- input$flood_t
    fin <- pairs[is.finite(pairs$death), ]
    lim <- c(0, max(fin$death, t) * 1.05)
    alive <- fin$birth <= t & fin$death > t
    plot(fin$birth, fin$death, pch = 19, cex = ifelse(alive, 1.2, 0.7),
         col = c("black", "red", "blue")[fin$dim + 1],
         xlim = lim, ylim = lim, xlab = "Birth", ylab = "Death")
    abline(0, 1, lty = 2, col = "grey")
    abline(v = t, h = t, lty = 3, col = "darkgreen")
    legend("bottomright", legend = paste0("H", 0:2),
           col = c("black", "red", "blue"), pch = 19)
  })

  output$pd_summary <- renderPrint({
    st <- flood_state()
    pairs <- st$pairs
    t <- input$flood_t
    alive <- (pairs$birth <= t) & (pairs$death > t)
    b <- table(factor(pairs$dim[alive], levels = 0:2))
    cat(sprintf("Alive at t = %.3f:  β₀ = %d, β₁ = %d, β₂ = %d\n",
                t, b[1], b[2], b[3]))
    cat(sprintf("(essential classes: %d)\n", sum(is.infinite(pairs$death))))
  })

  # shared tabs (all modes)
  output$faces_out <- renderPrint({
    simplices <- simplices_current()
    cat("Faces in dim", input$dim, ":\n")
    print(faces(simplices, target_dim = input$dim))
  })

  output$boundary_out <- renderPrint({
    simplices <- simplices_current()
    cat("Boundary matrix for ∂", input$dim, ":\n")
    print(boundary(simplices, input$dim))
  })

  output$betti_out <- renderPrint({
    simplices <- simplices_current()
    res <- sapply(0:(input$dim + 1), function(d) {
      paste0("β_", d, " = ", betti_number(simplices, d, input$tol))
    })
    cat(paste(res, collapse = "\n"))
  })

  output$euler_out <- renderPrint({
    simplices <- simplices_current()
    cat("Euler Characteristic (χ):", euler_characteristic(simplices, input$tol))
  })

  output$abstract_out <- renderPrint({
    simplices <- simplices_current()
    cat("Abstract Simplicial Complex (dim", input$dim, "):\n")
    print(abstract_simplicial_complex(simplices, input$dim, input$tol))
  })

  # original VR plot (pad_y bug fixed: was undefined)
  output$vr_plot <- renderPlot({
    req(input$mode == "pts")
    pts <- parse_points(input$points_input)
    vr <- VietorisRipsComplex(pts, input$epsilon_pts)
    g <- vr$network
    req(!is.null(g))
    xr <- range(pts[, 1]); yr <- range(pts[, 2])
    pad_y <- 0.08 * diff(yr)
    par(xpd = NA)
    plot(
      g,
      layout = as.matrix(pts),
      vertex.label = 1:nrow(pts),
      vertex.size = 12,
      edge.arrow.mode = 0,
      rescale = FALSE,
      xlim = xr,
      ylim = c(yr[1] - pad_y, yr[2] + pad_y),
      asp = 0
    )
  })
}

shinyApp(ui, server)

