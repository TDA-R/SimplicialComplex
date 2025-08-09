library(Matrix)
library(shiny)
library(gtools)
library(igraph)

source("R/faces.R")
source("R/Betti.R")
source("R/Boundary.R")
source("R/EulerCharacteristic.R")
source("R/AbstractSimplicialComplex.R")

source("R/VRcomplex.R")

parse_points <- function(txt) {
  lines <- strsplit(txt, "\n")[[1]]
  lines <- trimws(lines)
  lines <- lines[nchar(lines) > 0]
  mat <- do.call(rbind, lapply(lines, function(ln) as.numeric(strsplit(ln, "[,\\s]+")[[1]][1:2])))
  matrix(mat, ncol = 2, byrow = FALSE)
}

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
        choices = c("Simplices (manual)" = "simp", "Data points → VR simplices" = "pts"),
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

      tags$hr(),
      numericInput("dim", "Target Dimension k", value = 0, min = 0, max = 10),
      numericInput("tol", "tolerance (for Betti & Euler)", value = 0.1, step = 0.05),
      actionButton("compute", "Compute", class = "btn-primary", width = "100%")
    ),

    mainPanel(
      conditionalPanel(
        condition = "input.mode == 'pts'",
        h4("Plot"),
        plotOutput("vr_plot", height = "300px"),
        tags$hr()
      ),
      tabsetPanel(
        tabPanel("Faces",
                 h4("Faces Interpretation"),
                 withMathJax(
                   HTML("simplicial complex outputs the faces of dimension \\( k \\)(points, edge,...)")
                 ),
                 verbatimTextOutput("faces_out")
        ),

        tabPanel("Boundary Matrix",
                 h4("Boundary Matrix \\( \\partial_k \\)"),
                 withMathJax(
                   HTML("$$\\partial_k \\sigma = \\sum_i (-1)^i [v_0\\ v_1\\ \\ldots\\ \\hat{v}_i\\ \\ldots\\ v_k]$$ <br/>
                        Boundary function returns a matrix, the shape is the dimension of the faces and it's lower one.
                        For example, if you compute 1-dimensional boundary with the data above, you'll get a 4 (lower dimensional face)
                        x 5 (dimensiional face) matrix. <br/>
                        Noew you get [1,2,3,4] and [(1,2), (3,4), (1,3), (2,3), (2,4)], which you can then get whether there's a
                        pair for 1 or -1 by using the formula above, or just empty.
                        In detail, the formula will be $${(-1)^1[1, (1,2)], (-1)^2[2, (1,2)], (-1)^3[3, (3,4)], (-1)^4[4, (3,4)],
                        (-1)^5[1, (1,3)], (-1)^6[3, (1,3)], (-1)^7[2, (2,3)], (-1)^8[3, (2,3)], (-1)^9[2, (2,4)], (-1)^10[4, (2,4)]}$$")
                 ),
                 verbatimTextOutput("boundary_out"),
                 withMathJax(
                   HTML("The reason of using the formula is that it gives us the orientation of the face.<br/>
                        So the main reason of 1 and -1 is to let boundaries of boundaries equal zero :
                        $$\\partial_{k-1}\\partial_k=0$$<br/>
                        ")
                 ),
        ),

        tabPanel("Betti Numbers",
                 h4("Betti Number \\( \\beta_k \\)"),
                 withMathJax(
                   HTML("Betti number represent the number of holes in k dimension.
                        \\( \\beta_0 \\) is the number of connected components,
                        \\( \\beta_1 \\) is the number of 1-dimensional loops, and so on.<br/>
                        The formula of betti number is : <br/>
                        \\(\\beta_k = \\text{rank} (\\text{ker} \\partial_k) - \\text{rank} (\\text{Im} \\partial_{k+1})\\),
                        then by replacing \\(\\text{ker} \\partial_k\\) with rank-nullity formula:<br/>
                        \\(\\text{rank} (\\text{ker} \\partial_k) +\\text{rank} (\\text{Im} \\partial_k) = \\text{dim} (\\partial_k)\\), you get:<br/>
                        $$\\beta_k = \\text{dim} (\\partial_k) - \\text{rank}(\\text{Im} \\partial_k)) - \\text{rank} (\\text{Im} \\partial_{k+1})$$<br/>
                        ker: How many closed structure, all k-cycles (no obundary k-chains)<br/>
                        Im: Which closed structures are the boundaries of high-dimensional objects. (k+1 dimensional boundaries)<br/>
                        Betti: The number of no boundaries k-chain, minus those that are actually some higher-dimensional boundary (quotient space)")
                 ),
                 verbatimTextOutput("betti_out")
        ),

        tabPanel("Euler Characteristic",
                 h4("Euler Characteristic \\( \\chi \\)"),
                 withMathJax(
                   HTML("Euler formula :\\[ \\chi = \\sum_{k=0}^{n} (-1)^k \\beta_k \\] <br/>
                        It returns the sum of the Betti numbers with alternating signs.
                        Different shapes could have the same Euler characteristic,
                        but it gives a global property of the simplicial complex.")
                 ),
                 verbatimTextOutput("euler_out")
        ),

        tabPanel("Abstract Simplicial Complex",
                 h4("Abstract Simplicial Complex"),
                 withMathJax(
                   HTML("Give a summary of the global structure of simplicial complex")
                 ),
                 verbatimTextOutput("abstract_out")
        )
      )
    )
  )
)

server <- function(
    input, output, session
) {
  simplices_reactive <- eventReactive(input$compute, {
    if (input$mode == "simp") {
      eval(parse(text = input$simplices_input))
    } else {
      pts <- parse_points(input$points_input)
      vr  <- VietorisRipsComplex(pts, input$epsilon_pts)
      vr$simplices
    }
  }, ignoreInit = TRUE)

  observeEvent(input$compute, {
    tryCatch({
      simplices <- simplices_reactive()
      k <- input$dim
      tol <- input$tol

      if (input$mode == "pts") {
        output$simplices_generated <- renderPrint({
          print(simplices)
        })
      }

      output$faces_out <- renderPrint({
        cat("Faces in dim", k, ":\n")
        print(faces(simplices, target_dim = k))
      })

      output$boundary_out <- renderPrint({
        cat("Boundary matrix for ∂", k, ":\n")
        print(boundary(simplices, k))
      })

      output$betti_out <- renderPrint({
        res <- sapply(0:(k + 1), function(d) {
          paste0("β_", d, " = ", betti_number(simplices, d, tol))
        })
        cat(paste(res, collapse = "\n"))
      })

      output$euler_out <- renderPrint({
        chi <- euler_characteristic(simplices, tol)
        cat("Euler Characteristic (χ):", chi)
      })

      output$abstract_out <- renderPrint({
        cat("Abstract Simplicial Complex (dim", k, "):\n")
        print(abstract_simplicial_complex(simplices, k, tol))
      })
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  })

  output$vr_plot <- renderPlot({
    req(input$mode == "pts")
    pts <- parse_points(input$points_input)
    vr  <- VietorisRipsComplex(pts, input$epsilon_pts)
    g <- vr$network
    req(!is.null(g))
    plot(
      g,
      layout = as.matrix(pts),
      vertex.label = 1:nrow(pts),
      vertex.size = 12,
      edge.arrow.mode = 0,
      asp = 1
    )
  })
}

shinyApp(ui, server)
