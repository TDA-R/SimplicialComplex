library(shiny)
library(Matrix)
library(gtools)

source("R/faces.R")
source("R/Betti.R")
source("R/Boundary.R")
source("R/EulerCharacteristic.R")
source("R/AbstractSimplicialComplex.R")

# default simplices
default_simplices <- list(c(2, 1, 3), c(4, 2), c(5), c(2, 3, 5, 4))

# text format for user input
default_input <- paste0(
  "list(\n",
  paste(sapply(default_simplices, function(x) paste0("  c(", paste(x, collapse=", "), ")")), collapse = ",\n"),
  "\n)"
)

ui <- fluidPage(
  titlePanel("Simplicial Complex Explorer"),

  sidebarLayout(
    sidebarPanel(
      textAreaInput("simplices_input", "Input Simplices", value = default_input, rows = 10),
      numericInput("dim", "Target Dimension k", value = 0, min = 0, max = 10),
      numericInput("eps", "Epsilon (for Betti & Euler)", value = 0.1, step = 0.05),
      actionButton("compute", "Compute")
    ),

    mainPanel(
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
                   HTML("$$\\partial_k \\sigma = \\sum_i (-1)^i [v_0\\ v_1\\ \\ldots\\ \\hat{v}_i\\ \\ldots\\ v_k]$$")
                 ),
                 verbatimTextOutput("boundary_out")
        ),

        tabPanel("Betti Numbers",
                 h4("Betti Number \\( \\beta_k \\)"),
                 withMathJax(
                   HTML("Betti number represent the number of holes in k dimension. \\( \\beta_0 \\) is the number of connected components, \\( \\beta_1 \\) is the number of 1-dimensional loops, and so on.")
                 ),
                 verbatimTextOutput("betti_out")
        ),

        tabPanel("Euler Characteristic",
                 h4("Euler Characteristic \\( \\chi \\)"),
                 withMathJax(
                   HTML("Euler formula :\\[ \\chi = \\sum_{k=0}^{n} (-1)^k \\beta_k \\]")
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

server <- function(input, output, session) {

  observeEvent(input$compute, {
    tryCatch({
      simplices <- eval(parse(text=input$simplices_input))
      k <- input$dim
      eps <- input$eps

      # faces
      output$faces_out <- renderPrint({
        cat("Faces in dim", k, ":\n")
        print(faces(simplices, target_dim=k))
      })

      # boundary
      output$boundary_out <- renderPrint({
        cat("Boundary matrix for ∂", k, ":\n")
        print(boundary(simplices, k))
      })

      # betti
      output$betti_out <- renderPrint({
        res <- sapply(0:(k+1), function(d) {
          paste0("β_", d, " = ", betti_number(simplices, d, eps))
        })
        cat(paste(res, collapse = "\n"))
      })

      # euler characteristic
      output$euler_out <- renderPrint({
        χ <- euler_characteristic(simplices, eps)
        cat("Euler Characteristic (χ):", χ)
      })

      # abstract simplicial complex
      output$abstract_out <- renderPrint({
        cat("Abstract Simplicial Complex (dim", k, "):\n")
        print(abstract_simplicial_complex(simplices, k, eps))
      })

    }, error = function(e) {
      showNotification(paste("Error:", e$message), type="error")
    })
  })
}

shinyApp(ui, server)
