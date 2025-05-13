library(shiny)
library(shinythemes)
library(heatmaply)
library(plotly)
library(ggplot2)
library(SummarizedExperiment)
library(DT)
library(limma)

ui <- fluidPage(
  theme = shinytheme("flatly"),
  titlePanel("TCGA Explorer: Interactive Genomic Analysis"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("rse_file", "Upload SummarizedExperiment (.rds):"),
      uiOutput("subtype_col_ui"),
      uiOutput("subtype_ui"),
      selectizeInput("genes", "Select Genes (Optional):", choices = NULL, multiple = TRUE),
      downloadButton("downloadPlot", "Download Heatmap (PNG)"),
      helpText("App by Renuka Reddy Namala")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Heatmap", plotlyOutput("heatmap", height = "700px")),
        tabPanel("PCA Plot", plotOutput("pcaPlot", height = "600px")),
        tabPanel("DE Analysis", dataTableOutput("deTable")),
        tabPanel("Volcano Plot", plotOutput("volcanoPlot", height = "600px"))
      )
    )
  )
)

server <- function(input, output, session) {
  rse <- reactiveVal(NULL)
  z_matrix <- reactiveVal(NULL)
  subtype_choices <- reactiveVal("All")
  subtype_colname <- reactiveVal(NULL)
  
  observeEvent(input$rse_file, {
    req(input$rse_file)
    rse_data <- readRDS(input$rse_file$datapath)
    rse(rse_data)
    
    expr <- assay(rse_data)
    expr_z <- t(scale(t(expr)))
    z_matrix(expr_z)
    
    updateSelectizeInput(session, "genes", choices = rownames(expr_z))
  })
  
  output$subtype_col_ui <- renderUI({
    req(rse())
    selectInput("subtype_col", "Select Subtype Column:", choices = c("None", colnames(colData(rse()))))
  })
  
  observeEvent(input$subtype_col, {
    req(input$subtype_col, rse())
    if (input$subtype_col != "None") {
      subtypes <- colData(rse())[[input$subtype_col]]
      subtype_choices(c("All", sort(unique(as.character(subtypes)))))
      subtype_colname(input$subtype_col)
    } else {
      subtype_choices("All")
      subtype_colname(NULL)
    }
  })
  
  output$subtype_ui <- renderUI({
    selectInput("subtype", "Filter by Subtype:", choices = subtype_choices())
  })
  
  filtered_matrix <- reactive({
    req(z_matrix())
    mat <- z_matrix()
    samples <- colnames(mat)
    
    if (!is.null(input$subtype) && input$subtype != "All" && !is.null(rse()) && !is.null(subtype_colname())) {
      subtypes <- colData(rse())[[subtype_colname()]]
      selected_samples <- colnames(rse())[which(as.character(subtypes) == input$subtype)]
      samples <- intersect(samples, selected_samples)
      if (length(samples) == 0) return(NULL)
      mat <- mat[, samples, drop = FALSE]
    }
    
    if (length(input$genes) > 0) {
      mat <- mat[rownames(mat) %in% input$genes, , drop = FALSE]
    }
    if (ncol(mat) == 0 || nrow(mat) == 0) return(NULL)
    mat
  })
  
  output$heatmap <- renderPlotly({
    mat <- filtered_matrix()
    req(!is.null(mat))
    req(nrow(mat) >= 2, ncol(mat) >= 2)
    
    heatmaply(
      mat,
      xlab = "Samples",
      ylab = "Genes",
      main = paste("Heatmap -", input$subtype),
      scale_fill_gradient_fun = scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0),
      fontsize_row = 9,
      fontsize_col = 8,
      showticklabels = c(TRUE, FALSE),
      margins = c(80, 120, 60, 20),
      key.title = "Z-Score"
    )
  })
  
  output$downloadPlot <- downloadHandler(
    filename = function() paste0("heatmap_", input$subtype, ".png"),
    content = function(file) {
      png(file, width = 1200, height = 900)
      mat <- filtered_matrix()
      if (is.null(mat)) {
        dev.off()
        return()
      }
      heatmaply(
        mat,
        xlab = "Samples",
        ylab = "Genes",
        main = paste("Heatmap -", input$subtype),
        scale_fill_gradient_fun = scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0),
        fontsize_row = 9,
        fontsize_col = 8,
        showticklabels = c(TRUE, FALSE),
        margins = c(80, 120, 60, 20),
        key.title = "Z-Score"
      )
      dev.off()
    }
  )
  
  output$pcaPlot <- renderPlot({
    mat <- filtered_matrix()
    req(!is.null(mat))
    req(ncol(mat) >= 2)
    subtypes <- if (!is.null(subtype_colname())) colData(rse())[[subtype_colname()]] else factor(rep("All", ncol(mat)))
    
    pca <- prcomp(t(mat), scale. = TRUE)
    df <- data.frame(
      PC1 = pca$x[, 1],
      PC2 = pca$x[, 2],
      Subtype = factor(subtypes[colnames(mat)])
    )
    
    ggplot(df, aes(x = PC1, y = PC2, color = Subtype)) +
      geom_point(size = 3, alpha = 0.8) +
      theme_minimal() +
      labs(title = "PCA of Samples", x = "PC1", y = "PC2")
  })
  
  output$deTable <- renderDataTable({
    req(rse(), subtype_colname())
    subtypes <- colData(rse())[[subtype_colname()]]
    if (length(unique(subtypes)) < 2) return(NULL)
    
    design <- model.matrix(~ 0 + factor(subtypes))
    colnames(design) <- levels(factor(subtypes))
    fit <- lmFit(assay(rse()), design)
    cont.matrix <- makeContrasts(contrasts = paste(colnames(design)[1], "-", colnames(design)[2]), levels = design)
    fit2 <- contrasts.fit(fit, cont.matrix)
    fit2 <- eBayes(fit2)
    topTable(fit2, number = 50)
  })
  
  output$volcanoPlot <- renderPlot({
    req(rse(), subtype_colname())
    subtypes <- colData(rse())[[subtype_colname()]]
    validate(need(length(unique(subtypes)) >= 2, "Need at least two groups"))
    
    design <- model.matrix(~ 0 + factor(subtypes))
    colnames(design) <- levels(factor(subtypes))
    fit <- lmFit(assay(rse()), design)
    cont.matrix <- makeContrasts(contrasts = paste(colnames(design)[1], "-", colnames(design)[2]), levels = design)
    fit2 <- contrasts.fit(fit, cont.matrix)
    fit2 <- eBayes(fit2)
    tt <- topTable(fit2, number = Inf, adjust = "BH")
    
    tt$Significance <- "NS"
    tt$Significance[tt$logFC > 1 & tt$adj.P.Val < 0.05] <- "Up"
    tt$Significance[tt$logFC < -1 & tt$adj.P.Val < 0.05] <- "Down"
    
    ggplot(tt, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
      geom_point(alpha = 0.6) +
      scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "gray")) +
      theme_minimal() +
      labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-log10 Adjusted P-Value")
  })
}

shinyApp(ui = ui, server = server)
