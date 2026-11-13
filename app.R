library(shiny)
library(bslib)
library(dplyr)
library(ggplot2)
library(DT)

# ------------------------------------------------------------------------------
# 1. Load Processed Multi-Omics Data (with fallback demo generator)
# ------------------------------------------------------------------------------
load_omics_data <- function() {
  if (file.exists("data/processed_omics.rds")) {
    return(readRDS("data/processed_omics.rds"))
  }
  
  set.seed(42)
  n_genes <- 1500
  genes <- paste0("Gene_", sprintf("%04d", 1:n_genes))
  genes[1:6] <- c("Fev", "Lmx1b", "Tph2", "Slc6a4", "Ddc", "Gata2")
  
  de_df <- tibble(
    gene = genes,
    rna_log2FC = rnorm(n_genes, 0, 1.3),
    rna_padj = 10^(-runif(n_genes, 0.1, 8)),
    atac_log2FC = rnorm(n_genes, 0, 1.1),
    atac_padj = 10^(-runif(n_genes, 0.1, 6)),
    tf_motif = sample(c("Pet-1 / FEV", "LMX1B", "GATA2/3", "REST", "None"), n_genes, replace = TRUE)
  ) %>%
    mutate(
      rna_log2FC = ifelse(gene %in% c("Fev", "Lmx1b", "Tph2", "Slc6a4"), abs(rna_log2FC) + 1.8, rna_log2FC),
      rna_padj = ifelse(gene %in% c("Fev", "Lmx1b", "Tph2", "Slc6a4"), 1e-12, rna_padj),
      atac_log2FC = ifelse(gene %in% c("Fev", "Lmx1b", "Tph2", "Slc6a4"), abs(atac_log2FC) + 1.5, atac_log2FC),
      status = case_when(
        rna_padj < 0.05 & rna_log2FC > 1  ~ "Up-regulated",
        rna_padj < 0.05 & rna_log2FC < -1 ~ "Down-regulated",
        TRUE ~ "Not Significant"
      )
    )
  
  pathway_df <- tibble(
    pathway = c(
      "Serotonergic Synapse Maturation", "Axon Guidance & Arborization",
      "Chromatin Remodeling Complex", "Synaptic Vesicle Exocytosis",
      "Neurotransmitter Transport", "Oxidative Phosphorylation",
      "RNA Splicing & Processing", "Cell Cycle Regulation"
    ),
    NES = c(2.45, 2.21, 1.98, 1.87, 1.76, -1.65, -1.82, -2.15),
    padj = c(1.2e-6, 4.5e-5, 3.1e-4, 1.1e-3, 4.2e-3, 8.5e-3, 2.1e-3, 3.4e-4),
    gene_count = c(48, 64, 39, 52, 31, 74, 61, 45)
  )
  
  list(de = de_df, pathways = pathway_df)
}

omics_data <- load_omics_data()

# ------------------------------------------------------------------------------
# 2. User Interface
# ------------------------------------------------------------------------------
ui <- page_navbar(
  title = "Multi-Omics & Study Design Workbench",
  theme = bs_theme(version = 5, bootswatch = "flatly", primary = "#1B4965"),
  
  # TAB 1: Integrative RNA-seq & ATAC-seq Explorer
  nav_panel(
    title = "1. RNA + ATAC Explorer",
    layout_sidebar(
      sidebar = sidebar(
        width = 300,
        sliderInput("fc_thresh", "Abs(Log2 Fold Change) Cutoff:", min = 0.5, max = 3, value = 1.0, step = 0.25),
        selectInput("padj_thresh", "Adjusted P-Value Cutoff:", choices = c("0.05" = 0.05, "0.01" = 0.01, "0.001" = 0.001), selected = 0.05),
        selectInput("motif_filter", "Filter by TF Footprint Motif:", choices = c("All", unique(omics_data$de$tf_motif))),
        hr(),
        downloadButton("download_de", "Download Filtered CSV", class = "btn-primary w-100")
      ),
      layout_columns(
        col_widths = c(6, 6),
        card(
          card_header("Transcriptomic Differential Expression (RNA-seq)"),
          plotOutput("volcano_plot", height = "340px")
        ),
        card(
          card_header("Chromatin Accessibility vs. Expression Concordance"),
          plotOutput("concordance_plot", height = "340px")
        )
      ),
      card(
        card_header("Interactive Gene-Level Summary Table"),
        DTOutput("de_table")
      )
    )
  ),
  
  # TAB 2: Pathway & Enrichment Explorer
  nav_panel(
    title = "2. Pathway Enrichment (GSEA)",
    layout_sidebar(
      sidebar = sidebar(
        width = 300,
        sliderInput("nes_cutoff", "Minimum |NES| Score:", min = 1.0, max = 2.5, value = 1.5, step = 0.1),
        helpText("Normalized Enrichment Scores (NES) derived from fgsea / clusterProfiler workflows.")
      ),
      card(
        card_header("Enriched Biological Processes & Regulatory Programs"),
        plotOutput("pathway_plot", height = "420px")
      )
    )
  ),
  
  # TAB 3: Upstream Study Design & Statistical Power Simulator
  nav_panel(
    title = "3. Study Design & Power Consult",
    layout_sidebar(
      sidebar = sidebar(
        width = 320,
        h5("Experimental Parameters"),
        sliderInput("bio_cv", "Biological Coefficient of Variation (CV):", min = 0.15, max = 0.65, value = 0.30, step = 0.05),
        sliderInput("target_fc", "Target Fold Change to Detect:", min = 1.25, max = 3.0, value = 1.5, step = 0.25),
        sliderInput("seq_depth", "Sequencing Depth (Million Reads/Sample):", min = 10, max = 60, value = 30, step = 5),
        numericInput("alpha_level", "False Discovery / Alpha Level:", value = 0.05, min = 0.001, max = 0.1, step = 0.01),
        helpText("Use during pre-sequencing consultations to balance biological replicate count against sequencing depth.")
      ),
      card(
        card_header("Statistical Power vs. Biological Replicates per Group"),
        plotOutput("power_plot", height = "420px"),
        card_footer("Recommendation: Adding biological replicates increases statistical power substantially more than increasing sequencing depth beyond 25–30M reads.")
      )
    )
  )
)

# ------------------------------------------------------------------------------
# 3. Server Logic
# ------------------------------------------------------------------------------
server <- function(input, output, session) {
  
  filtered_de <- reactive({
    df <- omics_data$de %>%
      mutate(
        sig_status = case_when(
          rna_padj < as.numeric(input$padj_thresh) & rna_log2FC >= input$fc_thresh  ~ "Up-regulated",
          rna_padj < as.numeric(input$padj_thresh) & rna_log2FC <= -input$fc_thresh ~ "Down-regulated",
          TRUE ~ "Not Significant"
        )
      )
    if (input$motif_filter != "All") {
      df <- df %>% filter(tf_motif == input$motif_filter)
    }
    df
  })
  
  output$volcano_plot <- renderPlot({
    ggplot(filtered_de(), aes(x = rna_log2FC, y = -log10(rna_padj), color = sig_status)) +
      geom_point(alpha = 0.65, size = 1.8) +
      scale_color_manual(values = c("Up-regulated" = "#D95F02", "Down-regulated" = "#1B9E77", "Not Significant" = "grey75")) +
      geom_vline(xintercept = c(-input$fc_thresh, input$fc_thresh), linetype = "dashed", color = "grey40") +
      geom_hline(yintercept = -log10(as.numeric(input$padj_thresh)), linetype = "dashed", color = "grey40") +
      labs(x = "Log2 Fold Change (RNA-seq)", y = "-Log10 Adjusted P-Value", color = NULL) +
      theme_minimal(base_size = 13) +
      theme(legend.position = "top")
  })
  
  output$concordance_plot <- renderPlot({
    ggplot(filtered_de(), aes(x = atac_log2FC, y = rna_log2FC, color = sig_status)) +
      geom_point(alpha = 0.65, size = 1.8) +
      geom_smooth(method = "lm", se = FALSE, color = "#1B4965", linetype = "dotted") +
      scale_color_manual(values = c("Up-regulated" = "#D95F02", "Down-regulated" = "#1B9E77", "Not Significant" = "grey75")) +
      labs(x = "Promoter/Enhancer Accessibility Log2FC (ATAC-seq)", y = "Transcript Log2FC (RNA-seq)", color = NULL) +
      theme_minimal(base_size = 13) +
      theme(legend.position = "top")
  })
  
  output$de_table <- renderDT({
    filtered_de() %>%
      select(Gene = gene, `RNA Log2FC` = rna_log2FC, `RNA Padj` = rna_padj,
             `ATAC Log2FC` = atac_log2FC, `TF Footprint` = tf_motif, Status = sig_status) %>%
      mutate(across(where(is.numeric), ~ signif(.x, 3))) %>%
      datatable(options = list(pageLength = 6, scrollX = TRUE), rownames = FALSE)
  })
  
  output$download_de <- downloadHandler(
    filename = function() paste0("filtered_multiomics_", Sys.Date(), ".csv"),
    content = function(file) write.csv(filtered_de(), file, row.names = FALSE)
  )
  
  output$pathway_plot <- renderPlot({
    omics_data$pathways %>%
      filter(abs(NES) >= input$nes_cutoff) %>%
      ggplot(aes(x = NES, y = reorder(pathway, NES), size = gene_count, color = -log10(padj))) +
      geom_point() +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      scale_color_viridis_c(option = "plasma", end = 0.85) +
      labs(x = "Normalized Enrichment Score (NES)", y = NULL, size = "Leading Edge Genes", color = "-Log10 Padj") +
      theme_minimal(base_size = 14)
  })
  
  output$power_plot <- renderPlot({
    reps <- 2:12
    # Approximate total dispersion combining biological CV and Poisson shot noise from sequencing depth
    tech_var <- 1 / (input$seq_depth * 0.8)
    total_sd <- sqrt(log(1 + input$bio_cv^2 + tech_var))
    delta <- log(input$target_fc)
    
    power_vals <- sapply(reps, function(n) {
      power.t.test(n = n, delta = delta, sd = total_sd, sig.level = input$alpha_level, type = "two.sample")$power
    })
    
    tibble(Replicates = reps, Power = power_vals) %>%
      ggplot(aes(x = Replicates, y = Power)) +
      geom_line(color = "#1B4965", linewidth = 1.3) +
      geom_point(color = "#1B4965", size = 3) +
      geom_hline(yintercept = 0.80, linetype = "dashed", color = "#D95F02", linewidth = 1) +
      annotate("text", x = 10.5, y = 0.83, label = "80% Target Power", color = "#D95F02", fontface = "bold") +
      scale_x_continuous(breaks = reps) +
      scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
      labs(x = "Biological Replicates per Condition (n)", y = "Statistical Power (1 - Beta)") +
      theme_minimal(base_size = 14)
  })
}

shinyApp(ui, server)
