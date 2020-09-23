# Define UI for app that draws a histogram ----
setwd("/hpcnfs/scratch/DP/frossi/ANALYSIS/frossi/Rmd_Hackathon2")

ui <- fluidPage(
  tabsetPanel(
  tabPanel("VolcanoPlot",
           sidebarLayout(
             sidebarPanel(helpText("Select a log2 fold-change value"),
                          sliderInput(inputId = "l2fc", label = "l2fc:", min = 1, max = 5, value = 2, step = 0.5),
                          hr(),
                          helpText("Select an adjusted p-value"),
                          sliderInput(inputId = "adj_pvalue", label = "adj_pvalue:", min = 0.01, max = 0.1, value = 0.05, step = 0.01)),
                 mainPanel(plotOutput(outputId = "VolcanoPlot")))),
  
  tabPanel("UpSetPlot",
           sidebarLayout(
             sidebarPanel(helpText("Select a log2 fold-change value"),
                          sliderInput(inputId = "l2fcUpSet", label = "l2fc:", min = 1, max = 5, value = 2,
                                      step = 0.5),
                          hr(),
                          helpText("Select an adjusted p-value"),
                          sliderInput(inputId = "adj_pvalueUpSet", label = "adj_pvalue:", min = 0.01,
                                      max = 0.1, value = 0.05, step = 0.01),
                          hr(),
                          helpText("Select a protein complex to intersect with \n differentially expressed genes"),
                          radioButtons("radio", label = "Gene intersection", choices = c("ComplexA" = "TFA", "ComplexB" = "TFB"))
             ),
             mainPanel(plotOutput(outputId = "UpSetPlot")))),
  
  tabPanel("ChromosomePlot",
           sidebarLayout(
             sidebarPanel(helpText("Select a chromosome"),
                          selectInput("chromosome", "Chromosome:", choices = unique(tab[, 1])),
                          hr(),
                          helpText("Select a sample"),
                          selectInput("sample", "Sample:", choices = c("WT", "KO_rep1", "KO_rep2")),
                          hr(),
                          helpText("Select an histone modification"),
                          selectInput("histMod", "Histone modification:", choices = c("H3K4me3", "H3K27ac"))
           ),
           mainPanel(plotOutput(outputId = "ChrPlot"))))
  
))

   


server <- function(input, output) {
  
  output$VolcanoPlot <- renderPlot({
    
    setwd("/hpcnfs/scratch/DP/frossi/ANALYSIS/frossi/Rmd_Hackathon2")
    
    l2fc <- input$l2fc
    adj_pvalue <- input$adj_pvalue
    
    df <- read.table("DEG.tab2", header = T)
    
    df <- df[complete.cases(df), ]
    
    df <- df %>%
      dplyr::select(log2FoldChange, padj) %>%
      dplyr::mutate(DEG = case_when(padj < adj_pvalue & log2FoldChange < -l2fc ~ 'Downregulated',
                                    padj < adj_pvalue & log2FoldChange > l2fc ~ 'Upregulated',
                                    log2FoldChange <= l2fc & log2FoldChange >= -l2fc ~ 'NotDE',
                                    padj > adj_pvalue ~ 'NotDE'))
    
    ggplot(data = df, aes(x = log2FoldChange, y = -log10(padj), colour = DEG)) +
      geom_point(size = 1.5) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5)) +
      coord_cartesian(ylim = c(0, max(-log10(df$padj))+8)) +
      scale_colour_manual(values = c("darkgreen", "gray", "red")) +
      annotate("text", x = 3.5, y = max(-log10(df$padj))+8, label = paste("log2FoldChange =", l2fc)) +
      annotate("text", x = 5, y = max(-log10(df$padj))+4, label = nrow(df[df$DEG=='Upregulated',]), colour = "red") +
      annotate("text", x = -5, y = max(-log10(df$padj))+4, label = nrow(df[df$DEG=='Downregulated',]), colour =
                 "darkgreen") +
      annotate("text", x = -3.5, y = max(-log10(df$padj))+8, label = paste("adjusted p-value", adj_pvalue)) +
      labs(x = "log2FoldChange",
           y = "-log10 adj p-value",
           title = "Treated-vs-Ctrl")
    
  })
  
  output$UpSetPlot <- renderPlot({
    
    l2fc = input$l2fcUpSet
    adj_pvalue <- input$adj_pvalueUpSet
    
    df <- read.table("DEG.tab2", header = T)
    
    df <- df[complete.cases(df), ]
    
    df <- df %>%
      dplyr::select(Geneid, log2FoldChange, padj) %>%
      dplyr::mutate(DEG = case_when(padj < adj_pvalue & log2FoldChange < - l2fc ~ 'Downregulated',
                                    padj < adj_pvalue & log2FoldChange > l2fc ~ 'Upregulated',
                                    log2FoldChange <= l2fc & log2FoldChange >= -l2fc ~ 'NotDE',
                                    padj > adj_pvalue ~ 'NotDE'))
    
    UP <- df %>%
      dplyr::filter(DEG=="Upregulated") %>%
      dplyr::select(Geneid)
    
    DOWN <- df %>%
      dplyr::filter(DEG=="Downregulated") %>%
      dplyr::select(Geneid)
    
    TFA_1 <- read.table("TFA_1")
    TFA_1$id1 <- "TFA_1"
    TFA_1$id2 <- "TFA"
    
    TFA_2 <- read.table("TFA_2")
    TFA_2$id1 <- "TFA_2"
    TFA_2$id2 <- "TFA"
    
    TFB_1 <- read.table("TFB_1")
    TFB_1$id1 <- "TFB_1"
    TFB_1$id2 <- "TFB"
    
    TFB_2 <- read.table("TFB_2")
    TFB_2$id1 <- "TFB_2"
    TFB_2$id2 <- "TFB"
    
    tab <- rbind(TFA_1, TFA_2, TFB_1, TFB_2)
    
    tab_TF <- tab[tab$id2==input$radio,] %>%
      dplyr::select(id1) %>%
      unique()
    
    TF1_ids <- tab[tab$id1==tab_TF[1,], 1]
    TF2_ids <- tab[tab$id1==tab_TF[2,], 1]
    
    combined <- reduce(list(data.frame(gene = as.character(UP$Geneid), UP = 1),
                            data.frame(gene = as.character(DOWN$Geneid), DOWN = 1),
                            data.frame(gene = as.character(TF1_ids), TF1_ids = 1),
                            data.frame(gene = as.character(TF2_ids), TF2_ids = 1)
    ), full_join)
    
    combined[is.na(combined)] <- 0
    
    upset(combined,
          group.by = "degree",
          order.by = "degree",
          nsets = 4,
          #empty.intersections = "on",
          matrix.color = "blue",
          sets.bar.color=rev(c("red", "darkgreen", "orange", "blue")))
    
})
  
  output$ChrPlot <- renderPlot({
    
    H3K4me3 <- read.table("/hpcnfs/scratch/DP/frossi/ANALYSIS/econway/ChIPseq/mouse/06bigwig/noSubtract/histone_mod/tabA")
    tabA$id <- "H3K4me3"
    
    H3K27ac <- read.table("/hpcnfs/scratch/DP/frossi/ANALYSIS/econway/ChIPseq/mouse/06bigwig/noSubtract/histone_mod/tabB")
    tabB$id <- "H3K27ac"
    
    tab <- rbind(tabA, tabB)
    
    colnames(tab) <- c("chr", "start", "end", "WT",
                       "KO_rep1", "KO_rep2", "id")
    
    tabChr <- tab %>%
      dplyr::filter(chr == input$chromosome & id == input$histMod) %>%
      dplyr::select(WT, KO_rep1, KO_rep2) %>%
      reshape2::melt() %>%
      dplyr::filter(variable == input$sample)
    
    
    sampleMean <- rollapply(tabChr[, 2], width = 5, by = 5, FUN = mean)
    
    df = data.frame(pos = 1:length(sampleMean), mean = sampleMean)
    
    df$id = input$histMod
    
    ggplot(data = df, aes(x = pos, y = mean, colour = id)) +
      geom_line(size = 1) +
      coord_cartesian(ylim = c(0, max(df$mean)+0.5)) +
      theme_classic(base_size=13) +
      labs(x = "50 kbp window",
           y = "mean (CPM)",
           title = paste(input$chromosome, input$histMod, input$sample)) +
      theme(legend.position = "none",
            axis.text = element_text(size = 16),
            axis.title = element_text(size = 18),
            legend.title = element_blank(),
            plot.title = element_text(size = 22, hjust = 0.5),
            legend.text = element_text(size = 18)) +
      scale_colour_manual(values = "blue")
  })
  
}

shinyApp(ui = ui, server = server)
