# Load necessary libraries
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)  # Mouse gene annotation
library(enrichplot)
library(DOSE)

file_path <- "E:/seq/Neuron/Neuron seq lipao.csv"

# Load the data
rna_data <- read.csv(file_path)

# Filter DEGs based on p-value threshold (p < 0.05)
deg <- rna_data[rna_data$diffexp_deseq2_pvalue_Control.vs.cKO < 0.05, ]

# Generate the volcano plot
deg$color <- ifelse(deg$diffexp_log2fc_Control.vs.cKO > 0, "#A50F15", "#08519C")
ggplot(deg, aes(x = diffexp_log2fc_Control.vs.cKO, y = -log10(diffexp_deseq2_pvalue_Control.vs.cKO), color = color)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_identity() +
  theme_minimal() +
  labs(title = "Volcano Plot of DEGs",
       x = "Log2 Fold Change",
       y = "-Log10 P-value") +
  theme(plot.title = element_text(hjust = 0.5))

# Prepare gene list for enrichment analysis
# Assumes 'gene_symbol' column has gene symbols
gene_list <- deg$gene_symbol

# Perform GO enrichment analysis
go_enrich <- enrichGO(
  gene = gene_list,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",  # Change to "ENTREZID" if IDs are Entrez
  ont = "ALL",
  pvalueCutoff = 0.05,
  readable = TRUE
)

# Visualize GO enrichment results
dotplot(go_enrich, showCategory = 20) +
  ggtitle("GO Enrichment Analysis")

# Perform KEGG pathway enrichment analysis
# Map gene symbols to Entrez IDs if needed
gene_entrez <- bitr(
  gene_list,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

kegg_enrich <- enrichKEGG(
  gene = gene_entrez$ENTREZID,
  organism = "mmu",  # Mouse organism code
  pvalueCutoff = 0.05
)

# Visualize KEGG enrichment results
dotplot(kegg_enrich, showCategory = 20) +
  ggtitle("KEGG Pathway Enrichment Analysis")

# Save the volcano plot to a file
ggsave("E:/seq/volcano_plot_new_data.png", width = 8, height = 6, dpi = 300)

# Save the GO and KEGG enrichment results
write.csv(as.data.frame(go_enrich), "E:/seq/GO_enrichment_results_new_data.csv", row.names = FALSE)
write.csv(as.data.frame(kegg_enrich), "E:/seq/KEGG_enrichment_results_new_data.csv", row.names = FALSE)

# Filter GO enrichment results for new specific terms
specific_go_terms <- c(
  "regulation of exocytosis",
  "lipid localization",
  "lipid transport",
  "regulation of lipid metabolic process",
  "neuron death",
  "cellular senescence",
  "mitotic nuclear division",
  "neural precursor cell proliferation",
  "developmental cell growth",
  "lysosome",
  "vesicle membrane",
  "exocytic vesicle",
  "lysosomal membrane",
  "endocytic vesicle membrane",
  "amyloid-beta binding",
  "Cysteine-type endopeptidase activity involved in apoptotic process",
  "Phospholipid binding",
  "GTPase activity",
  "Cellular senescence"
)

filtered_go <- go_enrich@result[go_enrich@result$Description %in% specific_go_terms, ]

# Plot specific GO terms
if (nrow(filtered_go) > 0) {
  filtered_go_result <- new("enrichResult",
                            result = filtered_go,
                            gene = go_enrich@gene,
                            geneSets = go_enrich@geneSets,
                            pvalueCutoff = go_enrich@pvalueCutoff,
                            pAdjustMethod = go_enrich@pAdjustMethod,
                            qvalueCutoff = go_enrich@qvalueCutoff,
                            organism = go_enrich@organism,
                            keytype = go_enrich@keytype,
                            readable = go_enrich@readable,
                            ontology = go_enrich@ontology
  )
  
  dotplot(filtered_go_result, showCategory = length(specific_go_terms)) +
    ggtitle("Specific GO Terms Enrichment (Filtered)") +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave("E:/seq/specific_go_dotplot_new_data.png", width = 8, height = 6, dpi = 300)
}

# Save filtered GO results
write.csv(filtered_go, "E:/seq/filtered_GO_terms_new_data.csv", row.names = FALSE)
