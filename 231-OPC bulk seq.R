# Load necessary libraries
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)  # Mouse gene annotation
library(enrichplot)
library(DOSE)

# File path to RNA-seq data
file_path <- "E:/seq/OPC/OPC seq lipao.csv"

# Load the data
rna_data <- read.csv(file_path)

# Filter DEGs based on p-value threshold
deg <- rna_data[rna_data$pvalue < 0.02, ]

# Generate the volcano plot
deg$color <- ifelse(deg$log2FoldChange > 0, "#A50F15", "#08519C")
ggplot(deg, aes(x = log2FoldChange, y = -log10(pvalue), color = color)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_identity() +
  theme_minimal() +
  labs(title = "Volcano Plot of DEGs",
       x = "Log2 Fold Change",
       y = "-Log10 P-value") +
  theme(plot.title = element_text(hjust = 0.5))

# Prepare gene list for enrichment analysis
# Assumes 'gene_id' column has gene symbols or Entrez IDs
gene_list <- deg$gene_id

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
ggsave("E:/seq/OPC/volcano_plot.png", width = 8, height = 6, dpi = 300)

# Save the GO and KEGG enrichment results
write.csv(as.data.frame(go_enrich), "E:/seq/OPC/GO_enrichment_results.csv", row.names = FALSE)
write.csv(as.data.frame(kegg_enrich), "E:/seq/OPC/KEGG_enrichment_results.csv", row.names = FALSE)

# Filter GO enrichment results for specific terms
specific_go_terms <- c(
  "gliogenesis",
  "positive regulation of cell adhesion",
  "negative regulation of nervous system development",
  "negative regulation of neurogenesis",
  "glial cell differentiation",
  "negative regulation of cell development",
  "positive regulation of cytokine production",
  "synapse organization",
  "positive regulation of response to external stimulus",
  "negative regulation of cellular component movement",
  "neuron to neuron synapse",
  "receptor complex",
  "membrane microdomain",
  "postsynaptic membrane",
  "collagen-containing extracellular matrix",
  "organelle subcompartment",
  "transport vesicle",
  "exocytic vesicle",
  "phagocytic vesicle",
  "early endosome",
  "Ras GTPase binding",
  "transmembrane receptor protein kinase activity",
  "phospholipid binding",
  "integrin binding",
  "actin filament binding",
  "metal ion transmembrane transporter activity",
  "anion transmembrane transporter activity",
  "passive transmembrane transporter activity",
  "acidic amino acid transmembrane transporter activity",
  "lipid transporter activity"
)

# Subset GO enrichment results to only the terms of interest
filtered_go <- go_enrich@result[go_enrich@result$Description %in% specific_go_terms, ]

# Check if any terms matched
if (nrow(filtered_go) == 0) {
  stop("No GO terms matched your specified terms.")
}

# Convert filtered GO results to enrichResult object for plotting
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

# Plot specific GO terms
dotplot(filtered_go_result, showCategory = length(specific_go_terms)) +
  ggtitle("Specific GO Terms Enrichment") +
  theme(plot.title = element_text(hjust = 0.5))

# Filter KEGG enrichment results for specific pathways
specific_kegg_terms <- c(
  "gliogenesis",
  "positive regulation of cell adhesion",
  "negative regulation of nervous system development",
  "negative regulation of neurogenesis",
  "glial cell differentiation",
  "negative regulation of cell development"
)  # Add your KEGG pathway terms here

filtered_kegg <- kegg_enrich@result[kegg_enrich@result$Description %in% specific_kegg_terms, ]

# Check if any pathways matched
if (nrow(filtered_kegg) == 0) {
  stop("No KEGG pathways matched your specified terms.")
}

# Convert filtered KEGG results to enrichResult object for plotting
filtered_kegg_result <- new("enrichResult",
                            result = filtered_kegg,
                            gene = kegg_enrich@gene,
                            geneSets = kegg_enrich@geneSets,
                            pvalueCutoff = kegg_enrich@pvalueCutoff,
                            pAdjustMethod = kegg_enrich@pAdjustMethod,
                            qvalueCutoff = kegg_enrich@qvalueCutoff,
                            organism = kegg_enrich@organism,
                            keytype = kegg_enrich@keytype,
                            readable = kegg_enrich@readable
)

# Plot specific KEGG pathways
dotplot(filtered_kegg_result, showCategory = length(specific_kegg_terms)) +
  ggtitle("Specific KEGG Pathways Enrichment") +
  theme(plot.title = element_text(hjust = 0.5))

# Save plots as images
ggsave("E:/seq/OPC/specific_go_dotplot.png", width = 8, height = 6, dpi = 300)
ggsave("E:/seq/OPC/specific_kegg_dotplot.png", width = 8, height = 6, dpi = 300)

# Save filtered results to CSV files
write.csv(filtered_go, "E:/seq/OPC/filtered_GO_terms.csv", row.names = FALSE)
write.csv(filtered_kegg, "E:/seq/OPC/filtered_KEGG_pathways.csv", row.names = FALSE)
