library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(grid)

# Load Seurat object
ea <- readRDS("D:/111/MRS-MGCCA_combined.rds")

# Extract and format gene list from file
gene_list <- c("Prrx1","Cntn3", "Rarb", "Adamtsl4", "Hspg2", "Slc9a9", "Clu", "Plod2", "Ctnnal1",
               "Irak2","Lgr5", "Cap2", "Fbln2", "Slc28a3", "Plod2", "Rarb", "Fabp5", "Tnfrsf19", "Chka", "Fam20a", "Lmo7")#正选择
gene_list <- c(
  "Dpt", "Slc24a5", "Chrm3", "Prkcg", "Dlk1", "Lama2", "Ogn", "Limch1", "Lsamp", "Thsd4", 
  "Asb9", "Cacna1b", "Edil3", "Ctss", "Irf4", "Pde9a", "Epgn", "St8sia4", "Cd274", "Epsti1", 
  "Klk9", "Ccdc3", "Atp6v1g3", "Ccdc190", "Slc35f3", "Wif1", "Hmgcs2", "Fabp9", "Fsip1", "Plppr1"
)#MG中的DEG


gene_list <- c("Prrx1","Cntn3" ,"Adamtsl4", "Hspg2","Ctnnal1","Lgr5", "Cap2", "Fam20a")#正选择且表达模式存在差异的

gene_list <- c("Clu", "Slc28a3","Irak2","Plod2","Rarb" , "Fabp5", "Tnfrsf19")#正选择且表达模式一致组

gene_list <- c(
  "Mylk", "Ctsc", "Trps1", "Hspb1", "Kalrn", "Ptprk", "Matn2", "Slit2", "Nfatc2", "Moxd1",
  "Col3a1", "Igfbp2", "Cdh13", "Prune2", "Slc35f1", "Cacnb2", "Cd14", "Rhoj", "Synpo2", "Krt17",
  "Hp", "Lama3", "Smtn", "Tnc", "Gask1b", "Htra1", "Cd44", "Zeb2", "Ism1", "Col6a2",
  "Aff3", "Csrp1", "Sox5", "Robo2", "Sorbs1", "Nrcam", "Actg2", "S100a10", "Mia", "Filip1l",
  "Mfge8", "Lalba", "Cblb", "Gfpt2", "Msrb3", "Mpped2", "Epha4", "Atp6v1c2", "Arhgap20", "Lpgat1",
  "Dsp", "Sparc", "Fabp5", "Gli3", "Atp2b2", "Thsd7a", "Galnt18", "Fhod3", "Cav1", "Slc28a3",
  "Sema3a", "Antxr1", "Igfbp6", "Lpl", "Cd74", "Ltf", "Dach1", "Lcp1", "Rbms3", "Lmo7",
  "Ptn", "Tnfrsf19", "Lgals3", "Homer1", "Dpyd", "Bicc1", "Lcn2", "Lama4", "Btn1a1", "Mpz",
  "Ccn2", "Pigr"
)



# Define time and sample columns
time_column <- "time"
sample_column <- "merge_dataset"

# Extract expression data and metadata
expression_data <- GetAssayData(ea, slot = "data")
metadata <- ea@meta.data

# Get all unique time points and samples
all_times <- sort(unique(metadata[, time_column]))
all_samples <- unique(metadata[, sample_column])

# Function: Plot gene expression over time for a single gene
plot_gene_expression <- function(gene_name) {
  # Print the gene name
  print(paste("Plotting gene:", gene_name))
  
  # Merge expression data and metadata
  data <- data.frame(expression = expression_data[gene_name, ], 
                     time = metadata[, time_column],
                     sample = metadata[, sample_column])
  
  # Summarize mean expression by sample and time
  mean_expression <- data %>%
    group_by(sample, time) %>%
    summarize(mean_expression = mean(expression))
  
  # Plot line chart with legend
  p <- ggplot(mean_expression, aes(x = factor(time), y = mean_expression, color = sample, group = sample)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    ggtitle(gene_name) +
    xlab("Time") +
    ylab("Average Expression") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +  # Centered title
    scale_color_brewer(palette = "Set1") +
    scale_x_discrete(breaks = as.character(all_times))
  
  return(p)
}

# Plot expression charts for all genes
plots <- lapply(gene_list, plot_gene_expression)

# Arrange all plots in a grid with four plots per row
num_cols <- 4
num_rows <- ceiling(length(plots) / num_cols)

# Save the combined plot as a high-resolution file
output_file <- "MRS_MG_line.png"

# Create a canvas for saving the image with high resolution
png(output_file, width = 4000, height = num_rows * 1000, res = 300)

# Draw all plots and add a main title
grid.arrange(grobs = plots, ncol = num_cols, 
             top = textGrob("shared gene", gp = gpar(fontsize = 20, fontface = "bold")))

# Close the device
dev.off()
