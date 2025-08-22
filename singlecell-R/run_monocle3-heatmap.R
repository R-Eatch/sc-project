############################################
##  0)  Load packages
############################################
pkgs <- c("Seurat", "monocle3", "ClusterGVis",
          "ComplexHeatmap", "org.Mm.eg.db",
          "dplyr", "ggplot2")
stopifnot(all(sapply(pkgs, require, character.only = TRUE)))

############################################
##  1)  Global settings
############################################
dataset          <- "R-MG"
cds_file         <- paste0("D:/Data/", dataset, "_subset_cds.rds")
target_types     <- c(#"Basal",
                      "LumHR",
                      #"LumSEC-Lac",
                      #"LumSEC-Lip",
                      #"LumSEC-Mgp",
                      #"LumSEC-Vcam1",
                      "MaSC",
                      "MaSC-Pro"
)
cluster_method   <- "kmeans"
cluster_num      <- 3
go_topn          <- 8
kegg_topn        <- 8
orgdb            <- org.Mm.eg.db
kegg_code        <- "mmu"
cores_graph_test <- 4
plot_w           <- 10
plot_h           <- 10
outdir           <- "cluster_results"; dir.create(outdir, showWarnings = FALSE)
markGenes        <- c('Esr1','Prlr','Pgr','Kit','Stc2','Lef1','Krt5')
lineage          <-"HR"
############################################
##  2)  Load CDS  &  subset cell types
############################################
cds <- load_monocle_objects(directory_path = cds_file)
plot_cells(cds,
                 color_cells_by = "newcelltype",
                 label_cell_groups=TRUE,
                 label_leaves=TRUE,
                 label_branch_points=TRUE,
                 graph_label_size=3)
cds_sub <- cds[, cds@colData@listData[["newcelltype"]] %in% target_types]

############################################
##  3)  graph_test  →  CSV
############################################
gt <- graph_test(cds_sub, neighbor_graph = "principal_graph",
                 cores = cores_graph_test)
write.csv(gt,
          file.path(outdir, paste0(dataset, "_graph_test.csv")),
          row.names = TRUE)

genes_use <- rownames(gt)[gt$q_value == 0 & gt$morans_I > 0.25]

############################################
##  4)  Pseudotime matrix  →  clusterData
############################################
mat <- pre_pseudotime_matrix(cds_obj  = cds_sub,
                             gene_list = genes_use) |>
  as.matrix()

ck <- clusterData(obj            = mat,
                  cluster.method = cluster_method,
                  cluster.num    = cluster_num)

## export gene ↔ cluster
write.csv(ck$wide.res[, c("gene","cluster")],
          file.path(outdir, paste0(dataset, "_gene_cluster.csv")),
          row.names = FALSE)

############################################
##  5)  Enrichment (GO / KEGG)
############################################
ego <- enrichCluster(ck, type = "BP",
                     OrgDb = orgdb,
                     id.trans = TRUE,
                     fromType = "SYMBOL", toType = "ENTREZID",
                     pvalueCutoff = 0.05, topn = go_topn)

ekegg <- enrichCluster(ck, type = "KEGG",
                       OrgDb = orgdb, organism = kegg_code,
                       id.trans = TRUE,
                       fromType = "SYMBOL", toType = "ENTREZID",
                       pvalueCutoff = 0.05, topn = kegg_topn)

#termanno <- ego   |> transmute(id   = group,
#                               term = Description,
#                               pval = pvalue,
#                               ratio = ratio)
#kegganno <- ekegg |> transmute(id   = as.integer(group),
#                               term = Description,
#                               pval = pvalue,
#                               ratio = ratio)

############################################
##  6)  Helper to save plots
############################################
save_ht <- function(ht, stub){
  pdf (file.path(outdir, paste0(stub, ".pdf")),
       width = plot_w, height = plot_h); ComplexHeatmap::draw(ht); dev.off()
  png (file.path(outdir, paste0(stub, ".png")),
       width = 2000, height = 2400, res = 300); ComplexHeatmap::draw(ht); dev.off()
}

############################################
##  7-1)  PLOT 1 : trend-only (trend on right)
############################################

p_trend <- visCluster(ck,
                      plot.type     = "both",
                      #line.side     = "right",   # trend ► right
                      #show_row_dend = FALSE,
                      markGenes     = markGenes,
                      add.sampleanno = F)      # no GO/KEGG

save_ht(p_trend, paste0(dataset, "_TrendOnly"))

############################################
##  7-2)  PLOT 2 : GO  (trend left, GO right)
############################################
p_go <- visCluster(ck,
                   plot.type       = "both",
                   line.side = "left",
                   show_row_dend   = FALSE,
                   annoTerm.data   = ego,
                   add.sampleanno = F,
                   markGenes     = markGenes,
                   markGenes.side = "left",
                   go.col = rep(ggsci::pal_d3()(cluster_num),each = go_topn,go.size = "pvalue")
                   )

save_ht(p_go, paste0(dataset, lineage,"_GO_top", go_topn))

############################################
##  7-3)  PLOT 3 : KEGG (trend left, KEGG right)
############################################
p_kegg <- visCluster(ck,
                     plot.type       = "both",
                     line.side = "left",
                     show_row_dend   = FALSE,
                     annoTerm.data   = ekegg,
                     add.sampleanno = F,
                     markGenes     = markGenes,
                     markGenes.side = "left",
                     go.col = rep(ggsci::pal_d3()(cluster_num),each = go_topn,go.size = "pval")
)

save_ht(p_kegg, paste0(dataset,lineage, "_KEGG_top", kegg_topn))

message("✓ All done for ", dataset,
        ".  Files saved in: ", normalizePath(outdir))

