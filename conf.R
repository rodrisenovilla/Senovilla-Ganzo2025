## PhyloBrain cross-species config (public version)
## Configure paths via environment variables if needed:
##   PB_BASE_DIR, PB_ORTHO_DIR, PB_OBJECTS_DIR, PB_MARKERS_FILE, PB_GECKO_MAP_FILE, PB_OUT_DIR

suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(SeuratObject)
  library(patchwork)
  library(Matrix)
  library(glmGamPoi)
  library(sctransform)
  library(Hmisc)
  library(ComplexHeatmap)
  library(stringr)
  library(ggplot2)
  library(biomaRt)
  library(RColorBrewer)
  library(circlize)   # colorRamp2 lives here
  library(corrplot)
  library(pheatmap)
  library(RANN)
})

## ---------- Paths (portable) ----------
PB_BASE_DIR <- Sys.getenv("PB_BASE_DIR", ".")
PB_PATHS <- list(
  orthologues = Sys.getenv("PB_ORTHO_DIR",  file.path(PB_BASE_DIR, "orthologues")),
  objects     = Sys.getenv("PB_OBJECTS_DIR", file.path(PB_BASE_DIR, "pb_objects")),
  markers     = Sys.getenv("PB_MARKERS_FILE", file.path(PB_BASE_DIR, "markers_091024.csv")),
  gecko_map   = Sys.getenv("PB_GECKO_MAP_FILE", file.path(PB_BASE_DIR, "gecko", "map_symbol.csv")),
  out         = Sys.getenv("PB_OUT_DIR", PB_BASE_DIR)
)

.ensure_file <- function(path, msg) {
  if (!file.exists(path)) stop(msg, ": ", path, call. = FALSE)
}
.ensure_dir  <- function(path) if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)

## ---------- Species metadata ----------
objects <- list(mouse="me9", human="cs13", chicken="ce3", bgalgal1="ce3", gecko="ge5", zebrafish="zbf161820")
animal  <- list(mouse="steelblue", human="orange2", chicken="firebrick", bgalgal1="firebrick", gecko="gold3", zebrafish="forestgreen")
general_ct <- c("gold4","steelblue4","green4","firebrick4","tan4","black")

## Palettes per species
colors_cs13 <- c("#ffd000","#eac435","#f6eb1e","#03045e","#023e8a","#8fe1ef","#bbd0ff","#1982c4","#1a749e","#0077b6",
                 "#48cae4","#90e0ef","#1a759f","#55a630","#2b9348","#307d19","#124004","#a9d117","#e85d04","#d68706",
                 "#f15025","#9e5f5f","#dc2f02","#f42b03","#df8080","#e35053","#ef233c","#da627d","#942b3b","#ec8385",
                 "#f50202","firebrick","indianred2","coral3","tan4","darkorange3","tan","#c08552","#895737","#5e3023",
                 "#9e9fba","#b29eba","#BDBDBD","#6c757d","#252525")

colors_ce3 <- c("#eeefa8","#dddf00","#f6eb1e","#eac435",
                "#023e8a","#90e0ef","#03045d","#1982c4","#48cae4","#1a749e","#3d8eff","#1a759f","#0077b6","#ade8f4",
                "#656d4a","#2b9348","#55a630","#007f5f",
                "#e85d04","#d68706","#dc2f02",
                "#da627d","#942b3b","#f29479","#f42b03",
                "#c08552","tan","tan4","peru","darkgoldenrod3",
                "#737373","#BDBDBD","#6c757d","grey20","#bbd0ff","#474b4e","#BDBDBD","black")

colors_bgalgal1 <- colors_ce3

colors_me9 <- c("#eeefa8","#f6eb1e","#eac435","#f2eb1e","#03045e","#90e0ef","#bbd0ff","#1982c4","#3d8eff","#1a759f",
                "#0077b6","#48cae4","#55a630","#2b9348","#e85d04","#d68706","#dc2f02","#ead5dc","#f42b03","#df8080",
                "#e35053","#ef233c","#da627d","#fcbba1","#942b3b","#ec8385","#f15025","#c08552","#895737","#5e3023",
                "grey","#737373","#BDBDBD","#e5e5e5","#6c757d","#252525")

colors_zbf161820 <- c("#eac435","#f6eb1e","#dddf00","#3d8eff","#023f89","#90e0ef","#1a759f","#023e8a","#ade8f4","#1982c4",
                      "steelblue","#656d4a","#2b9348","#55a630","#007f5f","#e85d04","#db2f02","#e35053","#f42b03","#ec8385",
                      "#f29479","#f15025","#ef233c","#e2294f","#da627d","#dc2f02","#942b3b","#c08552","#895737","#5e3023",
                      "#ead5dc","#dee2ff","#bbd0ff","#d3d5de","grey75","grey40","#343a40")

colors_ge5 <- c("#ffd000","#eac435","#f6eb1e","#03045e","#023e8a","#90e0ef","#1982c4","#1a759f","#0077b6","#48cae4",
                "#3d8eff","#bbd0ff","steelblue","#55a630","#2b9348","#307d19","#124004","#d68706","#e85d04","#dc2f02",
                "#f42b03","#df8080","#e35053","#da627d","indianred2","#c08552","#895737","#9e9fba","#b29eba","#BDBDBD",
                "#6c757d","#252525","grey","grey30","black")

color <- list(
  mouse = colors_me9,
  chicken = colors_ce3,
  bgalgal1 = colors_bgalgal1,
  human = colors_cs13,
  gecko = colors_ge5,
  zebrafish = colors_zbf161820
)

## ---------- Cell type order (091024) ----------
ce3_ct <- c("AMP","Hem","Die-Mes-RP","Rho-RP","dTel","vTel","OV","POA","THT-PHT","P3","aP2","bP2","P1","PG","aMes","bMes","a'lMes","plMes","a'It","pIt","RL","aR1-2","aR3-8","lR1-8","bR1-8","Die-FP","Mes-FP","It-FP","Rho-dFP","Rho-vFP","DLL1-N","TAL1-N","TAL2-N","LHX9-N","LHX3-N","DLX-N","TBX20-N","GAP43-N")
bgalgal1_ct <- ce3_ct
ge5_ct <- c("AMP","Mes-RP","Rho-RP","OV","dTel","vTel","POA","ATO","THT","PHT","P3","aP2-1","bP2-1","aMes","a'aMes","plMes","bMes","a'It","pIt","RL","aR1-6","l-bR1-2","lR3-6","bR3","bR4-6","Rho-FP","Die-Mes-It-FP","DLL1-N","DLX-N","OLF-GAP43-N","OLF-GADD45G-N","VSX1-N","TAL-N","EGL-N","GAP43-N")
me9_ct <- c("AMP","Hem","Die-Mes-RP","Rho-RP","Pal","SPal","OV","ATO","POA-THT","PHT","P3","P2-1","aMes","bMes","a'It","pIt","RL","aR1","aR2","aR3","aR4","aR5","lR1-2","lR3","lR4","bR0-2","bR3-5","Mes-FP","It-FP","Rho-FP","Phox2-N","Tal2-N","Msx3-N","Lhx-N","EGL-N","Gap43-N")
cs13_ct <- c("AMP","Hem","Mes-RP","Pal","SPal","POA","THT","PHT","PThE","PTh","bPTh","P2","P1","amMes","a'a-lMes","paMes","plMes","bMes","a'It","bIt","pIt","RL","aR1-2","aR3-4","aR5-6","aR7-8","lR1-2","lR3-4","lR5-6","lR7-8","bR1-2","bR3-4","bR5-6","bR7-8","Die-FP","ZLI","Mes-It-dFP","It-vFP","Mes-vFP","Rho-FP","HT-N","SPal-N","MSX-N","NGN1-N","GAP43-N")
zbf161820_ct <- c("AMP","Hem-PG-Mes-RP","Rho-RP","Pal","SPal","OV","ATO","THT","PHT","P3","P2-1","a'aMes","paMes","lMes","bMes","It","aR1-2","aR3","aR5","aR7-8","lR2-4","lR3-4","lR7-8","bR1-4","bR5","bR7-8","R6","Die-FP","Mes-FP","Rho-FP","dla-N","Pal-N","SPal-N","dlb-N","EGL-N","phox2-N","gap43-N")

order <- list(
  mouse = me9_ct,
  chicken = ce3_ct,
  bgalgal1 = bgalgal1_ct,
  human = cs13_ct,
  gecko = ge5_ct,
  zebrafish = zbf161820_ct
)

## ---------- Vesicle / region blocks ----------
ce3_v <-       c(rep("DO",4),rep("FB",10),rep("MB",5),rep("HB",6),rep("FP",5),rep("N",8))
bgalgal1_v <-  ce3_v
me9_v <-       c(rep("DO",4),rep("FB",8),rep("MB",3),rep("HB",12),rep("FP",3),rep("N",6))
zbf161820_v <- c(rep("DO",3),rep("FB",8),rep("MB",5),rep("HB",11),rep("FP",3),rep("N",7))
cs13_v <-      c(rep("DO",3),rep("FB",10),rep("MB",7),rep("HB",14),rep("FP",6),rep("N",5))
ge5_v <-       c(rep("DO",3),rep("FB",10),rep("MB",5),rep("HB",7),rep("FP",2),rep("N",8))

vesicle <- list(
  mouse = me9_v,
  chicken = ce3_v,
  bgalgal1 = bgalgal1_v,
  human = cs13_v,
  gecko = ge5_v,
  zebrafish = zbf161820_v
)

## ---------- Regression variables per species ----------
vars.out <- list(
  mouse = c("percent.mt","G2M.Score","S.Score"),
  chicken = c("S.Score","G2M.Score","percent.mt"),
  bgalgal1 = c("S.Score","G2M.Score","percent.mt"),
  human = c("S.Score","G2M.Score","percent.mt"),
  gecko = c("percent.mt","percent.rpl","S.Score","G2M.Score"),
  zebrafish = c("S.Score","G2M.Score","percent.mt")
)

## ---------- File prefixes / species tags ----------
obj <- list(mouse="me9", chicken="ce3", bgalgal1="bgalgal1", human="cs13", gecko="ge5", zebrafish="zbf161820")
siglas <- list(mouse="mm_", chicken="gg_", bgalgal1="bgalgal1_", human="hs_", gecko="pp_", zebrafish="dr_")

## ---------- Utilities ----------
calculate_sparsity <- function(sparse_matrix) {
  total <- prod(dim(sparse_matrix))
  nonz  <- length(sparse_matrix@x)
  (total - nonz) / total
}

introduce_random_sparsity_sparse <- function(data, sparsity_level) {
  mask <- matrix(runif(length(data)), nrow = nrow(data), ncol = ncol(data)) < sparsity_level
  sparse_data <- data
  sparse_data[mask] <- 0
  sparse_data
}

## ---------- Marker search / plotting ----------
search_marker <- function(species, query, assay = "SCT") {
  name <- deparse(substitute(species))
  col_code <- animal[[names(objects)[objects == name]]]
  species_code <- stringr::str_to_title(names(objects)[objects == name])

  if (species_code == "Gecko") {
    .ensure_file(PB_PATHS$gecko_map, "Gecko map_symbol file not found")
    map_symbol <- read.csv(PB_PATHS$gecko_map)
    map_gecko <- function(genes) {
      mapped <- sapply(genes, function(g) map_symbol$V4[grep(g, map_symbol$GeneSymbol)]) |> unlist()
      sapply(mapped, function(g) rownames(species)[grep(g, rownames(species))]) |> unlist()
    }
    gene_gecko <- map_gecko(query)
    genes <- map_symbol[map_symbol$V4 %in% gene_gecko, c("GeneSymbol","V4")]
    DefaultAssay(species) <- assay
    if (nrow(genes) > 0) {
      for (i in seq_len(nrow(genes))) {
        print(
          FeaturePlot(species, features = genes$V4[i]) +
            ggtitle(paste(genes$GeneSymbol[i], genes$V4[i], sep=" : "), subtitle = species_code) +
            theme(plot.subtitle = element_text(colour = col_code))
        )
      }
    } else {
      message("Marker not present in the dataset (check nomenclature).")
    }
  } else {
    genes <- sapply(query, function(g) rownames(species)[grep(g, rownames(species))]) |> unlist()
    DefaultAssay(species) <- assay
    if (length(genes) > 0) {
      for (g in genes) {
        print(
          FeaturePlot(species, features = g) +
            ggtitle(g, subtitle = species_code) +
            theme(plot.subtitle = element_text(colour = col_code))
        )
      }
    } else {
      message("Marker not present in the dataset (check nomenclature).")
    }
  }
}

plot_markers <- function(species, cluster, assay = "SCT") {
  .ensure_file(PB_PATHS$markers, "Markers file not found")
  markers_df <- read.csv(PB_PATHS$markers, sep = ",", dec = ".")
  name <- deparse(substitute(species))
  col_code <- animal[[names(objects)[objects == name]]]
  species_code <- stringr::str_to_title(names(objects)[objects == name])

  plot_genes <- function(genes, clust_code, species_code) {
    neg_genes <- genes[grepl("-$", genes)] |> stringr::str_replace("-$", "")
    pos_genes <- genes[!grepl("-$", genes)]

    if (species_code == "Gecko") {
      .ensure_file(PB_PATHS$gecko_map, "Gecko map_symbol file not found")
      map_symbol <- read.csv(PB_PATHS$gecko_map)
      map_gecko <- function(g) {
        mapped <- sapply(g, function(x) map_symbol$V4[grep(x, map_symbol$GeneSymbol)]) |> unlist()
        sapply(mapped, function(x) rownames(species)[grep(x, rownames(species))]) |> unlist()
      }
      neg_gecko <- map_gecko(neg_genes)
      pos_gecko <- map_gecko(pos_genes)
      pos <- map_symbol[map_symbol$V4 %in% pos_gecko, c("GeneSymbol","V4")]
      neg <- map_symbol[map_symbol$V4 %in% neg_gecko, c("GeneSymbol","V4")]
      DefaultAssay(species) <- assay
      if (nrow(pos) > 0) for (i in seq_len(nrow(pos))) {
        print(FeaturePlot(species, features = pos$V4[i]) +
                ggtitle(paste(pos$GeneSymbol[i], pos$V4[i], sep=" : "), subtitle = clust_code) +
                theme(plot.subtitle = element_text(colour = col_code)))
      }
      if (nrow(neg) > 0) for (i in seq_len(nrow(neg))) {
        print(FeaturePlot(species, features = neg$V4[i], cols = c("lightgrey","red4")) +
                ggtitle(paste(neg$GeneSymbol[i], neg$V4[i], sep=" : "), subtitle = clust_code) +
                theme(plot.subtitle = element_text(colour = col_code)))
      }
    } else {
      neg <- sapply(neg_genes, function(g) rownames(species)[grep(g, rownames(species))]) |> unlist()
      pos <- sapply(pos_genes, function(g) rownames(species)[grep(g, rownames(species))]) |> unlist()
      DefaultAssay(species) <- assay
      for (g in pos) {
        print(FeaturePlot(species, features = g) +
                ggtitle(g, subtitle = clust_code) +
                theme(plot.subtitle = element_text(colour = col_code)))
      }
      for (g in neg) {
        print(FeaturePlot(species, features = g, cols = c("lightgrey","red4")) +
                ggtitle(g, subtitle = clust_code) +
                theme(plot.subtitle = element_text(colour = col_code)))
      }
    }
  }

  if (cluster == "ALL") {
    cluster_order <- order[[tolower(species_code)]]
    for (clust in cluster_order) {
      clust_code <- paste(species_code, clust, sep=" : ")
      species_markers <- markers_df[grepl(name, colnames(markers_df))]
      genes <- species_markers[species_markers[,1] == clust, 2] |> stringr::str_split(", ") |> unlist()
      plot_genes(genes, clust_code, species_code)
    }
  } else {
    clust_code <- paste(species_code, cluster, sep=" : ")
    species_markers <- markers_df[grepl(name, colnames(markers_df))]
    genes <- species_markers[species_markers[,1] == cluster, 2] |> stringr::str_split(", ") |> unlist()
    plot_genes(genes, clust_code, species_code)
  }
}

## ---------- Orthologue filtering ----------
orthologue_filtering <- function(species1, species2, main = TRUE, sparsity = FALSE) {
  options(future.globals.maxSize = 10 * 1024^3)

  files <- list.files(PB_PATHS$orthologues)
  files <- files[grep(species1, files)]
  files <- files[grep(species2, files)]
  .ensure_file(file.path(PB_PATHS$orthologues, files[1]), "Orthologues file not found")
  orth <- read.table(file.path(PB_PATHS$orthologues, files[1]))
  orth <- orth[, c(species1, species2)]

  rna_files <- list.files(PB_PATHS$objects)
  rna_files <- rna_files[grep("rna", rna_files)]
  file_sp1 <- rna_files[grep(obj[[species1]], rna_files)]
  .ensure_file(file.path(PB_PATHS$objects, file_sp1[1]), "Object file not found")
  sp1 <- readRDS(file.path(PB_PATHS$objects, file_sp1[1]))

  sp1.data <- sp1@assays$RNA@counts
  sp1.data_orth <- merge(x = orth, y = sp1.data, by.x = species1, by.y = "row.names")
  rownames(sp1.data_orth) <- if (isTRUE(main)) sp1.data_orth[, species1] else sp1.data_orth[, species2]
  sp1.data_orth <- sp1.data_orth[, -c(1,2)]
  sp1.data_orth <- Matrix(as.matrix(sp1.data_orth), sparse = TRUE)

  sp1_orth <- CreateSeuratObject(sp1.data_orth, meta.data = sp1@meta.data)

  sp1_orth.list <- SplitObject(sp1_orth, split.by = "orig.ident")
  if (length(sp1_orth.list) > 1) {
    sp1_orth.list <- lapply(sp1_orth.list, SCTransform, vst.flavor = "v1",
                            method = "glmGamPoi", vars.to.regress = vars.out[[species1]])
    features <- SelectIntegrationFeatures(sp1_orth.list, nfeatures = 2000)
    sp1_orth.list <- PrepSCTIntegration(sp1_orth.list, anchor.features = features)
    sp1_orth.list <- lapply(sp1_orth.list, RunPCA, features = features)
    anchors <- FindIntegrationAnchors(sp1_orth.list, reduction = "rpca", dims = 1:50,
                                      normalization.method = "SCT", anchor.features = features)
    sp1_orth <- IntegrateData(anchors, normalization.method = "SCT", dims = 1:50)
    sp1_orth <- RunPCA(sp1_orth, npcs = 50, verbose = FALSE, VariableFeatures(sp1_orth), assay = "integrated")
    sp1_orth <- RunUMAP(sp1_orth, dims = 1:50)
  }

  if (isTRUE(sparsity)) {
    if (!requireNamespace("scWGCNA", quietly = TRUE)) stop("Package 'scWGCNA' required for pseudocells.")
    Idents(sp1_orth) <- "cell_types"
    sp1_orth <- scWGCNA::calculate.pseudocells(
      s.cells = sp1_orth, seeds = 0.3, nn = 30, reduction = "pca", dims = 1:30
    )
    sp1_orth$cell_types <- sp1_orth$orig.cluster
    sp1_orth <- FindVariableFeatures(sp1_orth)
    sp1_orth <- RunPCA(sp1_orth, npcs = 50)
    sp1_orth <- RunUMAP(sp1_orth, dims = 1:50)
  }

  Idents(sp1_orth) <- "cell_types"
  sp1_orth
}

## ---------- GSI correlation ----------
gsicorr_g <- function(species1) {
  DefaultAssay(species1) <- "RNA"
  species1.markers <- FindAllMarkers(species1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  top <- species1.markers %>% group_by(cluster) %>%
    filter(p_val_adj < 1e-9) %>%
    arrange(abs(avg_log2FC)) %>%
    top_n(n = 400) %>% dplyr::select(gene)

  species1_cluster.averages <- AverageExpression(species1, assays = "SCT")
  Idents(species1) <- "SCT_snn_res.0"
  species1_averages <- AverageExpression(species1, assays = "SCT")
  g <- species1_cluster.averages$SCT / species1_averages$SCT[,1]
  g <- g[top$gene, ]
  g
}

gsicorr_plot <- function(species1_g, species2_g) {
  species1 <- substitute(species1_g) |> deparse() |> stringr::str_replace("_g","")
  species2 <- substitute(species2_g) |> deparse() |> stringr::str_replace("_g","")

  s2 <- as.matrix(species2_g[rownames(species2_g) %in% rownames(species1_g), ])
  s2 <- s2[, order[[species2]]]
  colnames(s2) <- paste("S1", colnames(s2), sep="-")
  s1 <- as.matrix(species1_g[rownames(species1_g) %in% rownames(species2_g), ])
  s1 <- s1[, order[[species1]]]
  colnames(s1) <- paste("S2", colnames(s1), sep="-")

  geTable <- merge(s2, s1, by = "row.names")
  rownames(geTable) <- geTable$Row.names
  geTable <- geTable[, -1]

  gsi.cor <- cor(geTable, method = "spearman")

  ## Permutation p-values (compact)
  set.seed(1)
  nPerm <- 100
  shuffled <- replicate(nPerm, {
    sA <- apply(geTable[, 1:ncol(s2)], 1, sample)
    sB <- apply(geTable[, (ncol(s2)+1):ncol(geTable)], 1, sample)
    cor(cbind(t(sA), t(sB)), method = "spearman")
  }, simplify = FALSE)

  ptab <- matrix(NA_real_, ncol = ncol(geTable), nrow = ncol(geTable),
                 dimnames = list(colnames(geTable), colnames(geTable)))
  idx <- combn(seq_len(ncol(geTable)), 2)
  for (i in seq_len(ncol(idx))) {
    cs <- sapply(shuffled, "[", idx[1,i], idx[2,i])
    p <- mean(abs(cs) >= abs(gsi.cor[idx[1,i], idx[2,i]]))
    ptab[idx[1,i], idx[2,i]] <- p
    ptab[idx[2,i], idx[1,i]] <- p
  }
  p_inv <- 1 - ptab

  colnames(gsi.cor) <- stringr::str_sub(colnames(gsi.cor), 4)
  rownames(gsi.cor) <- stringr::str_sub(rownames(gsi.cor), 4)

  mat <- gsi.cor[(ncol(s2)+1):ncol(gsi.cor), 1:ncol(s2)]
  mat_p <- p_inv[(ncol(s2)+1):ncol(gsi.cor), 1:ncol(s2)]

  ## Quick corrplot (visible markers for p > .95)
  corrplot::corrplot(mat, order = "original", tl.pos = "lt", method = "color",
                     tl.col = "black", is.corr = FALSE, tl.cex = 0.7,
                     sig.level = 0.95, insig = "pch", pch = 5, p.mat = mat_p,
                     pch.cex = 0.15, pch.col = "black",
                     mar = c(3,1,5,1),
                     col = colorRampPalette(rev(brewer.pal(10, "RdBu")))(200))

  ## ComplexHeatmap with annotations
  row_colors <- setNames(color[[species1]][order[[species1]]], order[[species1]])
  row_colors <- row_colors[rownames(mat)]
  names(row_colors) <- row_colors
  col_colors <- setNames(color[[species2]][order[[species2]]], order[[species2]])
  col_colors <- col_colors[colnames(mat)]
  names(col_colors) <- col_colors

  col_fun <- colorRamp2(seq(-1, 1, length.out = 200), colorRampPalette(rev(brewer.pal(10,"RdBu")))(200))

  col_ves <- factor(vesicle[[species2]], levels = unique(vesicle[[species2]]))
  row_ves <- factor(vesicle[[species1]], levels = unique(vesicle[[species1]]))
  mk_lab <- function(x) {
    labs <- unique(x)
    labs_list <- apply(data.frame(text = labs, col = c("gold3","steelblue","forestgreen","firebrick","peru","grey50")[seq_along(labs)]),
                       1, as.list)
    setNames(lapply(labs_list, as.data.frame), labs)
  }
  species1_left <- rowAnnotation(
    v1 = anno_textbox(row_ves, mk_lab(row_ves), background_gp = gpar(fill = "white", col = "white")),
    sp1 = color[[species1]], col = list(sp1 = row_colors), show_legend = FALSE
  )
  species2_top <- HeatmapAnnotation(
    v2 = anno_textbox(col_ves, mk_lab(col_ves), background_gp = gpar(fill = "white", col = "white")),
    sp2 = color[[species2]], col = list(sp2 = col_colors), show_legend = FALSE
  )

  annot <- Heatmap(
    mat, name = "GSI Corr", col = col_fun, cluster_columns = FALSE, cluster_rows = FALSE,
    row_split = row_ves, column_split = col_ves,
    top_annotation = species2_top, left_annotation = species1_left,
    row_names_side = "left", column_names_side = "top",
    row_title = NULL, column_title = NULL,
    column_names_gp = gpar(fontsize = 9.25), row_names_gp = gpar(fontsize = 9.25),
    cell_fun = function(j, i, x, y, w, h, fill) {
      if (!is.na(mat_p[i, j]) && mat_p[i, j] > 0.95) grid.points(x, y, pch = 16, size = unit(1, "mm"))
    }
  )
  print(annot)

  out_dir <- file.path(PB_PATHS$out, "corr")
  .ensure_dir(out_dir)
  pdf(file.path(out_dir, paste0("corr_", species1, "_", species2, "_heat.pdf"))); 
  corrplot::corrplot(mat, order="original", tl.pos="lt", method="color", tl.col="black",
                     is.corr=FALSE, tl.cex=0.7, sig.level=0.95, insig="pch", pch=5, p.mat=mat_p,
                     pch.cex=0.15, pch.col="black", mar=c(3,1,5,1),
                     col=colorRampPalette(rev(brewer.pal(10,"RdBu")))(200))
  dev.off()
  pdf(file.path(out_dir, paste0("corr_", species1, "_", species2, "_annot.pdf"))); print(annot); dev.off()

  mat
}

## ---------- Label transfer ----------
labtransfer_plot <- function(query, ref) {
  species1 <- substitute(query) %>% deparse()
  species2 <- substitute(ref) %>% deparse()

  Idents(query) <- "cell_types"; Idents(ref) <- "cell_types"
  query_ct <- order[[species1]]; ref_ct <- order[[species2]]

  pb_anchors <- FindTransferAnchors(query = query, reference = ref, dims = 1:50,
                                    reference.reduction = "pca", normalization.method = "SCT")
  predictions <- TransferData(anchorset = pb_anchors, refdata = ref$cell_types, dims = 1:50)
  query <- AddMetaData(query, metadata = predictions)
  ref <- RunUMAP(ref, dims = 1:50, reduction = "pca", return.model = TRUE)
  query <- MapQuery(anchorset = pb_anchors, reference = ref, query = query,
                    refdata = list(celltype = "cell_types"),
                    reference.reduction = "pca", reduction.model = "umap")

  Idents(ref) <- "cell_types"; levels(ref) <- ref_ct
  ref_plot <- DimPlot(ref, reduction = "umap", label = TRUE, label.size = 3, repel = TRUE) +
    ggtitle("Reference Annotations") + NoLegend() + scale_color_manual(values = color[[species2]])

  Idents(query) <- "predicted.celltype"
  ref_filt_ct <- ref_ct[levels(query) %in% ref_ct]
  levels(query) <- ref_filt_ct
  pred_plot <- DimPlot(query, reduction = "ref.umap", label = TRUE, label.size = 3, repel = TRUE) +
    ggtitle("Predicted Ref Cell Types for Query Cells") + NoLegend() + scale_color_manual(values = color[[species2]])

  Idents(query) <- "cell_types"; levels(query) <- query_ct
  query_plot <- DimPlot(query, reduction = "ref.umap", label = TRUE, label.size = 3, repel = TRUE) +
    ggtitle("Original Cell Types Transferred") + NoLegend() + scale_color_manual(values = color[[species1]])

  mat <- as.matrix(table(query$cell_types, query$predicted.id))
  for (i in seq_len(nrow(mat))) mat[i, ] <- mat[i, ] / sum(mat[i, ])

  if (!identical(colnames(mat), ref_ct)) {
    lack <- ref_ct[!ref_ct %in% colnames(mat)]
    if (length(lack)) {
      lack_df <- matrix(0, nrow = nrow(mat), ncol = length(lack), dimnames = list(rownames(mat), lack))
      mat <- cbind(mat, lack_df)
    }
    mat <- mat[query_ct, ref_ct]
  } else {
    mat <- mat[query_ct, ref_ct]
  }

  pred_data <- merge(query@meta.data %>% dplyr::select(cell_types), predictions, by = 0)
  rownames(pred_data) <- pred_data[,1]
  pred_data <- pred_data %>% dplyr::select(cell_types, predicted.id, prediction.score.max)
  pred_mat <- mat
  for (i in seq_len(ncol(pred_mat))) for (j in seq_len(nrow(pred_mat))) {
    m <- pred_data %>% filter(cell_types == rownames(pred_mat)[j], predicted.id == colnames(pred_mat)[i]) %>%
      summarise(mean(prediction.score.max)) %>% pull()
    pred_mat[j, i] <- ifelse(is.na(m), 0, m)
  }

  row_colors <- setNames(color[[species1]][order[[species1]]], order[[species1]])[rownames(mat)]
  names(row_colors) <- row_colors
  col_colors <- setNames(color[[species2]][order[[species2]]], order[[species2]])[colnames(mat)]
  names(col_colors) <- col_colors

  col_ves <- factor(vesicle[[species2]], levels = unique(vesicle[[species2]]))
  row_ves <- factor(vesicle[[species1]], levels = unique(vesicle[[species1]]))
  mk_lab <- function(x) {
    labs <- unique(x)
    labs_list <- apply(data.frame(text = labs, col = general_ct[seq_along(labs)]), 1, as.list)
    setNames(lapply(labs_list, as.data.frame), labs)
  }
  col_fun_list <- setNames(vector("list", length(levels(row_ves))), levels(row_ves))
  for (i in seq_along(col_fun_list)) col_fun_list[[i]] <- colorRamp2(c(0,1), c("white", general_ct[i]))

  species1_left <- rowAnnotation(
    v1 = anno_textbox(row_ves, mk_lab(row_ves), background_gp = gpar(fill = "white", col = "white")),
    sp1 = color[[species1]], col = list(sp1 = row_colors), show_legend = FALSE
  )
  species2_top <- HeatmapAnnotation(
    v2 = anno_textbox(col_ves, mk_lab(col_ves), background_gp = gpar(fill = "white", col = "white")),
    sp2 = color[[species2]], col = list(sp2 = col_colors), show_legend = FALSE
  )

  annot <- Heatmap(
    mat, name = "% Transf", col = col_fun_list[["N"]], cluster_columns = FALSE, cluster_rows = FALSE,
    row_split = row_ves, column_split = col_ves,
    top_annotation = species2_top, left_annotation = species1_left,
    row_names_side = "left", column_names_side = "top",
    row_title = NULL, column_title = NULL,
    column_names_gp = gpar(fontsize = 9.25), row_names_gp = gpar(fontsize = 9.25),
    rect_gp = gpar(type = "none"),
    cell_fun = function(j, i, x, y, width, height, fill) {
      ind <- row_ves[i]
      grid.rect(x = x, y = y, width = width, height = height,
                gp = gpar(fill = col_fun_list[[ind]](mat[i, j]), col = NA))
    }
  )

  circles <- Heatmap(
    mat, name = "% Transf", col = col_fun_list[["N"]], cluster_columns = FALSE, cluster_rows = FALSE,
    row_split = row_ves, column_split = col_ves,
    top_annotation = species2_top, left_annotation = species1_left,
    row_names_side = "left", column_names_side = "top",
    row_title = NULL, column_title = NULL,
    column_names_gp = gpar(fontsize = 9.25), row_names_gp = gpar(fontsize = 9.25),
    rect_gp = gpar(type = "none"),
    cell_fun = function(j, i, x, y, width, height, fill) {
      ind <- row_ves[i]
      grid.circle(x = x, y = y, r = abs(pred_mat[i, j]) * mean(unit.c(width, height)),
                  gp = gpar(fill = col_fun_list[[ind]](mat[i, j]), col = NA))
    }
  )

  print(annot); print(circles)

  out_dir <- file.path(PB_PATHS$out, "labeltransfer", paste0("ref_", species2))
  .ensure_dir(out_dir)

  pdf(file.path(out_dir, paste0("labeltransfer_", species1, "_", species2, "_dimplot.pdf"))); 
  print(ref_plot); print(pred_plot); print(query_plot); 
  dev.off()

  pdf(file.path(out_dir, paste0("labeltransfer_", species1, "_", species2, "_annot.pdf"))); 
  print(annot); print(circles); 
  dev.off()
}

## ---------- Integration (orthologue-based) ----------
integration <- function(species = c("mouse","chicken","bgalgal1","human","gecko","zebrafish"),
                        reference = "mouse", method = "rpca", fs = 2000) {
  species_names <- species

  species_orth <- lapply(species_names, function(sp) {
    files <- list.files(PB_PATHS$orthologues)
    files <- files[grep(sp, files)]
    files <- files[grep(reference, files)]
    if (sp == reference) files <- files[1]
    .ensure_file(file.path(PB_PATHS$orthologues, files[1]), "Orthologues file not found")
    orth <- read.table(file.path(PB_PATHS$orthologues, files[1]))
    if (sp == reference) orth <- orth[c("Orthogroup", reference)] else orth <- orth[, c(sp, reference)]
    orth
  })
  names(species_orth) <- species_names

  species_objects <- lapply(species_names, function(sp) {
    files <- list.files(PB_PATHS$objects)
    files <- files[grep("rna", files)]
    file_sp <- files[grep(obj[[sp]], files)]
    .ensure_file(file.path(PB_PATHS$objects, file_sp[1]), "Object file not found")
    s <- readRDS(file.path(PB_PATHS$objects, file_sp[1]))
    orth <- species_orth[[sp]]
    mat <- s@assays$RNA@counts
    filt <- merge(x = orth, y = mat, by.x = sp, by.y = "row.names")
    rownames(filt) <- filt[, reference]
    filt <- filt[, -c(1,2)]
    filt <- Matrix(as.matrix(filt), sparse = TRUE)
    so <- CreateSeuratObject(filt, meta.data = s@meta.data)
    so[["RNA"]] <- as(so[["RNA"]], "Assay")
    so
  })
  names(species_objects) <- species_names

  min_cells <- min(vapply(species_objects, function(x) ncol(x), integer(1)))
  set.seed(1)
  species_objects <- lapply(species_objects, function(x) subset(x, cells = sample(colnames(x), min_cells)))

  integration_obj <- list()
  for (nm in names(species_objects)) {
    sp <- species_objects[[nm]]
    sp_list <- SplitObject(sp, split.by = "species")
    sp_list <- lapply(sp_list, SCTransform, method = "glmGamPoi", vst.flavor = "v1",
                      vars.to.regress = vars.out[[nm]])
    integration_obj <- c(integration_obj, sp_list)
  }
  features <- SelectIntegrationFeatures(integration_obj, nfeatures = fs)
  integration_obj <- PrepSCTIntegration(integration_obj, anchor.features = features)
  integration_obj <- lapply(integration_obj, RunPCA, features = features)
  anchors <- FindIntegrationAnchors(integration_obj, reduction = method, dims = 1:50,
                                    normalization.method = "SCT", anchor.features = features)
  integrated <- IntegrateData(anchors, normalization.method = "SCT", dims = 1:50)
  integrated <- RunPCA(integrated, npcs = 50, verbose = FALSE, VariableFeatures(integrated), assay = "integrated")
  integrated <- RunUMAP(integrated, dims = 1:50)
  integrated
}

## ---------- Simple pseudocell helper ----------
pseudocells <- function(rna_obj, nn = 3) {
  all_barcodes <- rna_obj@meta.data %>% tibble::rownames_to_column("cell") %>% dplyr::select(cell, cell_types)
  dicts <- list(cell_types = dplyr::pull(all_barcodes, cell_types, cell))

  emb <- rna_obj@reductions$umap
  nn_obj <- RANN::nn2(emb@cell.embeddings[, 1:2], k = min(nn, nrow(emb@cell.embeddings)), eps = 0)
  set.seed(1)
  smp <- seq_len(nrow(nn_obj$nn.idx))

  df <- purrr::map_dfr(smp, function(idx) {
    tibble::tibble(
      i = idx,
      self_barcode = rownames(emb)[idx],
      self_cell_types = dicts$cell_types[rownames(emb)[idx]],
      others_barcode = rownames(emb)[nn_obj$nn.idx[idx, ]],
      others_cell_types = dicts$cell_types[rownames(emb)[nn_obj$nn.idx[idx, ]]],
      distance = nn_obj$nn.dists[idx, ]
    )
  })

  agg_list <- list()
  df2 <- df
  for (idx in unique(df2$i)) {
    barcodes <- df2$others_barcode[df2$i == idx]
    cell_types <- df2$others_cell_types[df2$i == idx]
    agg_list[[paste0("agg_", idx)]] <- list(barcodes = barcodes, cell_types = cell_types)
  }

  purrr::imap_dfr(agg_list, function(agg, name) {
    c_type <- names(which.max(table(agg$cell_types)))
    tibble::tibble(name = name, cell_types = c_type, barcodes = agg$barcodes)
  })
}

## ---------- Orthologue filtering (aggregated) ----------
orthologue_filtering_agg <- function(species1, species2, main = TRUE, nn = 3) {
  files <- list.files(PB_PATHS$orthologues)
  files <- files[grep(species1, files)]
  files <- files[grep(species2, files)]
  .ensure_file(file.path(PB_PATHS$orthologues, files[1]), "Orthologues file not found")
  orth <- read.table(file.path(PB_PATHS$orthologues, files[1]))[, c(species1, species2)]

  rna_files <- list.files(PB_PATHS$objects)
  rna_files <- rna_files[grep("rna", rna_files)]
  file_sp1 <- rna_files[grep(obj[[species1]], rna_files)]
  .ensure_file(file.path(PB_PATHS$objects, file_sp1[1]), "Object file not found")
  sp1 <- readRDS(file.path(PB_PATHS$objects, file_sp1[1]))

  pcells <- pseudocells(sp1, nn = nn)
  counts <- sp1@assays$RNA@counts
  agg_counts <- NULL
  for (ag in unique(pcells$name)) {
    cells <- pcells %>% dplyr::filter(name == ag) %>% dplyr::pull(barcodes)
    agg <- Matrix::rowSums(counts[, cells, drop = FALSE]) %>% as.matrix()
    colnames(agg) <- ag
    agg_counts <- if (is.null(agg_counts)) agg else cbind(agg_counts, agg)
  }

  sp1.data_orth <- merge(x = orth, y = agg_counts, by.x = species1, by.y = "row.names")
  rownames(sp1.data_orth) <- if (isTRUE(main)) sp1.data_orth[, species1] else sp1.data_orth[, species2]
  sp1.data_orth <- sp1.data_orth[, -c(1,2)]
  sp1.data_orth <- Matrix(as.matrix(sp1.data_orth), sparse = TRUE)

  ## build metadata from pseudocells
  meta_ct <- t(as.data.frame(setNames(lapply(unique(pcells$name), function(nm) {
    names(which.max(table(pcells %>% dplyr::filter(name == nm) %>% dplyr::pull(cell_types))))
  }), unique(pcells$name))))
  colnames(meta_ct) <- "cell_types"

  pcells$orig.ident <- vapply(pcells$barcodes, function(b) strsplit(b, "_")[[1]][1], "")
  meta_oi <- t(as.data.frame(setNames(lapply(unique(pcells$name), function(nm) {
    names(which.max(table(pcells %>% dplyr::filter(name == nm) %>% dplyr::pull(orig.ident))))
  }), unique(pcells$name))))
  colnames(meta_oi) <- "orig.ident"

  sp1_orth <- CreateSeuratObject(sp1.data_orth)
  sp1_orth <- AddMetaData(sp1_orth, metadata = cbind(meta_ct, meta_oi))

  sp1_orth.list <- SplitObject(sp1_orth, split.by = "orig.ident")
  if (length(sp1_orth.list) > 1) {
    sp1_orth.list <- lapply(sp1_orth.list, SCTransform, method = "glmGamPoi")
    features <- SelectIntegrationFeatures(sp1_orth.list, nfeatures = nrow(sp1_orth))
    sp1_orth.list <- PrepSCTIntegration(sp1_orth.list, anchor.features = features)
    sp1_orth.list <- lapply(sp1_orth.list, RunPCA, features = features)
    anchors <- FindIntegrationAnchors(sp1_orth.list, reduction = "rpca", dims = 1:50,
                                      normalization.method = "SCT", anchor.features = features)
    sp1_orth <- IntegrateData(anchors, normalization.method = "SCT", dims = 1:50)
    sp1_orth <- RunPCA(sp1_orth, npcs = 50, verbose = FALSE, VariableFeatures(sp1_orth), assay = "integrated")
    sp1_orth <- RunUMAP(sp1_orth, dims = 1:50)
  }

  Idents(sp1_orth) <- "cell_types"
  sp1_orth
}

## ---------- SAMap plotting (precomputed CSV) ----------
samap_plot <- function(species1, species2) {
  sp1_ct <- order[[species1]]; sp2_ct <- order[[species2]]
  sp1 <- siglas[[species1]];   sp2 <- siglas[[species2]]

  path <- file.path(PB_PATHS$out, "samap", paste0("samap_", species1, "_", species2, "_coef.csv"))
  .ensure_file(path, "SAMap coefficient file not found")
  corr <- read.csv(path, row.names = 1)

  corr_tx <- corr[grep(sp1, rownames(corr)), grep(sp2, colnames(corr))]
  rownames(corr_tx) <- str_replace(rownames(corr_tx), sp1, "")
  colnames(corr_tx) <- str_replace(colnames(corr_tx), sp2, "")
  mat <- as.matrix(corr_tx[sp1_ct, sp2_ct])

  row_colors <- setNames(color[[species1]][order[[species1]]], order[[species1]])[rownames(mat)]
  names(row_colors) <- row_colors
  col_colors <- setNames(color[[species2]][order[[species2]]], order[[species2]])[colnames(mat)]
  names(col_colors) <- col_colors

  col_ves <- factor(vesicle[[species2]], levels = unique(vesicle[[species2]]))
  row_ves <- factor(vesicle[[species1]], levels = unique(vesicle[[species1]]))
  mk_lab <- function(x) {
    labs <- unique(x)
    labs_list <- apply(data.frame(text = labs, col = general_ct[seq_along(labs)]), 1, as.list)
    setNames(lapply(labs_list, as.data.frame), labs)
  }
  col_fun_list <- setNames(vector("list", length(levels(row_ves))), levels(row_ves))
  for (i in seq_along(col_fun_list)) col_fun_list[[i]] <- colorRamp2(c(0,1), c("white", general_ct[i]))

  species1_left <- rowAnnotation(
    v1 = anno_textbox(row_ves, mk_lab(row_ves), background_gp = gpar(fill = "white", col = "white")),
    sp1 = color[[species1]], col = list(sp1 = row_colors), show_legend = FALSE
  )
  species2_top <- HeatmapAnnotation(
    v2 = anno_textbox(col_ves, mk_lab(col_ves), background_gp = gpar(fill = "white", col = "white")),
    sp2 = color[[species2]], col = list(sp2 = col_colors), show_legend = FALSE
  )

  annot <- Heatmap(
    mat, name = "Similarity Coef", col = col_fun_list[["N"]],
    cluster_columns = FALSE, cluster_rows = FALSE,
    row_split = row_ves, column_split = col_ves,
    top_annotation = species2_top, left_annotation = species1_left,
    row_names_side = "left", column_names_side = "top",
    row_title = NULL, column_title = NULL,
    column_names_gp = gpar(fontsize = 9.25), row_names_gp = gpar(fontsize = 9.25),
    rect_gp = gpar(type = "none"),
    cell_fun = function(j, i, x, y, width, height, fill) {
      ind <- row_ves[i]
      grid.rect(x = x, y = y, width = width, height = height,
                gp = gpar(fill = col_fun_list[[ind]](mat[i, j]), col = NA))
    }
  )
  print(annot)

  out_dir <- file.path(PB_PATHS$out, "samap")
  .ensure_dir(out_dir)
  pdf(file.path(out_dir, paste0("samap_", species1, "_", species2, "_annot.pdf"))); print(annot); dev.off()
}

## ---------- Auto-correlation (within species) ----------
autocorr <- function(species1, assay = "integrated") {
  rna_files <- list.files(PB_PATHS$objects)
  rna_files <- rna_files[grep("rna", rna_files)]
  file_sp1 <- rna_files[grep(obj[[species1]], rna_files)]
  .ensure_file(file.path(PB_PATHS$objects, file_sp1[1]), "Object file not found")
  sp1 <- readRDS(file.path(PB_PATHS$objects, file_sp1[1]))

  Idents(sp1) <- "cell_types"
  avg_exp <- AverageExpression(sp1, assays = assay, return.seurat = FALSE)[[assay]]
  cor_matrix <- cor(avg_exp, method = "pearson")
  mat <- cor_matrix[order[[species1]], order[[species1]]]

  row_colors <- setNames(color[[species1]][order[[species1]]], order[[species1]])[rownames(mat)]
  names(row_colors) <- row_colors
  col_colors <- setNames(color[[species1]][order[[species1]]], order[[species1]])[colnames(mat)]
  names(col_colors) <- col_colors

  col_fun <- colorRamp2(seq(-1, 1, length.out = 200), colorRampPalette(rev(brewer.pal(10,"RdBu")))(200))

  col_ves <- factor(vesicle[[species1]], levels = unique(vesicle[[species1]]))
  row_ves <- factor(vesicle[[species1]], levels = unique(vesicle[[species1]]))
  mk_lab <- function(x) {
    labs <- unique(x)
    labs_list <- apply(data.frame(text = labs, col = c("gold3","steelblue","forestgreen","firebrick","peru","grey50")[seq_along(labs)]),
                       1, as.list)
    setNames(lapply(labs_list, as.data.frame), labs)
  }

  species1_left <- rowAnnotation(
    v1 = anno_textbox(row_ves, mk_lab(row_ves), background_gp = gpar(fill = "white", col = "white")),
    sp1 = color[[species1]], col = list(sp1 = row_colors), show_legend = FALSE
  )
  species2_top <- HeatmapAnnotation(
    v2 = anno_textbox(col_ves, mk_lab(col_ves), background_gp = gpar(fill = "white", col = "white")),
    sp2 = color[[species1]], col = list(sp2 = col_colors), show_legend = FALSE
  )

  annot <- Heatmap(
    mat, name = "AutoCorr", col = col_fun,
    cluster_columns = FALSE, cluster_rows = FALSE,
    row_split = row_ves, column_split = col_ves,
    top_annotation = species2_top, left_annotation = species1_left,
    row_names_side = "left", column_names_side = "top",
    row_title = NULL, column_title = NULL,
    column_names_gp = gpar(fontsize = 9.25), row_names_gp = gpar(fontsize = 9.25)
  )
  print(annot)

  out_dir <- file.path(PB_PATHS$out, "autocorr")
  .ensure_dir(out_dir)
  pdf(file.path(out_dir, paste0(species1, "_autocorr.pdf"))); print(annot); dev.off()
}
