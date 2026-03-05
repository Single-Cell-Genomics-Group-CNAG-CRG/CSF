suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(magrittr)
  library(scRepertoire)
  library(patchwork)
  library(rscrublet) # devtools::install_github("iaaaka/Rscrublet")
  library(Matrix)
  library(glue)
  library(ggpubr)
  library(cowplot)
  library(ggExtra)
  library(harmony)
})

# -------------------------------

# function to run R scrublet on a Seurat object
# assumes you have RNA counts on your object (@assays$RNA@counts)
# returns the same Seurat object with a "doublet_score" metadata column

run_scrublet <- function(seu_obj){
  count_matrix = t(as(seu_obj@assays$RNA@counts,'dgTMatrix'))
  scrr = scrub_doublets(E_obs = count_matrix,
                        expected_doublet_rate=0.06,
                        min_counts=3,
                        min_cells=3,
                        min_gene_variability_pctl=85,
                        n_prin_comps=30)

  scrr=call_doublets(scrr)
  #plot_doublet_histogram(scrr)
  seu_obj$doublet_score = scrr$doublet_scores_obs
  seu_obj
}

# -------------------------------

# function to load a CellRanger library
# input: the expression directories and T/BCR if any
# also project and sample names to label/identify the object
# sp parameter is H or m to indicate mouse or human
# run scrublet can be skipped with the run_scrubblet arg
# it returns a Seurat object

load_seurat <- function(gex_dir, tcr_dir = NULL, bcr_dir = NULL, project, sample, sp = "H", run_scrublet = TRUE){

  print(paste0("Processing sample ", sample, "..."))

  #### create Seurat object -----------------

  obj <- Seurat::CreateSeuratObject(
    Read10X(
      data.dir = gex_dir,
      strip.suffix = TRUE),
    min.cells = 5,
    min.features = 5
  )

  # add metadata to object
  obj$project <- project
  obj$sample <- sample
  # compute mitochondrial percentage
  obj[["percent.mt"]] <- ifelse(sp == "H", PercentageFeatureSet(obj, pattern = "^MT-"), PercentageFeatureSet(obj, pattern = "^mt-"))
  # rename cells, in this case by appending the sample id
  obj <- RenameCells(obj, add.cell.id = paste0(project, "_", sample))

  print(obj)

  #### add TCR info (metadata) -----------------

  if (!is.null(tcr_dir)){
    if (file.exists(glue::glue("{tcr_dir}filtered_contig_annotations.csv"))){
      print(paste0("Sample ", sample, " has TCRs!"))
      print(paste0("Reading ", glue::glue("{tcr_dir}filtered_contig_annotations.csv")))
      tcr <- read.csv(glue::glue("{tcr_dir}filtered_contig_annotations.csv"))
      tcr <- scRepertoire::combineTCR(tcr, samples = project, ID = sample)
      obj <- scRepertoire::combineExpression(
        tcr,
        obj,
        proportion = FALSE,
        cloneTypes = c(Single = 1, Small = 10, Medium = 100, Large = 1000, Hyperexpanded = 10000)
      )
    }
    else {print(paste0(glue::glue("{tcr_dir}filtered_contig_annotations.csv"), " does not exist!"))}
  }
  else {print(paste0("Sample ", sample, " does NOT have TCRs!"))}


  #### add BCR info (metadata) -----------------

  if (!is.null(bcr_dir)){
    if (file.exists(glue::glue("{bcr_dir}filtered_contig_annotations.csv"))){
      print(paste0("Sample ", sample, " has BCRs!"))
      print(paste0("Reading ", glue::glue("{bcr_dir}filtered_contig_annotations.csv")))
      bcr <- read.csv(glue::glue("{bcr_dir}filtered_contig_annotations.csv"))
      bcr <- scRepertoire::combineBCR(bcr, samples = project, ID = sample)
      obj <- scRepertoire::combineExpression(
        bcr,
        obj,
        proportion = FALSE,
        cloneTypes = c(Single = 1, Small = 10, Medium = 100, Large = 1000, Hyperexpanded = 10000)
      )
    }
    else {print(paste0(glue::glue("{bcr_dir}filtered_contig_annotations.csv"), " does not exist!"))}
  }
  else {print(paste0("Sample ", sample, " does NOT have BCRs!"))}

  #### run scrublet (R port, from function I created) -----------------

  if (run_scrublet){

    print(paste0("Running Scrubblet for sample ", sample, "..."))
    obj <- run_scrublet(obj)

  }
  else {print(paste0("Skipping Scrubblet for sample ", sample))}

  print("Done!")

  obj
}

# -------------------------------

# COLOR PALETTE DEFINITION

pal_lv2 <- list(
  "B cells" = "#a0cfd9",
  "Plasma cells" = "#578797",
  "CD8 T cells" = "#4f8a65",
  "CD4 T cells" = "#9cd379",
  "NK cells" = "#b7b5e2",
  "Macrophages" = "#d96f6f",
  "Monocytes" = "#fcbe81",
  "DC" = "#fad57f",
  "Non-immune" = "#94735e"
)

pal_disease <- list(
  "Brain met" = "#456c2c",
  "Brain Metastasis" = "#456c2c",
  "Glioblastoma" = "#788fa3",
  "Inflammatory" = "#7d1517",
  "Lymphoma" = "#cc8630",
  "Healthy" = "black"
)

alt_pal_lv2 <- list(
  "B cells" = "#88CCEE",
  "Plasma cells" = "#332288",
  "CD8 T cells" = "#CC6677",
  "CD4 T cells" = "#AA4499",
  "NK cells" = "#882255",
  "Macrophages" = "#117733",
  "Monocytes" = "#44AA99",
  "DC" = "#DDCC77",
  "Non-immune" = "black"
)

pal_myeloid <- list('BAMs' = '#d7cd95',
                    'Microglia-like' = '#cc4566',
                    'Macrophages anti-inflammatory'='#b5da4c',
                    'Macrophages MT-RTM-like' = '#5d262a',
                    'Macrophages proliferative' = '#dabb43',
                    'Monocytes intermediate' = '#dc4733',
                    'Monocytes classical' = '#778632',
                    'DC mreg' = '#924026',
                    'DC1' = '#a38c73',
                    'DC2' = '#cf8238',
                    'DC5' = '#554825',
                    'pDC' = '#d48980'
)

pal_t <- list('CD4 CM' = '#61c271',
              'CD4 IFN Response' = '#6946c9',
              'CD4 Naive' = '#a1dc49',
              'CD4 T helper' = '#c454ca',
              'CD4 T reg' = '#cca83e',
              'CD4 Th17' = '#4a2c70',
              'CD8 Cytotoxic' = '#cad09a',
              'CD8 EM' = '#d64b83',
              'CD8 Exhausted' = '#7fd5cf',
              'CD8 Pre-exhausted' = '#d75332',
              'NK-gd' = '#6c82c4',
              'T cells Proliferative' = '#556931'
)

pal_clones <- c("#F0F921", "#F69441", "#CA4778", "#7D06A5", "#0D0887")

# gene expression
pal_gene_exp <- c("#ADD8E633", "#E46726")

pal_type <- list(
  "exp both" = "forestgreen",
  "exp AT" =  "orange1",
  "exp BT" = "darkblue",
  "NE" = "black"
)

pal_xenium <- list(
  "B and Plasma cells" = "#a0cfd9",
  "T cells" = "#9cd379",
  "Macrophages" = "#d96f6f",
  "Microglia" = "#fad57f",
  "Non-immune" = "#94735e",
  "Neutrophils" = "purple4",
  "Tumor" = 'grey'
)

pal_k_l <- list("IGKC+" = "#d6604d", "IGKC-" = "#4393c3", `NA` = "black")


# -------------------------------


marker_genes <- list(
  "CD4 Naive/CM" = c("CD4", "ANXA1", "PASK", "SELL", "LEF1", "NOSIP", "CCR7", "TCF7", "ACTN1", "FOXP1", "KLF2", "ITGA6", "CD8A-", "CD8B-", "GZMK-"),
  "CD4 Effector/Mem" = c("CD4", "ZNF683", "KLRB1", "PRDM1", "CX3CR1", "EOMES", "KLRG1", "TNFSF13B", "GZMK", "CCL5", "CCL4", "NKG7", "CD69", "ITGAE", "CD8A-", "CD8B-"),
  "T helper" = c("CD4", "CXCR3", "GATA3", "RORC", "RORA", "IL17F", "IL17A", "CCR6", "CXCR6", "IFNG", "IL4", "IL6ST", "CXCR5", "CXCL13", "PDCD1", "CD8A-", "CD8B-"),
  "CD4 IFN response" = c("CD4", "IFI16", "IFI35", "IFI44", "IFI44L", "IFI6", "IFIH1", "IFIT1", "IFIT2", "IFIT3", "IFIT5", "ISG15", "CD8A-", "CD8B-"),
  "CD4 Proliferative" = c("CD4", "MKI67", "TOP2A", "STMN1", "UBE2C", "PCLAF", "CENPF", "CDK1", "CD8A-", "CD8B-"),
  "T reg" = c("IL32", "CCR7", "LEF1", "TCF7", "FOXP3", "CTLA4", "IL2RA", "ICOS", "TIGIT", "TOX2", "IKZF2", "GATA3", "CD28", "CD8A-", "CD8B-"),
  "Gamma Delta" = c("TRGC1", "TRGC2", "TRDC", "CD8A-", "CD8B-", "CD4-"),
  # "MAIT" = c("KLRB1, IL7R", "SLC4A10"),
  "CD8 Naive/CM" = c("CD4-", "ANXA1", "PASK", "SELL", "LEF1", "NOSIP", "CCR7", "TCF7", "ACTN1", "FOXP1", "KLF2", "ITGA6", "CD8A", "CD8B", "GZMK-"),
  "CD8 Mem" = c("CD8A", "CD8B", "ZNF683", "KLRB1", "PRDM1", "CX3CR1", "EOMES", "KLRG1", "TNFSF13B", "CD4-"),
  "CD8 Cytotoxic" = c("CD8A", "CD8B", "GZMK", "GZMH", "CCL5", "CCL4", "CD69", "PRF1", "ITGAE", "CD4-", "CST7", "GZMA", "CCL4L2", "CTSW", "GZMH", "GZMM", "HLA-C"),
  "CD8 IFN response" = c("CD8A", "CD8B", "IFI16", "IFI35", "IFI44", "IFI44L", "IFI6", "IFIH1", "IFIT1", "IFIT2", "IFIT3", "IFIT5", "ISG15", "CD4-"),
  "CD8 Exhausted" = c("CD8A", "CD8B", "HAVCR2", "LAG3", "PDCD1", "TIGIT", "TOX", "TOX2", "LAYN", "CTLA4", "CD4-"),
  "CD8 Proliferative" = c("CD8A", "CD8B", "MKI67", "TOP2A", "STMN1", "UBE2C", "PCLAF", "CENPF", "CDK1", "CD4-"),
  # "ILC" = c("KIT", "NCR1", "KLRG1"),
  "NK" = c("NCAM1", "FCGR3A", "CX3CR1", "GNLY", "KLRC2", "KLRD1", "KLRC3", "KLRK1", "KLRC1", "NKG7", "XCL2", "KLRB1", "PRF1", "TRDC"),
  # "Immature B cell" = c("MS4A1", "CD79A", "CD19", "RAG1", "RAG2", "CD3E-", "CD3G-", "CD3D-", "CD4-", "CD8A-", "CD8B-"),
  "Naive B cell" = c("MS4A1", "IGHD", "IGHM", "CCR7", "SELL", "TCL1A", "CD79A", "VPREB3", "FCRL1", "NIBAN3", "CD79B", "HVCN1", "CD72", "FCER2", "CD83", "CD19", "CD3E-", "CD3G-", "CD3D-", "CD4-", "CD8A-", "CD8B-"),
  "Memory B cell" = c("CD79A", "MS4A1", "CD27", "TNFRSF13B", "ITGAX", "PRDM1", "CD24", "BANK1", "CD74", "HLA-DRA", "IGHA1", "BLK", "SPIB", "P2RX5", "IGHA2", "CD37", "CD3E-", "CD3G-", "CD3D-", "CD4-", "CD8A-", "CD8B-"),
  "Plasma cells" = c("MZB1", "SDC1", "IGHG1", "JCHAIN", "IGHA1", "IGHG3", "IGLC3", "IGLC1", "IGHGP", "DERL3", "IGHG4", "XBP1", "IRF4", "CD3E-", "CD3G-", "CD3D-", "CD4-", "CD8A-", "CD8B-"),
  "Monocytes" = c("CD14", "S100A8", "S100A9", "LYZ", "VCAN", "FCN1"),
  "M1 Macrophages" = c("HLA-DPB1", "HLA-DPA1", "HLA-DQA1", "HLA-DQB1", "HLA-DQA2", "HLA-DMA", "HLA-DRB5", "HLA-DRB1", "HLA-DRA", "HLA-DMB", "HLA-DQB2", "APOE", "APOC1", "CD68", "C1QA","C1QB", "C1QC","CCL2","IL1B","CCL4","CCL7","CCL8","NFKB","CD40", "CXCL2", "CXCL3", "CXCL9", "CXCL10","CXCL11","IDO1","NFKBIA", "TNF","CXCL8","G0S2","IL6","INHBA", "CD14-", "LYZ-", "VCAN-", "FCN1-"),
  "M2 Macrophages" =  c("APOE", "APOC1", "CD68", "C1QA","C1QB", "C1QC","CD68", "SELENOP", "MRC1", "CCL18","CD163", "CD209", "ARG1", "IL10", "CD274","CHIT1", "RNASE1", "TREM2", "IL10", "ITGA4", "LGALS9", "MARCO", "TGFB2", "TGFB1", "CSF1R", "CSF1", "SPP1","TREM2", "CD14-", "LYZ-", "VCAN-", "FCN1-"),
  "Myeloid proliferative" = c("CD68", "CD163", "MKI67", "TOP2A", "STMN1", "UBE2C", "PCLAF", "CENPF", "CDK1"),
  "Alveolar macrophages" = c("GPNMB", "SPP1", "CTSB", "C1QC", "C1QB", "APOC1", "APOE", "GLUL", "C1QA", "HMOX1", "FTL", "FN1", "PLTP", "MARCO", "CD163", "CD68", "CTSL",   "TREM2",  "TMIGD3", "FCGRT",  "CTSD"),
  "pDC" = c("IL3RA", "IRF7", "LILRA4", "IRF8", "JCHAIN", "GZMB"),
  "DC1" = c("CLEC9A", "XCR1", "IDO1", "CLNK", "ZNF366"),
  "DC2" = c("CD1C", "FCER1A", "CLEC10A"),
  "DC3" = c("CD1C", "S100A8", "S100A9", "ANXA1"),
  "DC4" = c("ITGAX", "FCGR3A", "SERPINA1", "LILRB2", "SIGLEC10"),
  "DC5" = c("AXL", "SIGLEC6", "CD22", "DAB2"),
  "Mesothelial cells" = c("UPK3B", "KRT7", "CDH2", "PECAM1", "PRG4")
)
