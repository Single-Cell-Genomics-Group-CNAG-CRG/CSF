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

