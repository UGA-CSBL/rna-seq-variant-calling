# Follow installation tutorials from github.com/PeeperLab/XenofilteR
# devtools::install_github("PeeperLab/XenofilteR")

# Environment setup -----
library(tidyverse)
library(XenofilteR)
bp.param <- SnowParam(workers = 20, type = "SOCK")

work_dir <- "/home/data/mutation"
res_dir <- str_glue("{work_dir}/bam/filtered/")
dir.create(res_dir, showWarnings = F, recursive = T)

# Sample matrix -----
design_mat <- data.table::fread(str_glue("{work_dir}/design_matrix.csv")) %>%
  filter(FileStatus == "Good")

sampleIDs <- design_mat$SampleName

sample.list <- design_mat %>%
  mutate(
    graft_bam_file = str_glue("{work_dir}/bam/human/{SampleName}Aligned.sortedByCoord.out.bam"),
    host_bam_file = str_glue("{work_dir}/bam/mouse/{SampleName}Aligned.sortedByCoord.out.bam")
  ) %>%
  dplyr::select(graft_bam_file, host_bam_file)

# Filter reads -----
XenofilteR(
  sample.list = sample.list,
  destination.folder = res_dir,
  bp.param = bp.param,
  output.names = sampleIDs
)

# sessionInfo() -----
# R version 4.0.3 (2020-10-10)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04 LTS
#
# Matrix products: default
# BLAS/LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.8.so
#
# locale:
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8
# [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=C
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C
# [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C
#
# attached base packages:
# [1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base
#
# other attached packages:
# [1] XenofilteR_1.6       Rsamtools_2.6.0      Biostrings_2.58.0    XVector_0.30.0
# [5] GenomicRanges_1.42.0 GenomeInfoDb_1.26.2  IRanges_2.24.1       S4Vectors_0.28.1
# [9] BiocGenerics_0.36.0  BiocParallel_1.24.1  forcats_0.5.1        stringr_1.4.0
# [13] dplyr_1.0.4          purrr_0.3.4          readr_1.4.0          tidyr_1.1.2
# [17] tibble_3.1.0         ggplot2_3.3.3        tidyverse_1.3.0
#
# loaded via a namespace (and not attached):
# [1] Rcpp_1.0.6                  lubridate_1.7.9.2           lattice_0.20-41
# [4] assertthat_0.2.1            utf8_1.1.4                  R6_2.5.0
# [7] cellranger_1.1.0            futile.options_1.0.1        backports_1.2.1
# [10] reprex_1.0.0                httr_1.4.2                  pillar_1.5.0
# [13] zlibbioc_1.36.0             rlang_0.4.10                readxl_1.3.1
# [16] rstudioapi_0.13             data.table_1.14.0           Matrix_1.3-2
# [19] RCurl_1.98-1.2              munsell_0.5.0               DelayedArray_0.16.1
# [22] broom_0.7.5                 compiler_4.0.3              modelr_0.1.8
# [25] pkgconfig_2.0.3             tidyselect_1.1.0            SummarizedExperiment_1.20.0
# [28] GenomeInfoDbData_1.2.4      matrixStats_0.58.0          fansi_0.4.2
# [31] withr_2.4.1                 crayon_1.4.1                dbplyr_2.1.0
# [34] GenomicAlignments_1.26.0    bitops_1.0-6                grid_4.0.3
# [37] jsonlite_1.7.2              gtable_0.3.0                lifecycle_1.0.0
# [40] DBI_1.1.1                   formatR_1.7                 magrittr_2.0.1
# [43] scales_1.1.1                cli_2.3.1                   stringi_1.5.3
# [46] fs_1.5.0                    futile.logger_1.4.3         xml2_1.3.2
# [49] ellipsis_0.3.1              generics_0.1.0              vctrs_0.3.6
# [52] lambda.r_1.2.4              tools_4.0.3                 Biobase_2.50.0
# [55] glue_1.4.2                  hms_1.0.0                   MatrixGenerics_1.2.1
# [58] colorspace_2.0-0            rvest_0.3.6                 haven_2.3.1
