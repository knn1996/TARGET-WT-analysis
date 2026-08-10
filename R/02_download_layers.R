library(TCGAbiolinks)
library(SummarizedExperiment)

project <- "TARGET-WT"
dir.create("data", showWarnings = FALSE)

q_expr <- GDCquery(
  project = project,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  access = "open"
)

q_meth <- GDCquery(
  project = project,
  data.category = "DNA Methylation",
  data.type = "Methylation Beta Value",
  platform = "Illumina Human Methylation 450",
  access = "open"
)

GDCdownload(q_expr, files.per.chunk = 20)
GDCdownload(q_meth, files.per.chunk = 5)

se_expr <- GDCprepare(q_expr, summarizedExperiment = TRUE)
se_meth <- GDCprepare(q_meth, summarizedExperiment = TRUE)

saveRDS(se_expr, "data/se_expr_raw.rds")
saveRDS(se_meth, "data/se_meth_raw.rds")

gdc_release <- getGDCInfo()
writeLines(
  c(
    paste("gdc_release:", gdc_release$data_release),
    paste("date_downloaded:", Sys.Date()),
    paste("expr_dim:", paste(dim(se_expr), collapse = " x ")),
    paste("meth_dim:", paste(dim(se_meth), collapse = " x "))
  ),
  "data/download_provenance.txt"
)

print(dim(se_expr))
print(dim(se_meth))
