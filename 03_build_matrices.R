library(SummarizedExperiment)
library(edgeR)
library(matrixStats)

se_expr <- readRDS("data/se_expr_raw.rds")
se_meth <- readRDS("data/se_meth_raw.rds")

N_TOP_GENES <- 5000
N_TOP_CPGS <- 5000
dir.create("results", showWarnings = FALSE)

bc_field <- function(x, i) vapply(strsplit(x, "-"), `[`, character(1), i)

case_id <- function(x) {
  paste(bc_field(x, 1), bc_field(x, 2), bc_field(x, 3), sep = "-")
}

select_primary_one_per_case <- function(se) {
  bc <- colnames(se)
  is_primary <- substr(bc_field(bc, 4), 1, 2) == "01"
  se <- se[, is_primary]
  cid <- case_id(colnames(se))
  ord <- order(cid, colnames(se))
  se <- se[, ord]
  cid <- cid[ord]
  se <- se[, !duplicated(cid)]
  colnames(se) <- cid[!duplicated(cid)]
  se
}

se_expr <- select_primary_one_per_case(se_expr)
se_meth <- select_primary_one_per_case(se_meth)

common_cases <- sort(intersect(colnames(se_expr), colnames(se_meth)))
se_expr <- se_expr[, common_cases]
se_meth <- se_meth[, common_cases]

message("primary-only cases: expr=", ncol(se_expr),
        " meth=", ncol(se_meth),
        " intersection=", length(common_cases))

counts <- assay(se_expr, "unstranded")
gene_meta <- rowData(se_expr)

keep_type <- gene_meta$gene_type == "protein_coding"
counts <- counts[keep_type, ]
gene_meta <- gene_meta[keep_type, ]

keep_expr <- rowSums(cpm(counts) > 1) >= 0.2 * ncol(counts)
counts <- counts[keep_expr, ]
gene_meta <- gene_meta[keep_expr, ]

logcpm <- cpm(counts, log = TRUE, prior.count = 2)
rownames(logcpm) <- gene_meta$gene_name
logcpm <- logcpm[!is.na(rownames(logcpm)) & !duplicated(rownames(logcpm)), ]

gene_var <- rowVars(logcpm)
expr_mat <- logcpm[order(gene_var, decreasing = TRUE)[seq_len(min(N_TOP_GENES, nrow(logcpm)))], ]

beta <- assay(se_meth)
probe_chr <- as.character(seqnames(rowRanges(se_meth)))

keep_probe <- !is.na(probe_chr) &
  !probe_chr %in% c("chrX", "chrY", "X", "Y") &
  rowSums(is.na(beta)) == 0

beta <- beta[keep_probe, ]
beta <- pmin(pmax(beta, 0.001), 0.999)
mval <- log2(beta / (1 - beta))

cpg_var <- rowVars(mval)
meth_mat <- mval[order(cpg_var, decreasing = TRUE)[seq_len(min(N_TOP_CPGS, nrow(mval)))], ]

stopifnot(identical(colnames(expr_mat), colnames(meth_mat)))

mofa_input <- list(expression = expr_mat, methylation = meth_mat)
saveRDS(mofa_input, "data/mofa_input.rds")

clin <- colData(se_expr)
saveRDS(clin, "data/clinical_colData.rds")

layer_summary <- data.frame(
  view = c("expression", "methylation"),
  n_features = c(nrow(expr_mat), nrow(meth_mat)),
  n_samples = c(ncol(expr_mat), ncol(meth_mat)),
  features_before_selection = c(nrow(logcpm), nrow(mval))
)
write.csv(layer_summary, "results/layer_summary.csv", row.names = FALSE)
print(layer_summary)
