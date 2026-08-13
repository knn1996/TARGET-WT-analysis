library(MOFA2)
library(SummarizedExperiment)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
source("barcodes.R")

model <- readRDS("models/mofa_wt_main.rds")
mofa_input <- readRDS("data/mofa_input.rds")
se_expr <- readRDS("data/se_expr_raw.rds")
se_meth <- readRDS("data/se_meth_raw.rds")
dir.create("results", showWarnings = FALSE)

N_TOP <- 100
VAR_MIN <- 1
TECH_CUTOFF <- 0.5
promoter_groups <- c("TSS200", "TSS1500", "5'UTR", "1stExon")

mofa_cases <- colnames(mofa_input$expression)

raw_detect <- function(se, assayname, count) {
  se <- se[, is_primary_sample(colnames(se))]
  cid <- case_id(colnames(se))
  ord <- order(cid, colnames(se))
  se <- se[, ord]; cid <- cid[ord]
  se <- se[, !duplicated(cid)]
  colnames(se) <- cid[!duplicated(cid)]
  se <- se[, mofa_cases]
  m <- assay(se, assayname)
  if (count) colSums(m > 0) else colSums(!is.na(m))
}

det <- cbind(detect_expr = raw_detect(se_expr, "unstranded", TRUE),
             detect_meth = raw_detect(se_meth, 1, FALSE))
Z <- get_factors(model)$group1
detcor <- round(cor(Z, det), 2)
write.csv(detcor, "results/factor_detection_correlation.csv")
print(detcor)

r2 <- get_variance_explained(model)$r2_per_factor$group1
detmax <- apply(abs(detcor), 1, max)[rownames(r2)]
technical <- detmax > TECH_CUTOFF

shared_mask  <- r2[, "expression"] >= VAR_MIN & r2[, "methylation"] >= VAR_MIN & !technical
private_mask <- r2[, "methylation"] >= VAR_MIN & r2[, "expression"] < VAR_MIN & !technical
stopifnot(any(shared_mask), any(private_mask))

shared_factor <- rownames(r2)[shared_mask][which.max(apply(r2[shared_mask, , drop = FALSE], 1, min))]
private_factor <- rownames(r2)[private_mask][which.max(r2[private_mask, "methylation"])]

selection <- data.frame(
  role = c("shared", "meth_private"),
  factor = c(shared_factor, private_factor),
  expr_r2 = c(r2[shared_factor, "expression"], r2[private_factor, "expression"]),
  meth_r2 = c(r2[shared_factor, "methylation"], r2[private_factor, "methylation"]),
  detect_cor = c(detmax[shared_factor], detmax[private_factor])
)
write.csv(selection, "results/factor_selection.csv", row.names = FALSE)
print(selection)

W_expr <- get_weights(model, views = "expression")[[1]]
W_meth <- get_weights(model, views = "methylation")[[1]]

top_feat <- function(W, factor, n) {
  v <- W[, factor]
  ord <- order(abs(v), decreasing = TRUE)[seq_len(n)]
  data.frame(feature = rownames(W)[ord], weight = round(v[ord], 4),
             stringsAsFactors = FALSE)
}

annotate_probes <- function(df) {
  a <- anno[df$feature, c("chr", "pos", "Relation_to_Island",
                          "UCSC_RefGene_Name", "UCSC_RefGene_Group")]
  cbind(df, as.data.frame(a, stringsAsFactors = FALSE))
}

is_promoter <- function(grp) {
  vapply(strsplit(grp, ";"), function(g) any(g %in% promoter_groups), logical(1))
}

anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

expr_shared  <- top_feat(W_expr, shared_factor, N_TOP)
meth_shared  <- annotate_probes(top_feat(W_meth, shared_factor, N_TOP))
meth_private <- annotate_probes(top_feat(W_meth, private_factor, N_TOP))
write.csv(expr_shared,  "results/f_shared_top_expression.csv",  row.names = FALSE)
write.csv(meth_shared,  "results/f_shared_top_methylation.csv",  row.names = FALSE)
write.csv(meth_private, "results/f_private_top_methylation.csv", row.names = FALSE)

island_levels <- c("Island", "N_Shore", "S_Shore", "N_Shelf", "S_Shelf", "OpenSea")
context_table <- rbind(
  shared  = table(factor(meth_shared$Relation_to_Island,  levels = island_levels)),
  private = table(factor(meth_private$Relation_to_Island, levels = island_levels))
)
write.csv(context_table, "results/probe_context_distribution.csv")
print(context_table)

expr_w <- W_expr[, shared_factor]
mp <- data.frame(
  probe  = rownames(W_meth),
  weight = W_meth[, shared_factor],
  island = anno[rownames(W_meth), "Relation_to_Island"],
  gene   = anno[rownames(W_meth), "UCSC_RefGene_Name"],
  group  = anno[rownames(W_meth), "UCSC_RefGene_Group"],
  stringsAsFactors = FALSE
)
mp <- mp[!is.na(mp$island) & mp$island == "Island" &
           !is.na(mp$group) & is_promoter(mp$group) & mp$gene != "", ]

pairs <- do.call(rbind, lapply(seq_len(nrow(mp)), function(i) {
  genes <- unique(strsplit(mp$gene[i], ";")[[1]])
  genes <- genes[genes %in% names(expr_w)]
  if (length(genes) == 0) return(NULL)
  data.frame(probe = mp$probe[i], gene = genes,
             meth_weight = mp$weight[i], expr_weight = expr_w[genes],
             stringsAsFactors = FALSE)
}))

write.csv(pairs, "results/f_shared_coupling_pairs.csv", row.names = FALSE)

shared_factors <- rownames(r2)[shared_mask]

coupling <- do.call(rbind, lapply(shared_factors, function(f) {
  ew <- W_expr[, f]
  mw <- W_meth[mp$probe, f]
  pr <- do.call(rbind, lapply(seq_len(nrow(mp)), function(i) {
    g <- unique(strsplit(mp$gene[i], ";")[[1]])
    g <- g[g %in% names(ew)]
    if (length(g) == 0) return(NULL)
    data.frame(meth_weight = mw[i], expr_weight = ew[g])
  }))
  ct <- cor.test(pr$meth_weight, pr$expr_weight)
  data.frame(factor = f, n_pairs = nrow(pr),
             expr_r2 = r2[f, "expression"], meth_r2 = r2[f, "methylation"],
             cor = round(unname(ct$estimate), 3), p = signif(ct$p.value, 3))
}))
write.csv(coupling, "results/f_shared_coupling_test.csv", row.names = FALSE)

png("results/figures/f_shared_coupling_scatter.png", width = 1400, height = 1200, res = 150)
plot(pairs$meth_weight, pairs$expr_weight, pch = 16, col = "#00000055",
     xlab = "Factor 2 methylation weight (promoter CpG island)",
     ylab = "Factor 2 expression weight")
abline(lm(expr_weight ~ meth_weight, data = pairs), col = "firebrick", lwd = 2)
abline(h = 0, v = 0, col = "grey70", lty = 2)
dev.off()

message("shared factor: ", shared_factor, "   meth-private factor: ", private_factor)
message("promoter-island probes matched to expressed genes: ", nrow(pairs))
message("weight coupling  r = ", coupling$cor, "  p = ", coupling$p_value)
print(coupling)
