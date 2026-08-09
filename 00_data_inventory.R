library(TCGAbiolinks)
library(dplyr)
library(readr)

project <- "TARGET-WT"
dir.create("results", recursive = TRUE, showWarnings = FALSE)

q_snv <- GDCquery(
  project = project,
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  access = "open"
)

res_snv <- getResults(q_snv)
count(res_snv, experimental_strategy, analysis_workflow_type)

res_snv %>%
  mutate(stype = substr(sapply(strsplit(cases, "-"), `[`, 4), 1, 2)) %>%
  count(stype)

q_expr <- GDCquery(
  project = project,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

q_meth <- GDCquery(
  project = project,
  data.category = "DNA Methylation",
  data.type = "Methylation Beta Value",
  platform = "Illumina Human Methylation 450"
)

queries <- list(snv = q_snv, expr = q_expr, meth = q_meth)

summarise_layer <- function(q, label) {
  res <- getResults(q)
  aliquots <- unique(as.character(res$cases))
  cases <- unique(substr(aliquots, 1, 16))
  size_gb <- if ("file_size" %in% names(res)) sum(as.numeric(res$file_size)) / 1e9 else NA_real_
  list(
    label = label,
    res = res,
    cases = cases,
    row = tibble(
      layer = label,
      n_files = nrow(res),
      n_aliquots = length(aliquots),
      n_cases = length(cases),
      total_size_gb = round(size_gb, 2)
    )
  )
}

layers <- Map(summarise_layer, queries, names(queries))

overview <- bind_rows(lapply(layers, function(x) x$row))

case_sets <- lapply(layers, function(x) x$cases)

pairwise <- tibble(
  comparison = c("snv_expr", "snv_meth", "expr_meth"),
  n_shared = c(
    length(intersect(case_sets$snv, case_sets$expr)),
    length(intersect(case_sets$snv, case_sets$meth)),
    length(intersect(case_sets$expr, case_sets$meth))
  )
)

three_way <- Reduce(intersect, case_sets)

membership <- tibble(
  case_id = sort(unique(unlist(case_sets)))
) %>%
  mutate(
    snv = case_id %in% case_sets$snv,
    expr = case_id %in% case_sets$expr,
    meth = case_id %in% case_sets$meth,
    n_layers = snv + expr + meth
  )

layer_profile <- membership %>%
  count(snv, expr, meth, name = "n_cases") %>%
  arrange(desc(n_cases))

sample_types <- bind_rows(lapply(layers, function(x) {
  tibble(
    layer = x$label,
    type_code = substr(as.character(x$res$cases), 18, 19)
  )
})) %>%
  count(layer, type_code, name = "n_aliquots") %>%
  arrange(layer, desc(n_aliquots))

write_csv(overview, "results/inventory_overview.csv")
write_csv(pairwise, "results/inventory_pairwise.csv")
write_csv(membership, "results/inventory_case_membership.csv")
write_csv(layer_profile, "results/inventory_layer_profile.csv")
write_csv(sample_types, "results/inventory_sample_types.csv")
write_lines(three_way, "results/cases_three_layer.txt")

print(as.data.frame(overview))
cat("\n")
print(as.data.frame(pairwise))
cat("\n")
cat("three-layer intersection:", length(three_way), "cases\n\n")
print(as.data.frame(layer_profile))
cat("\n")
print(as.data.frame(sample_types))
cat("\n")
print(sessionInfo())
