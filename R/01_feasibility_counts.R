library(TCGAbiolinks)
library(dplyr)
library(tidyr)
library(readr)

projects <- c("TARGET-WT", "TARGET-ALL-P2", "TARGET-ALL-P3", "TARGET-OS")
gdc_dir <- "GDCdata"
dir.create(gdc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("results", recursive = TRUE, showWarnings = FALSE)

sample_type_map <- c(
  "01" = "primary", "03" = "primary", "09" = "primary",
  "02" = "recurrent", "04" = "recurrent", "40" = "recurrent",
  "05" = "primary", "06" = "metastatic", "07" = "metastatic",
  "10" = "normal", "11" = "normal", "13" = "normal", "14" = "normal",
  "20" = "other", "41" = "recurrent", "50" = "other"
)

fetch_maf <- function(proj) {
  q <- GDCquery(
    project = proj,
    data.category = "Simple Nucleotide Variation",
    data.type = "Masked Somatic Mutation",
    access = "open"
  )
  GDCdownload(q, directory = gdc_dir)
  GDCprepare(q, directory = gdc_dir)
}

annotate <- function(maf, proj) {
  maf %>%
    mutate(
      barcode = as.character(Tumor_Sample_Barcode),
      sample_id = substr(barcode, 1, 19),
      patient_id = substr(barcode, 1, 16),
      type_code = substr(barcode, 18, 19),
      sample_class = dplyr::recode(type_code, !!!as.list(sample_type_map), .default = "unknown"),
      project = proj
    )
}

all_maf <- list()
for (p in projects) {
  message("=== ", p, " ===")
  res <- tryCatch(fetch_maf(p), error = function(e) {
    message("  failed: ", conditionMessage(e)); NULL
  })
  if (!is.null(res)) all_maf[[p]] <- annotate(res, p)
}

maf_all <- bind_rows(all_maf)

snv <- maf_all %>% filter(Variant_Type == "SNP")

per_sample <- snv %>%
  count(project, patient_id, sample_id, sample_class, name = "n_snv")

burden <- per_sample %>%
  group_by(project, sample_class) %>%
  summarise(
    n_samples = n(),
    median_snv = median(n_snv),
    q25 = quantile(n_snv, 0.25),
    q75 = quantile(n_snv, 0.75),
    min_snv = min(n_snv),
    max_snv = max(n_snv),
    frac_over_50 = mean(n_snv >= 50),
    frac_over_100 = mean(n_snv >= 100),
    .groups = "drop"
  ) %>%
  arrange(project, sample_class)

pairing <- per_sample %>%
  filter(sample_class %in% c("primary", "recurrent")) %>%
  distinct(project, patient_id, sample_class) %>%
  mutate(present = TRUE) %>%
  pivot_wider(names_from = sample_class, values_from = present, values_fill = FALSE) %>%
  group_by(project) %>%
  summarise(
    n_patients = n(),
    primary_only = sum(primary & !recurrent),
    recurrent_only = sum(!primary & recurrent),
    paired = sum(primary & recurrent),
    .groups = "drop"
  )

paired_burden <- per_sample %>%
  semi_join(
    per_sample %>%
      filter(sample_class %in% c("primary", "recurrent")) %>%
      distinct(project, patient_id, sample_class) %>%
      count(project, patient_id) %>%
      filter(n == 2),
    by = c("project", "patient_id")
  ) %>%
  filter(sample_class %in% c("primary", "recurrent")) %>%
  group_by(project, sample_class) %>%
  summarise(
    n = n(),
    median_snv = median(n_snv),
    frac_over_50 = mean(n_snv >= 50),
    .groups = "drop"
  )

write_csv(per_sample, "results/per_sample_snv_counts.csv")
write_csv(burden, "results/burden_by_project_sampleclass.csv")
write_csv(pairing, "results/primary_recurrent_pairing.csv")
write_csv(paired_burden, "results/paired_patient_burden.csv")

print(as.data.frame(burden))
cat("\n")
print(as.data.frame(pairing))
cat("\n")
print(as.data.frame(paired_burden))
cat("\n")
print(sessionInfo())
