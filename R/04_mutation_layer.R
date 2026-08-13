library(TCGAbiolinks)
library(maftools)
library(dplyr)
source("barcodes.R")

project <- "TARGET-WT"
EXOME_MB <- 38
PREVALENCE_CUTOFF <- 0.10
dir.create("results", showWarnings = FALSE)

q_snv <- GDCquery(
  project = project,
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  access = "open"
)

GDCdownload(q_snv)
maf_df <- GDCprepare(q_snv)

maf_df <- maf_df %>%
  filter(is_primary_sample(Tumor_Sample_Barcode)) %>%
  mutate(case = case_id(Tumor_Sample_Barcode))

message("rows: ", nrow(maf_df), "  cases: ", length(unique(maf_df$case)))
stopifnot(nrow(maf_df) > 0)

maf <- read.maf(maf_df)

coding <- c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
            "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
            "Splice_Site", "Translation_Start_Site", "Nonstop_Mutation")

burden <- maf_df %>%
  filter(Variant_Classification %in% coding) %>%
  dplyr::count(case, name = "n_coding_snv") %>%
  mutate(tmb = n_coding_snv / EXOME_MB) %>%
  arrange(n_coding_snv)

write.csv(burden, "results/wt_mutation_burden.csv", row.names = FALSE)

burden_stats <- data.frame(
  n_cases = nrow(burden),
  median_coding_snv = median(burden$n_coding_snv),
  q1 = quantile(burden$n_coding_snv, 0.25),
  q3 = quantile(burden$n_coding_snv, 0.75),
  min = min(burden$n_coding_snv),
  max = max(burden$n_coding_snv),
  median_tmb = median(burden$tmb)
)
write.csv(burden_stats, "results/wt_mutation_burden_stats.csv", row.names = FALSE)
print(burden_stats)

gene_hits <- maf_df %>%
  filter(Variant_Classification %in% coding) %>%
  distinct(case, Hugo_Symbol)

n_cases <- length(unique(maf_df$case))

gene_freq <- gene_hits %>%
  dplyr::count(Hugo_Symbol, name = "n_mutated") %>%
  mutate(freq = n_mutated / n_cases) %>%
  arrange(desc(n_mutated))

write.csv(gene_freq, "results/wt_gene_frequency.csv", row.names = FALSE)

genes_kept <- gene_freq %>% filter(freq >= PREVALENCE_CUTOFF) %>% pull(Hugo_Symbol)
message("genes surviving >=", PREVALENCE_CUTOFF * 100, "% filter: ", length(genes_kept))

pdf("results/wt_oncoplot.pdf", width = 9, height = 6)
oncoplot(maf, top = 20)
dev.off()

pdf("results/wt_burden_distribution.pdf", width = 6, height = 4)
hist(burden$n_coding_snv, breaks = 20,
     xlab = "coding SNVs per case", main = "TARGET-WT mutation burden (open access, n = 38)")
dev.off()

png("results/figures/wt_oncoplot.png", width = 1800, height = 1200, res = 150)
oncoplot(maf, top = 20)
dev.off()
