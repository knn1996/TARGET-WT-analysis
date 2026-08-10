library(TCGAbiolinks)
library(dplyr)

mofa_input <- readRDS("data/mofa_input.rds")
mofa_cases <- colnames(mofa_input$expression)

sum(mofa_cases %in% clin[[2]])
length(grep("^TARGET-50-C", clin[[2]]))

clin <- GDCquery_clinic("TARGET-WT", type = "clinical")
clin <- clin[, !duplicated(names(clin))]

clin_aligned <- data.frame(case = mofa_cases) %>%
  left_join(clin, by = c("case" = "submitter_id"))

message("cases matched to clinical: ", sum(!is.na(clin_aligned$vital_status)),
        " of ", length(mofa_cases))

clin_aligned <- clin_aligned %>%
  mutate(
    os_time = as.numeric(ifelse(vital_status == "Dead",
                                as.numeric(days_to_death),
                                as.numeric(days_to_last_follow_up))),
    os_event = as.integer(vital_status == "Dead"),
    age_years = as.numeric(age_at_diagnosis) / 365.25
  )

completeness <- data.frame(
  field = names(clin_aligned),
  n_missing = colSums(is.na(clin_aligned)),
  n_unique = sapply(clin_aligned, function(x) length(unique(x[!is.na(x)])))
)
completeness <- completeness[order(completeness$n_missing), ]
write.csv(completeness, "results/clinical_field_completeness.csv", row.names = FALSE)
print(completeness)

saveRDS(clin_aligned, "data/clinical_aligned.rds")

q_sup <- GDCquery(
  project = "TARGET-WT",
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  data.format = "BCR Biotab"
)
print(getResults(q_sup)[, c("file_name", "data_type")])
