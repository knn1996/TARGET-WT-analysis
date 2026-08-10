bc_field <- function(x, i) vapply(strsplit(x, "-"), `[`, character(1), i)

case_id <- function(x) paste(bc_field(x, 1), bc_field(x, 2), bc_field(x, 3), sep = "-")

is_primary_sample <- function(x) substr(bc_field(x, 4), 1, 2) == "01"
