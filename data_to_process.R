library(tidyverse)
library(jsonlite)

rm(list = ls())

path = "H:/counts/"

summary_files <- list.files(
  path = path,
  pattern = "*.txt.summary",
  full.names = T
)
summary_list <- lapply(summary_files, read_delim, delim = '\t')
summary_combined <- reduce(summary_list, left_join, by = "Status")

#rm(summary_list, summary_files)
summary_norm <-
  summary_combined |>
  mutate(across(!any_of("Status"), ~ . / sum(.))) |>
  rowwise() |>
  filter(sum(c_across(-Status)) != 0) |>
  ungroup() |>
  pivot_longer(-Status, names_to = "SRR_runs", values_to = "Value") |>
  pivot_wider(names_from = Status, values_from = Value) |>
  pivot_longer(
    cols = -c(SRR_runs, Assigned),
    names_to = "Unassigned_type",
    values_to = "Unassigned"
  )

ggplot(summary_norm) +
  geom_point(aes(x = Assigned, y = Unassigned)) +
  facet_wrap(~Unassigned_type)

counts_files <- list.files(path = path, pattern = "*.txt$", full.names = T)
counts_list <- lapply(counts_files, read_delim, delim = '\t', skip = 1)
counts_combined <- reduce(counts_list, left_join, by = "Geneid")

counts_scaled_100 <-
  counts_combined |>
  mutate(across(where(is.numeric), ~ .x / max(.x, na.rm = TRUE) * 100))

x <- counts_scaled_100 |> filter(Geneid == "LOC408608")
x <- x |>
  pivot_longer(cols = -Geneid, names_to = "SRR_run", values_to = "values")

ggplot(x) +
  geom_histogram(aes(x = values))

roi <- x |> filter(values > 10)
roi_runs <- roi$SRR_run |> str_replace_all("\\.SAM", "")

roi_runs_meta <-
  data_to_parse |>
  filter(acc %in% roi_runs)

# Extract keywords
parsed_jattr <- lapply(roi_runs_meta$jattr, fromJSON)
jattr_tibble <-
  lapply(parsed_jattr, function(x) {
    t(cbind(names(x), as.character(unname(x))))
  }) %>%
  lapply(function(x) {
    z <- as.data.frame(x[-1, ])
    row.names(z) <- x[1, ]
    as.data.frame(t(z))
  }) %>%
  do.call(plyr::rbind.fill, .) %>%
  as_tibble

df_wide <- as_tibble(
  cbind(
    roi_runs_meta[, !names(roi_runs_meta) %in% c("attributes", "jattr")],
    jattr_tibble
  ),
  .name_repair = "unique"
)

keywords <- c(
  "abdomen",
  "drone",
  "embryo",
  "female",
  "forager",
  "head",
  "larva",
  "male",
  "nurse",
  "pupa",
  "queen",
  "thorax",
  "whole body",
  "worker"
)

keywords_re <- keywords
keywords_re[-c(4, 7, 8, 10, 12, 13)] <- paste0(
  "\\b",
  keywords[-c(4, 7, 8, 10, 12, 13)],
  "(s)?\\b"
)
keywords_re[c(4, 8)] <- paste0("\\b", keywords[c(4, 8)], "\\b")
keywords_re[c(7, 10)] <- paste0("\\b", keywords[c(7, 10)], "(e)?\\b")
keywords_re[c(12, 13)] <- c(
  "\\bthora[xc](es)?\\b",
  "\\bwhole.?bod[yi](es)?\\b"
)


df_with_counts <- 
  roi_runs_meta |> 
  rowwise() %>%
  mutate(
    # Combine all character columns into one string per row
    combined_text = str_to_lower(paste(c_across(where(is.character)), collapse = " ")),
    
    # Find which keywords are present
    matched_keywords = list(keywords_re[str_detect(combined_text, keywords_re)]),
    
    # Count how many matched
    keyword_count = length(matched_keywords)
  ) %>%
  ungroup()

x <- read_csv("roi_runs_meta_with_sample_type.csv")
