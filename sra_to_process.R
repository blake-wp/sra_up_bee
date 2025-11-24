library(tidyverse)
library(bigrquery)
library(dotenv)
library(jsonlite)

rm(list = ls())
#### ---- Data import ----
# Use environment variables to mask private data.
load_dot_env()

# Read the csv bq table if available, otherwise redo the SQL query.
df <- try(read.csv("bquxjob_722b10fc_1971e82ea0f.csv"))
if (inherits(df, "try-error")) {
  # BigQeury SQL query
  query <-
    "
    SELECT *
    FROM `nih-sra-datastore.sra.metadata`
    WHERE
      organism = 'Apis mellifera' AND
      consent = 'public' AND
      assay_type = 'RNA-Seq'
    "
  tb <- bq_project_query(Sys.getenv("BIGQUERY_PROJECT"), query)
  df <- bq_table_download(tb)
}


#### ---- Data cleaning ----

glimpse(df)

# Any NA rows or duplicates for primary key 'acc'?
is.na(df$acc) %>% table() # None
duplicated(df$acc) %>% table() # None

# Check file sizes
is.na(df$mbytes) %>% table() # None
max(df$mbytes) # 39.7 Gb
min(df$mbytes) # 0 Mb

# What entries have 0 Mb?
df %>%
  filter(mbytes == 0) %>%
  group_by(bioproject) %>%
  count()

df %>%
  filter(bioproject == "PRJNA477521") %>%
  select(acc, mbases, mbytes, library_name) %>%
  glimpse

# One bioproject with samples with blank reads? Remove from list.
df <-
  df %>%
  filter(mbytes != 0)

nrow(df)

# Keep Illumina runs only
df <-
  df |> filter(platform == "ILLUMINA")

# Plot histo of run sizes
ggplot(df, aes(mbytes)) +
  geom_histogram(binwidth = 200) +
  ggtitle("Distribution of SRA run file size in the dataset (200 MB bins)") +
  theme_minimal()

# how many <5 Gb?
smaller_files <-
  df |> filter(mbytes < 5000)
smaller_files |> count() # 4069

ggplot(smaller_files, aes(mbytes)) +
  geom_histogram(binwidth = 200) +
  ggtitle("Distribution of SRA run file size in the dataset (200 MB bins)") +
  theme_minimal()

# Split into groups to better manage

for (i in 0:9) {
  assign(
    paste0("bin", i + 1),
    df |> filter(mbytes < (1001 + (1000 * i)) & mbytes > (1000 * i))
  )
  write_lines(
    get(paste0("bin", i + 1))[[1]],
    paste0("srr_", (1000 * i) + 1, "-", (1000 + (1000 * i)), "_MB.txt")
  )
}
