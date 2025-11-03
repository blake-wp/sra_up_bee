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
  geom_histogram(binwidth = 100)

# how many <5 Gb?
smaller_files <-
  df |> filter(mbytes < 5000)
smaller_files |> count()

ggplot(smaller_files, aes(mbytes)) +
  geom_histogram()

# Split into groups to better manage
bin1 <-
  df |> filter(mbytes < 1001)
write_lines(bin1$acc, "srr_0-1000_Gb.txt")
bin2 <-
  df |> filter(mbytes < 2001 & mbytes > 1000)
write_lines(bin2$acc, "srr_1001-2000_Gb.txt")
bin3 <-
  df |> filter(mbytes < 3001 & mbytes > 2000)
write_lines(bin3$acc, "srr_2001-3000_Gb.txt")
