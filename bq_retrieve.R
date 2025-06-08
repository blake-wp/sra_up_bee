library(tidyverse)
library(bigrquery)
library(dotenv)

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

# Any NA rows?
is.na(df$acc) %>% table() # None

# Check file sizes
is.na(df$mbytes) %>% table()  # None
max(df$mbytes)  # 39.7 Gb
min(df$mbytes)  # 0 Mb

# What entries have 0 Mb?
df %>%
  filter(mbytes == 0) %>% 
  group_by(bioproject) %>% 
  count()
# Three bioprojects, remove from list.
df <-
  df %>% 
  filter(mbytes != 0)

#### ---- Data analysis ----
# Factors affecting SRA run processing order:
#   Foraging bees = female, worker caste
#   Whole body samples
#   Head only samples
#   Thorax only samples
#   Abdomen only samples
#   Brief summary of above 4 categories before dividing further.

# First expand the 'jattr' nested table and add to the df
library(jsonlite)
parsed_jattr <- lapply(df$jattr, fromJSON)
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

jattr_df <- as_tibble(cbind(df[, !names(df) %in% c("attributes", "jattr")], jattr_tibble), .name_repair = "unique")

# Find instances of 'female' across expanded df.

y <-
  jattr_df %>%
  mutate(found_in = pmap_chr(across(everything()), ~ {
    cols_with_string <- names(.)[sapply(list(...), function(x) grepl(".*female.*", x))]
    if (length(cols_with_string) > 0) {
      paste(cols_with_string, collapse = ", ")
    } else {
      NA_character_
    }
  }))


y <-
  jattr_df %>% 
  select(where(is.character))
%>% 
  filter(row_number() %in% unlist(sapply(., function(x) { grep("female", x)}))) %>% 
  
  
  
  
  



z <- y[unlist(sapply(y, function(x) { grep("female", x)})), ]

jattr_df %>%
  group_by(acc) %>% 
  
  mutate(found_in = pmap_chr(across(everything()),

cols_with_string <- names(y)[sapply(tibble_as_list, function(x) { grepl(".*female.*", x, ignore.case = T) })]

jattr_df %>% 
  select(where(is.character)) %>% 
  mutate(female = rowSums(across(all_of(cols), `%in%`, "female")))

jattr_df %>% 
  select(where(is.character)) %>% 
  rowwise() %>% 
  mutate(female_count = sum(str_detect(c_across(everything()), fixed("female")))) %>% 
  ungroup() %>% 
  filter(female_count != NA)
  
  


jattr_df %>% mutate(d9 = rowSums(across(., `%in%`, "female")))


# Attributes column is in JSON format and contains the sample descriptions.
library(jsonlite)
parsed_jattr <- lapply(df$jattr, fromJSON)

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

jattr_df <- as_tibble(cbind(bioproject = df$bioproject, acc = df$acc, jattr_tibble))

jattr_tibble %>%
  summarise_all(~ sum(!is.na(.))) %>% 
  gather(column, value) %>% 
  arrange(desc(value)) %>% 
  print(n = 20)


as.data.frame(
  sapply(jattr_tibble, function(x) {
  sum(!is.na(x))
})
)

jattr_df %>% 
  mutate(tissue1 = str_to_lower(tissue_sam_ss_dpl145), .after = acc) %>% 
  filter(str_detect(tissue1, "thorax")) %>% 
  group_by(bioproject, tissue1) %>% 
  count() %>% 
  arrange(desc(n))
  

    





tissue_sam <-
  parsed_attributes %>% 
  do.call(rbind, .) %>% 
  do.call(rbind, .) %>%
  as_tibble() %>%
  filter(k == 'tissue_sam_ss_dpl145') %>% 
  mutate(v = str_to_lower(v)) %>% 
  count(v) %>% 
  arrange(desc(n))

isolate_sam <-
  parsed_attributes %>% 
  do.call(rbind, .) %>% 
  do.call(rbind, .) %>%
  as_tibble() %>%
  filter(k == 'isolate_sam_ss_dpl100') %>% 
  mutate(v = str_to_lower(v)) %>% 
  count(v) %>% 
  arrange(desc(n))

# Create df with 'acc' and 'tissue_sam_ss_dpl145' coumns
