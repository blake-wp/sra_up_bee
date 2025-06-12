library(tidyverse)
library(bigrquery)
library(dotenv)

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

# Any NA rows for primary key 'acc'?
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
# One bioproject, remove from list.
df <-
  df %>% 
  filter(mbytes != 0)

# First, looks like the 'attributes' and 'jattr' are nested as JSON objects.
# These each contain the same sample info, e.g. 'tissue_sam_ss_dpl145': 'thoracic muscle'.
# Expand 'jattr' and add it to the df. Drop 'jattr' and 'attributes' nested columns.
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


#### ---- Data analysis ----
# Factors affecting SRA run processing order:
#   Foraging bees = female, worker caste
#   Whole body samples
#   Head only samples
#   Thorax only samples
#   Abdomen only samples
#   Brief summary of above 4 categories before dividing further.


glimpse(df)


# Table occurrences of keywords by column.

# Get only the character columns
df_chr <-
  jattr_df %>% 
  select(where(is.character))

# Define keywords.
keywords <- c("female", "male", "drone", "queen", "worker", "nurse", "whole body", "head", "thorax", "abdomen")
keyword_list <- list()
# Get the columns the keywords appear in.








#### ---- RUBBISH ----
for (keyword in keywords) {
  keyword_list[keyword] <- 
      as_tibble(sapply(y, function(x) { grep(keyword, x, ignore.case = T)}) %>%
        lapply(., length) %>%
        unlist)
  
}
x <- as.data.frame(do.call(rbind, keyword_list))
names(x) <- names(y)
x<- x[ , colSums(x) != 0]
x

# intersections
# whole body and female
y %>% 
  filter(if_any(everything(), ~ grepl('female', .)) & if_any(everything(), ~ grepl('female', .))) %>% 
  count

# matrix of keywords
x <- c(1, 2, 3, 4)
k <- 2
combinations <- combn(keywords, k)
mf <- data.frame()
colnames(mf) <- keywords
row.names(mf) <- 
# The number of columns of the matrix is the number of elements in each combination (k).
num_cols <- length(x)
# The number of rows of the matrix is the total number of possible combinations.
num_rows <- nrow(combinations)
# For simplicity, reshape the matrix to have one column per combination.
resulting_matrix <- matrix(combinations, nrow = num_rows, ncol = num_cols)



y <-
  jattr_df %>%
  mutate(found_in = pmap_chr(across(everything()), ~ {
    cols_with_string <- names(.)[sapply(list(...), function(x) grepl("female", x))]
    if (length(cols_with_string) > 0) {
      paste(cols_with_string, collapse = ", ")
    } else {
      NA_character_
    }
  }))



%>% 
  filter(row_number() %in% unlist(sapply(., function(x) { grep("female", x)}))) %>% 
  
  
  
  
  





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
