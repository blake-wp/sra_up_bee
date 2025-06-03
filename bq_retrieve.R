library(tidyverse)
library(bigrquery)

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
  tb <- bq_project_query("mimetic-planet-458418-e4", query)
  df <- bq_table_download(tb)
}

glimpse(df)

max(df$mbytes)
min(df$mbytes)

# Attributes column is in JSON format and contains the sample descriptions.
library(jsonlite)
parsed_attributes <- lapply(df$attributes, fromJSON)
parsed_jattr <- lapply(df$jattr, fromJSON)

attribute_summary <- parsed_attributes %>%
  do.call(rbind, .) %>% 
  do.call(rbind, .) %>%
  as_tibble() %>% 
  count(k) %>% 
  arrange(desc(n))

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

jattr_df %>% 
  mutate(tissue1 = str_to_lower(tissue_sam_ss_dpl145), .after = acc) %>% 
  filter(str_detect(tissue1, "thorax")) %>% 
  group_by(bioproject, tissue1) %>% 
  count() %>% 
  arrange(desc(n))
  

jattr_df %>% 
  mutate(tissue1 = str_to_lower(tissue_sam_ss_dpl145), .after = acc) %>% 
  filter(tissue1 == "thorax") %>% 
  group_by(bioproject) %>% 
  count()
    


parsed_jattr %>%
  lapply(function(x) { tibble(Name = names(x), Value = as.character(unname(x))) } ) %>%
  do.call(rbind, .) %>% 
  count(Name) %>% 
  arrange(desc(n))

# Above two parsed JSON arrays are identical except for 'primary_search' in attributes.

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
