library(tidyverse)
library(bigrquery)

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
glimpse(df)

count_attributes <-
  df %>% 
  select(attributes) %>%
  do.call(rbind, .) %>%
  do.call(rbind, .) %>%
  count(k) %>% 
  arrange(desc(n))

tissue_sam <- 
  df %>% 
  select(attributes) %>%
  unnest(cols = c(attributes)) %>%
  filter(k == 'tissue_sam_ss_dpl145') %>% 
  count(v) %>% 
  arrange(desc(n))

# for i in length of df obs. extract nested entry
df$attributes[[i]][df$attributes[[i]]$k == 'tissue_sam_ss_dpl145', 2]
