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

# Any NA rows or duplicates for primary key 'acc'?
is.na(df$acc) %>% table() # None
duplicated(df$acc) %>% table() # None

# Check file sizes
is.na(df$mbytes) %>% table()  # None
max(df$mbytes)  # 39.7 Gb
min(df$mbytes)  # 0 Mb

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

df_wide <- as_tibble(cbind(df[, !names(df) %in% c("attributes", "jattr")], jattr_tibble), .name_repair = "unique")


#### ---- Data analysis ----
# Factors affecting SRA run processing order:
#   Foraging bees = female, worker caste
#   Whole body samples
#   Head only samples
#   Thorax only samples
#   Abdomen only samples
#   Brief summary of above 4 categories before dividing further.


glimpse(df_wide)


# Table occurrences of keywords by column.

# Define keywords and exclusions. 
# Note: bounded words \\b;
#       plural forms included;
#       "whole.?body" is used so regex searches for any character in between the two words 0 or 1 time.
keywords <- c("abdomen", "drone", "female", "forager", "head", "male", "nurse", "queen", "thorax", "whole body", "worker")
keywords_re <- keywords
keywords_re[-c(3, 6, 10)] <- paste0("\\b", keywords[-c(3, 6, 10)], "(s)?\\b")
keywords_re[c(3, 6)] <- paste0("\\b", keywords[c(3,6)], "\\b")
keywords_re[c(9,10)] <- c("\\bthora[xc](es)?\\b", "\\bwhole.?bod[yi](es)?\\b")

exclusions_re <- 
  c("embryo(s)?", "larva(e)?", "pupa(e)?") %>% 
  paste0("\\b", ., "\\b")

# Exclude rows from the expanded column
rows_to_exclude <- list()
for (exclude in exclusions_re) {
  rows_to_exclude[[exclude]] <-
    df_wide %>% 
    select(where(is.character)) %>% 
    sapply(., function(x) { grep(exclude, x, ignore.case = T)}) %>%
    unlist %>% 
    as.vector %>% 
    unique
}
rows_to_exclude <- 
  rows_to_exclude %>% 
  unlist %>% 
  unique
df_wide_ex <- df_wide[-rows_to_exclude, ]

# Grep each row of the table for the keyword, producing a list of all columns
# with the row number under the column the keyword appeared in.
# Count the number of row indices for each column (length).
keyword_list <- list()
for (keyword in keywords_re) {
  keyword_list[[keyword]] <- 
    df_wide_ex %>% 
    select(where(is.character)) %>%
    sapply(., function(x) { grep(keyword, x, ignore.case = T)}) %>%
    lapply(., length) %>% 
    unlist
}

# Tidy.
keyword_find <- as.data.frame(do.call(rbind, keyword_list))
keyword_find <- keyword_find[ , colSums(keyword_find) != 0]
ncol(keyword_find) # the keywords appear in this many columns of the expanded metadata.

# Summary of keyword spread.
data.frame(freq = sort(colSums(keyword_find), decreasing = T)) # most frequent columns

###
# Exclude matches with >60 chars

# Holds the row numbers where a match is found
keyword_list_filter <- list()
for (keyword in keywords_re) {
  keyword_list_filter[[keyword]] <- 
    df_wide_ex %>% 
    select(where(is.character)) %>%
    sapply(., function(x) { grep(keyword, x, ignore.case = T)})
}
names(keyword_list_filter) <- keywords

# Function to test if the number of characters in the matched cell is <60 chars, returns row number or zero if >=60.
return_rows <- function(row_ref, field_name) {
  ifelse(nchar(df_wide_ex[row_ref, field_name]) < 60, row_ref, 0)
}

# Put the row references in a new list.
new_list <- list()
for(key in names(keyword_list_filter)) {
  for(field in names(keyword_list_filter[[key]])) {
    new_list[[key]][[field]] <- 
      if(length(keyword_list_filter[[key]][[field]]) > 0) {
        sapply(keyword_list_filter[[key]][[field]], return_rows, field)
      }
  }
}

# Filter out lists of all zeros and zeros within lists.
filtered_list <- 
  new_list %>% 
  lapply(., function(x) {
    lapply(x, function(y) {
      y[sapply(y, function(z) {z != 0})]
    })
  }) %>% 
  lapply(., function(x) {
    Filter(function(y) sum(y) != 0, x)
  })

# Summary of more accurate keyword spread.
filtered_list_counts <- 
  filtered_list %>%
  lapply(., function(x) {
    lapply(x, length)
  }) %>%
  unlist %>% 
  data.frame(freq = .)

# Summary of more accurate keyword spread.
row.names(filtered_list_counts) %>% 
  gsub(".*\\.", "", .) %>% 
  tibble(field = ., freq = filtered_list_counts$freq) %>% 
  group_by(field) %>% 
  summarise(sum = sum(freq)) %>% 
  arrange(desc(sum))
###
# Use the same process for keyword combinations
# first replace cells >100 char with NA in df_wide_ex?
# df_wide_ex2 <-
#   df_wide_ex %>% 
#   mutate(across(where(is.character), sapply, function(x){
#     if_else(nchar(x) > 100, NA, x, missing = NA)
#   }))
df_wide_ex2 <-
  df_wide_ex %>% 
  mutate(across(where(is.character), function(x){
    sapply(x, function(y) {
      if_else(nchar(y) > 100, NA, y, missing = NA)
    })
  }))

# All combinations (x2) of keywords with replacement and distinct items.
keywords_x2 <- arrangements::combinations(x = keywords_re, k = 2, replace = TRUE)


# Filter for both keywords, ignore case.
df_wide_in <- apply(keywords_x2, 1, function(x) {
  df_wide_ex %>% 
    select(where(is.character)) %>% 
    filter(if_any(everything(), ~ grepl(x[1], ., ignore.case = T)) & if_any(everything(), ~ grepl(x[2], ., ignore.case = T)))
})
df_wide_in_count <- lapply(df_wide_in, count)


# Combine keyword combinations and their counts.
keywords_x2_count <-
  as.data.frame(cbind(
    arrangements::combinations(x = keywords, k = 2, replace = TRUE),
    unlist(df_wide_in_count)))

###

# df_wide_in contains only the chr columns, need 'mbytes' column from df_wide
df_wide_sizes <-
  df_wide_in %>% 
  lapply(., select, acc) %>% 
  lapply(., function(x) {
    df_wide %>% 
      semi_join(x, by = "acc") %>% 
      summarise(sum = sum(mbytes)) %>% 
      pull(sum)
  })


# Add gb info to key_key_count table. 
keywords_x2_count$gb <- round(unlist(df_wide_sizes)/1000, digits = 0)

# Duplicate the counts for the reciprocal combinations to make the plot more readable.
recip_key_x2_count <- 
  rename(keywords_x2_count, V1 = V2, V2 = V1) %>% 
  bind_rows(., keywords_x2_count) %>% 
  distinct
  
recip_key_x2_count$V3 <- as.numeric(recip_key_x2_count$V3)
recip_key_x2_count$gb <- as.numeric(recip_key_x2_count$gb)


ggplot(recip_key_x2_count, aes(V1, V2)) +
  geom_tile(aes(fill = V3)) +
  geom_text(aes(label = V3)) +
  scale_fill_gradient(low = "white",
                      high = "orange",
                      name = "SRA runs") +
  labs(title = "Number of SRA runs given a combination of two keywords",
       subtitle = "Diagonal represents a single keyword") +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4, size = 12),
        axis.text.y = element_text(size = 12),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        panel.grid.major = element_blank(),
        plot.caption.position = "plot",
        plot.caption = element_text(hjust = 0),
        plot.title.position = "plot",
        panel.background = element_blank(),
        panel.border = element_rect(color = "grey60", fill = NA),
        legend.position = "right",
        strip.background = element_rect(color = "grey60", fill = "white"))


  
ggplot(recip_key_x2_count, aes(V1, V2)) +
  geom_tile(aes(fill = gb)) +
  geom_text(aes(label = gb)) +
  scale_fill_gradient(low = "white",
                      high = "lightblue",
                      name = "Size (Gb)") +
  labs(title = "Total size of SRA runs given a combination of two keywords",
       subtitle = "Diagonal represents a single keyword") +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4, size = 12),
        axis.text.y = element_text(size = 12),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        panel.grid.major = element_blank(),
        plot.caption.position = "plot",
        plot.caption = element_text(hjust = 0),
        plot.title.position = "plot",
        panel.background = element_blank(),
        panel.border = element_rect(color = "grey60", fill = NA),
        legend.position = "right",
        strip.background = element_rect(color = "grey60", fill = "white"))




# consider
# female: { queen, worker: { nurse, forager } }
# male:   { drone }
# tissue: { whole.body: { head, thorax, abdomen } }

# First subset worker, whole.body. 
# Seems to be a huge number where the caste is not defined: 721

### Function to calculate total mbytes of a subset
calc_mbytes <- function(dataframe, ...) {
  keywords <- list(...)
  keywords <- unlist(keywords)
  i <- dataframe
  # This is an AND loop for the keywords.
  for (keyword in keywords) {
    i <- 
      i %>% 
      select(where(is.character)) %>% 
      filter(if_any(everything(), ~ grepl(paste0("\\b", keyword, "\\b"), ., ignore.case = T)))
  }
  # Get the accessions left over.
  acc <-
    i %>% 
    select(acc)
  # Filter the dataframe arg with the acc left over and sum mbytes.
  n <-
    dataframe %>% 
    semi_join(acc, by = "acc") %>% 
    summarise(sum = sum(mbytes)) %>% 
    pull(sum)
  # Print in Gb.
  print(paste0(n/1000, " Gb")) 
}

# case1 where whole.body
wb_Fworker <-
  df_wide_ex2 %>% 
    select(where(is.character)) %>% 
    filter((if_any(everything(), ~ grepl("whole.body", ., ignore.case = T)))) %>% 
    filter((if_any(everything(), ~ grepl("\\bworker(s)?\\b", ., ignore.case = T))) | (if_any(everything(), ~ grepl("\\bfemale(s)?\\b", ., ignore.case = T))))

# df_wide_in contains only the chr columns, need 'mbytes' column from df_wide
wb_Fworker <-
  df_wide_in %>% 
  lapply(., select, acc) %>% 
  lapply(., function(x) {
    wb_Fworker %>% 
      semi_join(x, by = "acc") %>% 
      summarise(sum = sum(mbytes)) %>% 
      pull(sum)
  })

# wb_male
wb_male <-
  df_wide_ex2 %>% 
  select(where(is.character)) %>% 
  filter((if_any(everything(), ~ grepl("whole.body", ., ignore.case = T)))) %>% 
  filter((if_any(everything(), ~ grepl("\\bmale(s)\\b?", ., ignore.case = T))))
wb <-
  df_wide_ex2 %>% 
  select(where(is.character)) %>% 
  filter((if_any(everything(), ~ grepl("\\bwhole.?bod[yi](es)?\\b", ., ignore.case = T))))

#### ---- RUBBISH ----

z<-df_wide_ex %>% 
  select(where(is.character)) %>% 
  filter(if_any(everything(), ~ grepl("\\bthora[xc](es)?\\b", ., ignore.case = T)) &
           if_any(everything(), ~ grepl("\\bthora[xc](es)?\\b", ., ignore.case = T)))


colnames(mf) <- keywords
row.names(mf) <- 
# The number of columns of the matrix is the number of elements in each combination (k).
num_cols <- length(x)
# The number of rows of the matrix is the total number of possible combinations.
num_rows <- nrow(combinations)
# For simplicity, reshape the matrix to have one column per combination.
resulting_matrix <- matrix(combinations, nrow = num_rows, ncol = num_cols)



y <-
  df_wide %>%
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
  
  
  
  
  





df_wide %>%
  group_by(acc) %>% 
  
  mutate(found_in = pmap_chr(across(everything()),

cols_with_string <- names(y)[sapply(tibble_as_list, function(x) { grepl(".*female.*", x, ignore.case = T) })]

df_wide %>% 
  select(where(is.character)) %>% 
  mutate(female = rowSums(across(all_of(cols), `%in%`, "female")))

df_wide %>% 
  select(where(is.character)) %>% 
  rowwise() %>% 
  mutate(female_count = sum(str_detect(c_across(everything()), fixed("female")))) %>% 
  ungroup() %>% 
  filter(female_count != NA)
  
  


df_wide %>% mutate(d9 = rowSums(across(., `%in%`, "female")))


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

df_wide <- as_tibble(cbind(bioproject = df$bioproject, acc = df$acc, jattr_tibble))

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

df_wide %>% 
  mutate(tissue1 = str_to_lower(tissue_sam_ss_dpl145), .after = acc) %>% 
  filter(str_detect(tissue1, "thorax")) %>% 
  group_by(bioproject, tissue1) %>% 
  count() %>% 
  arrange(desc(n))
  

keyword_list_filter %>% 
  





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
