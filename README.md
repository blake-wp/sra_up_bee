Find bee SRA runs from metadata
================
Blake Paget
2025

## Background

RNA-Seq data will be analysed to determine the expression pattern of an
uncharacterised protein. The Sequence Read Archive (SRA) maintained by
NCBI will be used. The metadata will be searched for appropriate runs.
In a separate publication selected runs will be downloaded and
analysed.<br>

“SRA has deposited its metadata into BigQuery to provide the
bioinformatics community with programmatic access to this data. You can
now search across the entire SRA by sequencing methodologies and sample
attributes.”<br>

A simple SQL query will be used to retrieve a table for further analysis
in R. The sample attributes are an inconsistent mess and rather than
trying to sort that with complex SQL queries, I’ll just use R.

## Script

#### Load libraries

``` r
library(tidyverse)
library(bigrquery)
library(dotenv)
```

<br>

#### Import csv from BigQuery

``` r
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
```

<br>

#### Data cleaning

``` r
glimpse(df)
```

    ## Rows: 4,692
    ## Columns: 37
    ## $ acc                                 <chr> "SRR7770436", "SRR22681164", "SRR6…
    ## $ assay_type                          <chr> "RNA-Seq", "RNA-Seq", "RNA-Seq", "…
    ## $ center_name                         <chr> "GEO", "ZHEJIANG UNIVERSITY OF TEC…
    ## $ consent                             <chr> "public", "public", "public", "pub…
    ## $ experiment                          <chr> "SRX4626009", "SRX18643300", "SRX3…
    ## $ sample_name                         <chr> "GSM3362362", "R-TRZ-2", "RNA109_P…
    ## $ instrument                          <chr> "Illumina HiSeq 4000", "Illumina N…
    ## $ librarylayout                       <chr> "PAIRED", "PAIRED", "PAIRED", "SIN…
    ## $ libraryselection                    <chr> "cDNA", "Oligo-dT", "RANDOM", "cDN…
    ## $ librarysource                       <chr> "TRANSCRIPTOMIC", "TRANSCRIPTOMIC"…
    ## $ platform                            <chr> "ILLUMINA", "ILLUMINA", "ILLUMINA"…
    ## $ sample_acc                          <chr> "SRS3727542", "SRS16088054", "SRS2…
    ## $ biosample                           <chr> "SAMN09937971", "SAMN32035609", "S…
    ## $ organism                            <chr> "Apis mellifera", "Apis mellifera"…
    ## $ sra_study                           <chr> "SRP159176", "SRP412451", "SRP1307…
    ## $ releasedate                         <chr> "2019-10-31 00:00:00.000000 UTC", …
    ## $ bioproject                          <chr> "PRJNA488619", "PRJNA908448", "PRJ…
    ## $ mbytes                              <int> 6041, 2127, 2302, 336, 2082, 3246,…
    ## $ loaddate                            <lgl> NA, NA, NA, NA, NA, NA, NA, NA, NA…
    ## $ avgspotlen                          <int> 300, 300, 302, 100, 302, 302, 150,…
    ## $ mbases                              <int> 15048, 6976, 6517, 1081, 4952, 100…
    ## $ insertsize                          <lgl> NA, NA, NA, NA, NA, NA, NA, NA, NA…
    ## $ library_name                        <chr> "", "R-TRZ-2", "112_S9_L001_R1_001…
    ## $ biosamplemodel_sam                  <chr> "[]", "[Model organism or animal]"…
    ## $ collection_date_sam                 <chr> "", "", "2015-08-10", "", "2013-07…
    ## $ geo_loc_name_country_calc           <chr> "", "", "USA", "", "USA", "", "", …
    ## $ geo_loc_name_country_continent_calc <chr> "", "", "North America", "", "Nort…
    ## $ geo_loc_name_sam                    <chr> "[]", "[]", "[USA: NC]", "[]", "[U…
    ## $ ena_first_public_run                <chr> "[]", "[]", "[]", "[]", "[]", "[]"…
    ## $ ena_last_update_run                 <chr> "[]", "[]", "[]", "[]", "[]", "[]"…
    ## $ sample_name_sam                     <chr> "[]", "[]", "[]", "[]", "[]", "[]"…
    ## $ datastore_filetype                  <chr> "[fastq,run.zq,sra]", "[fastq,run.…
    ## $ datastore_provider                  <chr> "[gs,ncbi,s3]", "[gs,ncbi,s3]", "[…
    ## $ datastore_region                    <chr> "[gs.us-east1,ncbi.public,s3.us-ea…
    ## $ attributes                          <chr> "{\n  \"attributes\": [{\n    \"k\…
    ## $ run_file_version                    <int> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1…
    ## $ jattr                               <chr> "{\"geo_accession_exp\": [\"GSM336…

``` r
# Any NA rows for primary key 'acc'?
is.na(df$acc) %>% table()
```

    ## .
    ## FALSE 
    ##  4692

``` r
# Check file sizes
is.na(df$mbytes) %>% table()
```

    ## .
    ## FALSE 
    ##  4692

``` r
max(df$mbytes)
```

    ## [1] 39732

``` r
min(df$mbytes)
```

    ## [1] 0

``` r
# What entries have 0 Mb?
df %>%
  filter(mbytes == 0) %>% 
  group_by(bioproject) %>% 
  count()
```

    ## # A tibble: 1 × 2
    ## # Groups:   bioproject [1]
    ##   bioproject      n
    ##   <chr>       <int>
    ## 1 PRJNA477521    30

``` r
# One bioproject, remove from list.
df <-
  df %>% 
  filter(mbytes != 0)

# Looks like the 'attributes' and 'jattr' are nested as JSON objects. 
# Not shown here, I have extracted these and compared them.
# These each contain the same sample info, e.g. 'tissue_sam_ss_dpl145': 'thoracic muscle'.
# Expand 'jattr' and add it to the df. Drop 'jattr' and 'attributes' nested columns.
```

``` r
library(jsonlite)
```

``` r
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
```

    ## New names:
    ## • `collection_date_sam` -> `collection_date_sam...25`
    ## • `collection_date_sam` -> `collection_date_sam...49`

``` r
nrow(jattr_df)
```

    ## [1] 4662

<br>

#### Data analysis

``` r
# Get only the character columns
df_chr <-
  jattr_df %>% 
  select(where(is.character))
```
