library(tidyverse)
library(jsonlite)
library(dotenv)

rm(list = ls())

#### ---- Keyword Extraction Prep ----
load_dot_env()
df <- try(read.csv("bquxjob_722b10fc_1971e82ea0f.csv"))
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
keyword_raw <-
  jattr_tibble |>
  select(social_caste_sam:notes_sam) |>
  mutate(acc = df$acc) |>
  unite(
    "combined",
    social_caste_sam:notes_sam,
    sep = ",",
    remove = TRUE,
    na.rm = TRUE
  )

# Submit this file to Claude and ask to go through the 'combined'
# column and extract significant keywords. Claude should produce a
# text that has keywords grouped by different categories.
write_csv(keyword_raw, "runs_concat_meta.csv")
# Get Claude to create R code that will go through each row of 'df',
# matching keywords, creating a new table with categories as column
# headings and adding found keywords under those columns.

# Convert df to all character for the matching functions
df <-
  df |>
  as_tibble() |>
  mutate(across(everything(), as.character))

#### ---- Code Produced By Claude ----

# Define keyword categories
keywords <- list(
  anatomical_terms = c(
    "abdomen",
    "abdominal",
    "antenna",
    "antennae",
    "antennal",
    "brain",
    "brains",
    "body",
    "bodies",
    "cell",
    "cells",
    "chest",
    "cuticle",
    "exoskeleton",
    "eye",
    "eyes",
    "fatbody",
    "ganglia",
    "ganglion",
    "gland",
    "glands",
    "gonad",
    "gonadal",
    "gut",
    "guts",
    "haemolymph",
    "hemolymph",
    "head",
    "heads",
    "hindgut",
    "hypopharyngeal",
    "ileum",
    "integument",
    "intestine",
    "mandibular",
    "midgut",
    "midguts",
    "mouthparts",
    "muscle",
    "ovaries",
    "ovary",
    "rectum",
    "spermetheca",
    "spermethecia",
    "testicles",
    "testis",
    "thorax",
    "thoraxes",
    "thoracic",
    "tissue",
    "tubule",
    "tubules",
    "ventriculus",
    "wing",
    "wings"
  ),

  bee_species_types = c(
    "africanized",
    "apis",
    "capensis",
    "cerana",
    "domestic",
    "domesticated",
    "european",
    "honey",
    "honeybee",
    "honeybees",
    "italian",
    "korean",
    "koreana",
    "ligustica",
    "mellifera",
    "wild",
    "wildtype",
    "worker",
    "workers",
    "workerbee",
    "workerbees"
  ),

  life_stages_development = c(
    "adult",
    "adults",
    "brood",
    "capped",
    "egg",
    "eggs",
    "embryo",
    "embryonic",
    "embryos",
    "emerged",
    "emergence",
    "instar",
    "larva",
    "larvae",
    "larval",
    "mature",
    "oocyte",
    "pharate",
    "prepupa",
    "prepupal",
    "pupa",
    "pupae"
  ),

  bee_castes_roles = c(
    "drone",
    "forager",
    "foragers",
    "guard",
    "guards",
    "nurse",
    "queen",
    "queens",
    "queenless",
    "soldier",
    "virgin",
    "worker",
    "workers"
  ),

  parasites_pathogens = c(
    "ascosphaera",
    "destructor",
    "mite",
    "mites",
    "nosema",
    "varroa",
    "virus",
    "viral"
  ),

  diseases_conditions = c(
    "bqcv",
    "deformed",
    "disease",
    "diseased",
    "dwv",
    "iapv",
    "infected",
    "infection",
    "infested",
    "sbv"
  ),

  pesticides_chemicals = c(
    "acetmiprid",
    "amitraz",
    "carbaryl",
    "carbendazim",
    "chlorothalonil",
    "chlorphyrifos",
    "clothianidin",
    "coumophos",
    "curcumin",
    "deltamethrin",
    "dinotefuran",
    "dmso",
    "flupyradifurone",
    "fluvalinate",
    "glyfosphate",
    "halofuginone",
    "imidacloprid",
    "insecticide",
    "kresoxym",
    "pesticide",
    "quercetin",
    "rosmarinic",
    "thiacloprid",
    "thiamethoxam",
    "tunicamycin"
  ),

  treatments_interventions = c(
    "antibiotic",
    "crispr",
    "detoxificant",
    "inhibited",
    "injected",
    "knockdown",
    "knocked",
    "treated",
    "treatment",
    "untreatedbees"
  ),

  behaviors_traits = c(
    "aggressive",
    "highaggressioncolony",
    "lowaggressioncolony",
    "behavior",
    "grooming",
    "hygienic",
    "unresponsive",
    "vsh"
  ),

  colony_hive_terms = c(
    "apiary",
    "apiaries",
    "beehive",
    "cage",
    "caged",
    "colonies",
    "colony",
    "comb",
    "entrance",
    "frames",
    "hive",
    "hives",
    "nest"
  ),

  biological_processes = c(
    "colonization",
    "colonized",
    "development",
    "evolution",
    "feed",
    "feeding",
    "flight",
    "forage",
    "germline",
    "inseminated",
    "laying",
    "mated",
    "reproductive"
  ),

  molecular_biology_terms = c(
    "circrna",
    "cyp",
    "doublesex",
    "dsx",
    "gfp",
    "lncrna",
    "mirna",
    "mir",
    "mrna",
    "mrnas",
    "rnaseq",
    "scrna",
    "totalrna"
  ),

  research_methods = c(
    "arrayexpress",
    "blast",
    "bulkrna",
    "collected",
    "collection",
    "dissected",
    "experimental",
    "exposed",
    "fastq",
    "fasta",
    "homogenate",
    "homogenates",
    "library",
    "lysate",
    "pooled",
    "sample",
    "samples",
    "seq",
    "trimmed"
  ),

  genetics_genomics = c(
    "genotype",
    "genetics",
    "genomics",
    "haploid",
    "hybrid",
    "hybrids"
  ),

  organisms_microbiota = c(
    "alvi",
    "bifidobacterium",
    "gilliamella",
    "microbiota",
    "snodgrassella"
  ),

  substances_materials = c(
    "carbohydrates",
    "jelly",
    "nectar",
    "pheromone",
    "pollen",
    "sugar",
    "sucrose",
    "venom",
    "water"
  ),

  food_sources = c(
    "apple",
    "buckwheat",
    "cherry",
    "chestnut",
    "floral",
    "orange",
    "orchard",
    "pear",
    "rockrose",
    "tupelo"
  ),

  geographical_locations = c(
    "australia",
    "beijing",
    "brazil",
    "california",
    "cambridge",
    "china",
    "chinese",
    "edinburgh",
    "europe",
    "european",
    "fujian",
    "germany",
    "hawaii",
    "hawaiian",
    "illinois",
    "india",
    "israel",
    "israeli",
    "italy",
    "italian",
    "jiangxi",
    "korea",
    "korean",
    "scotland",
    "taiwan",
    "united kingdom"
  ),

  institutions_organizations = c(
    "university",
    "institute",
    "research",
    "college",
    "laboratory",
    "usda",
    "department"
  ),

  time_periods = c(
    "day",
    "days",
    "hour",
    "hours",
    "week",
    "weeks",
    "month",
    "months",
    "year",
    "spring",
    "summer",
    "fall",
    "morning",
    "afternoon"
  ),

  experimental_conditions = c(
    "cold",
    "control",
    "ctrl",
    "diet",
    "field",
    "stress",
    "temperature"
  )
)

# Function to extract keywords from a row
extract_keywords <- function(row_text, category_keywords) {
  # Convert row to lowercase text string
  text_lower <- tolower(paste(row_text, collapse = " "))

  # Find matching keywords
  matches <- category_keywords[str_detect(
    text_lower,
    regex(paste0("\\b", category_keywords, "\\b"), ignore_case = TRUE)
  )]

  # Return as comma-separated string or NA
  if (length(matches) > 0) {
    return(paste(unique(matches), collapse = ", "))
  } else {
    return(NA_character_)
  }
}

# Create df_meta by iterating through df, excluding long meta likely to be expt methods
df_meta <- df %>%
  filter(if_any(where(is.character), ~ nchar(.) < 100)) |>
  rowwise() %>%
  mutate(
    anatomical_terms = extract_keywords(
      c_across(everything()),
      keywords$anatomical_terms
    ),
    bee_species_types = extract_keywords(
      c_across(everything()),
      keywords$bee_species_types
    ),
    life_stages_development = extract_keywords(
      c_across(everything()),
      keywords$life_stages_development
    ),
    bee_castes_roles = extract_keywords(
      c_across(everything()),
      keywords$bee_castes_roles
    ),
    parasites_pathogens = extract_keywords(
      c_across(everything()),
      keywords$parasites_pathogens
    ),
    diseases_conditions = extract_keywords(
      c_across(everything()),
      keywords$diseases_conditions
    ),
    pesticides_chemicals = extract_keywords(
      c_across(everything()),
      keywords$pesticides_chemicals
    ),
    treatments_interventions = extract_keywords(
      c_across(everything()),
      keywords$treatments_interventions
    ),
    behaviors_traits = extract_keywords(
      c_across(everything()),
      keywords$behaviors_traits
    ),
    colony_hive_terms = extract_keywords(
      c_across(everything()),
      keywords$colony_hive_terms
    ),
    biological_processes = extract_keywords(
      c_across(everything()),
      keywords$biological_processes
    ),
    molecular_biology_terms = extract_keywords(
      c_across(everything()),
      keywords$molecular_biology_terms
    ),
    research_methods = extract_keywords(
      c_across(everything()),
      keywords$research_methods
    ),
    genetics_genomics = extract_keywords(
      c_across(everything()),
      keywords$genetics_genomics
    ),
    organisms_microbiota = extract_keywords(
      c_across(everything()),
      keywords$organisms_microbiota
    ),
    substances_materials = extract_keywords(
      c_across(everything()),
      keywords$substances_materials
    ),
    food_sources = extract_keywords(
      c_across(everything()),
      keywords$food_sources
    ),
    geographical_locations = extract_keywords(
      c_across(everything()),
      keywords$geographical_locations
    ),
    institutions_organizations = extract_keywords(
      c_across(everything()),
      keywords$institutions_organizations
    ),
    time_periods = extract_keywords(
      c_across(everything()),
      keywords$time_periods
    ),
    experimental_conditions = extract_keywords(
      c_across(everything()),
      keywords$experimental_conditions
    )
  ) %>%
  ungroup() %>%
  select(acc, everything(), -all_of(setdiff(names(df), "acc"))) # Remove original columns, keep only metadata

#### ---- Summary Data ----

path <- Sys.getenv("PATH_DATA")

summary_files <- list.files(
  path = path,
  recursive = T,
  pattern = "*.txt.summary",
  full.names = T
)
summary_list <- lapply(summary_files, read_delim, delim = '\t')
summary_combined <- reduce(summary_list, left_join, by = "Status")
write_csv(
  summary_combined,
  paste0(Sys.getenv("PATH_PARENT"), "summary_srr_0-4000_MB.csv"),
  quote = "none"
)

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

# Plot distributions of assigned and unassigned reads
ggplot(summary_norm) +
  geom_point(aes(x = Assigned, y = Unassigned), size = 0.3) +
  xlab("Proportion Assigned") +
  ylab("Proportion Unassigned") +
  ggtitle(paste0(
    "Summary of Reads for SRR Runs of Size 0-4000 MB (n = ",
    summary_norm$SRR_runs |> unique() |> length(),
    ")"
  )) +
  facet_wrap(~Unassigned_type)


#### ---- Counts Data ----

path <- "//wsl.localhost/Ubuntu-24.04/home/blake/sra/counts/"

counts_files <- list.files(
  path = path,
  recursive = T,
  pattern = "counts_[A-Z]+[0-9]+_*.*\\.txt$",
  full.names = T
)
counts_list <- lapply(counts_files, read_delim, delim = '\t', skip = 1)
counts_combined <- reduce(counts_list, left_join, by = "Geneid")

csv_file_name <- "counts_srr_0-4000_MB.csv"

write_csv(
  counts_combined,
  paste0(path, csv_file_name),
  quote = "none"
)

# Read in CSV to avoid re-processing data
counts_combined <- read_csv(
  paste0(path, csv_file_name),
)

counts_scaled_100 <-
  counts_combined |>
  mutate(across(
    where(is.numeric),
    ~ {
      if (
        all(. == 0, na.rm = TRUE) ||
          min(., na.rm = TRUE) == max(., na.rm = TRUE)
      ) {
        0
      } else {
        100 *
          (. - min(., na.rm = TRUE)) /
          (max(., na.rm = TRUE) - min(., na.rm = TRUE))
      }
    }
  )) |>
  rename_with(~ str_remove(., "\\.SAM"))

x <- counts_scaled_100 |> filter(Geneid == Sys.getenv("GENE1"))
x <- x |>
  pivot_longer(cols = -Geneid, names_to = "acc", values_to = "values")

ggplot(x) +
  geom_histogram(aes(x = values))


y <- crossing(df_meta, x, .name_repair = "unique") |>
  filter(str_detect(acc...24, acc...1))

z <-
  y |>
  select(
    Geneid,
    acc...1,
    acc...24,
    values,
    anatomical_terms,
    bee_castes_roles
  ) |>
  pivot_longer(
    cols = -c(Geneid, acc...1, acc...24, values),
    names_to = "Category",
    values_to = "Keyword"
  )

z_counts <-
  z |>
  group_by(Keyword) |>
  count()
z <- left_join(z, z_counts, by = "Keyword", keep = F)

# z <- y |>
#   pivot_longer(
#     cols = -c(Geneid, acc...1, acc...24, values),
#     names_to = "Category",
#     values_to = "Keyword"
#   )

max_y <- max(z$values, na.rm = TRUE)

library(ggh4x)

ggplot(data = z, aes(x = Keyword, y = values, col = n)) +
  geom_boxplot() +
  geom_text(
    aes(label = paste("n =", n), y = max_y),
    hjust = 1,
    size = 3,
    colour = "black"
  ) +
  scale_color_gradient(trans = "log", low = "red", high = "green") +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_wrap2(. ~ Category, ncol = 1, scales = "free_y") +
  force_panelsizes(rows = c(2, 1)) # first facet twice as tall

z |>
  summarise(
    na_count = sum(is.na(values)),
    inf_count = sum(is.infinite(values)),
    nan_count = sum(is.nan(values))
  )


counts_binned <-
  counts_scaled_100 |>
  pivot_longer(cols = -Geneid, names_to = "runs", values_to = "counts") |>
  pull(counts) |>
  data.frame()

ggplot(data = counts_binned) +
  geom_histogram(
    aes(
      x = pull.pivot_longer.counts_scaled_100..cols....Geneid..names_to....runs...
    ),
    binwidth = 1
  ) +
  xlab("count value (normalised to 100)") +
  ylab("frequency") +
  scale_y_continuous(
    #   transform = "log10",
    #   breaks = 10^(0:8), # gridline positions
    labels = scales::comma,
    #   minor_breaks = as.vector(outer(1:9, 10^(0:10), `*`))
  ) +
  ggtitle(
    label = paste0(
      "Histogram of all count values (sra runs = ",
      ncol(counts_scaled_100) - 1,
      ")"
    )
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_line(color = "grey70"),
    panel.grid.minor.y = element_line(color = "grey85")
  )


#### ---- RUBBISH ----

# Check list of srr runs agains data collected

srr_1 <- read_lines("srr_1-1000_MB.txt")
srr_2 <- read_lines("srr_1001-2000_MB.txt")
srr_3 <- read_lines("srr_2001-3000_MB.txt")
srr_4 <- read_lines("srr_3001-4000_MB.txt")
srr_collected <- names(counts_scaled_100)[-1]
length(srr_collected) -
  (length(srr_1) + length(srr_2) + length(srr_3) + length(srr_4)) # There were three that couldnt be processed

length(srr_1) - sum(srr_1 %in% srr_collected)
length(srr_2)
sum(srr_2 %in% srr_collected)
length(srr_3) - sum(srr_3 %in% srr_collected)

mean_by_keyword <-
  z %>%
  group_by(Keyword) %>%
  summarise(
    mean_value = mean(values, na.rm = TRUE),
    se = sd(values) / sqrt(n())
  ) %>%
  arrange(desc(mean_value)) |>
  head(n = 20)

ggplot(
  data = mean_by_keyword,
  aes(x = reorder(Keyword, -mean_value), y = mean_value)
) +
  geom_bar(stat = "identity") +
  geom_errorbar(
    aes(ymin = mean_value - se, ymax = mean_value + se),
    width = 0.2
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

group_by_keyword <-
  z %>%
  group_by(Keyword) %>%
  mutate(mean_value)
arrange(desc(mean_value)) |>
  head(n = 20)

ggplot(data = group_by_keyword) +
  geom_boxplot(aes(x = Keyword, y = values)) +
  coord_flip()
