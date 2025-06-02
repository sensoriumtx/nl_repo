#!/usr/bin/env Rscript

### Does not work. SPARQL query not extracting correct identifier.


options(warn = -1)  # suppress all warnings
suppressMessages(library(tidyverse))
suppressMessages(library(parallel))
suppressMessages(source("scripts/common.r"))
suppressMessages(source("scripts/SPARQL.R"))

options(scipen = 100, digits = 4)
args <- commandArgs(TRUE)

inFile <- NULL
outFile <- NULL
endpoint <- NULL
numCores <- 10  # Fixed core usage

# Parse command-line arguments
loop <- TRUE
while (loop) {
  if (args[1] == "--endpoint") {
    if (args[2] == "dev") {
      endpoint <- endpoint_dev
    }
  }
  if (args[1] == "--in") {
    inFile <- args[2]
  }
  if (args[1] == "--out") {
    outFile <- args[2]
  }
  if (length(args) > 1) {
    args <- args[2:length(args)]
  } else {
    loop <- FALSE
  }
}

# Input validation
if (is.null(inFile) || is.null(outFile)) stop("Both --in and --out arguments are required.")
if (!file.exists(inFile)) stop(paste("Input file not found:", inFile))
if (!dir.exists(dirname(outFile))) {
  dir.create(dirname(outFile), recursive = TRUE)
  message("Created output directory: ", dirname(outFile))
}

message("Reading input: ", inFile)
df_input <- read_csv(inFile, show_col_types = FALSE)

# Drop fully blank rows in input
df_input <- df_input %>%
  filter(!if_all(everything(), ~ is.na(.) || . == ""))

if (!"cmp" %in% names(df_input)) {
  stop("Input file must contain a 'cmp' column.")
}

# Extract unique cmp values
unique_cmp <- df_input$cmp %>%
  na.omit() %>%
  trimws() %>%
  discard(~ .x == "") %>%
  unique()

if (length(unique_cmp) == 0) {
  stop("No valid cmp identifiers found in input.")
}

total <- length(unique_cmp)
message("Found ", total, " unique cmp identifiers.")
message("Fetching Classification for Compounds")

# Function to query SPARQL classification per compound
getCompoundClassification <- function(cmp_uri, index = NA) {
  q <- paste(sparql_prefix, sprintf(
    '
    SELECT DISTINCT ?cmp ?cmpLabel ?kingdom ?superclass ?class
    WHERE {
      VALUES ?cmp { %s }
      ?cmp sen:lcLabel ?cmpLabel .
      OPTIONAL { ?cmp sen:ClassyFire_kingdom ?kingdom . }
      OPTIONAL { ?cmp sen:ClassyFire_superclass ?superclass . }
      OPTIONAL { ?cmp sen:ClassyFire_class ?class . }
    }
    ', cmp_uri))
  
  for (i in 1:3) {
    result <- tryCatch({
      raw <- SPARQL(endpoint, q, ns = prefix, extra = query_options, format = "json")
      return(raw$results)
    }, error = function(e) {
      message("Attempt ", i, " failed for ", cmp_uri, ": ", conditionMessage(e))
      Sys.sleep(2^i)
      return(NULL)
    })
    if (!is.null(result)) return(result)
  }
  
  message("All attempts failed for ", cmp_uri)
  return(data.frame())
}

# Parallel classification queries
results <- bind_rows(mclapply(seq_along(unique_cmp), function(i) {
  getCompoundClassification(unique_cmp[i], index = i)
}, mc.cores = numCores))

# Deduplicate and collapse classification
results_dedup <- results %>%
  select(any_of(c("cmp", "cmpLabel", "kingdom", "superclass", "class"))) %>%
  group_by(cmp) %>%
  summarize(
    cmpLabel = paste(sort(unique(na.omit(cmpLabel))), collapse = "; "),
    kingdom = first(na.omit(kingdom)),
    superclass = first(na.omit(superclass)),
    class = first(na.omit(class)),
    .groups = "drop"
  )

message("Retrieved classification for ", nrow(results_dedup), " unique compounds.")

# Merge classification with original input
df_merged <- df_input %>%
  distinct() %>%
  left_join(results_dedup, by = "cmp") %>%
  filter(!if_all(everything(), ~ is.na(.) || . == ""))

# Write output
message("Writing final output to: ", outFile)
write_csv(df_merged, outFile)
message("Complete. Final output has ", nrow(df_merged), " rows.")
