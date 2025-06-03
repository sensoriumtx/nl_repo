#!/usr/bin/env Rscript

suppressMessages(library(tidyverse))
suppressMessages(library(stringi))
suppressMessages(library(parallel))
suppressMessages(source("scripts/common.r"))
suppressMessages(source("scripts/SPARQL.R"))

options(scipen = 100, digits = 4)

# -------------------- Argument Parsing --------------------
args <- commandArgs(trailingOnly = TRUE)
parsed_args <- list(endpoint = "dev", cmp_file = NULL, cmp_column = NULL, chunk_size = 50, out = NULL)

i <- 1
while (i <= length(args)) {
  key <- args[i]
  val <- args[i + 1]
  if (key == "--endpoint") parsed_args$endpoint <- val
  if (key == "--cmp_file") parsed_args$cmp_file <- val
  if (key == "--cmp_column") parsed_args$cmp_column <- val
  if (key == "--chunk_size") parsed_args$chunk_size <- as.integer(val)
  if (key == "--out") parsed_args$out <- val
  i <- i + 2
}

# -------------------- Input Validation --------------------
if (is.null(parsed_args$cmp_file) || is.null(parsed_args$cmp_column) || is.null(parsed_args$out)) {
  stop("Missing required arguments: --cmp_file, --cmp_column, or --out", call. = FALSE)
}

endpoint <- if (parsed_args$endpoint == "dev") endpoint_dev else stop("Unknown endpoint")

cmp_list <- read_csv(parsed_args$cmp_file, show_col_types = FALSE) %>%
  pull(parsed_args$cmp_column) %>%
  unique() %>%
  discard(is.na)

if (length(cmp_list) == 0) stop("No valid cmp values found.", call. = FALSE)

# -------------------- Split Into Chunks --------------------
chunks <- split(cmp_list, ceiling(seq_along(cmp_list) / parsed_args$chunk_size))

# -------------------- SPARQL Queries --------------------
outcomes <- c(
  "\"Agonist\"", "\"Inhibitor\"", "\"Antagonist\"", "\"Activator\"",
  "\"Inhibition\"", "\"Blocker\"", "\"Channel blocker\"", "\"Blocker (channel blocker)\"",
  "\"Full agonist\"", "\"Partial agonist\"", "\"Activation\"", "\"Inhibitor (gating inhibitor)\"",
  "\"Stimulator\"", "\"Inducer\"", "\"Positive\"", "\"Pore blocker\"", "\"Suppressor\"",
  "\"Potentiation\"", "\"Gating inhibitor\"", "\"Opener\"", "\"Active\"", "\"Unspecified\""
)

get_cmp_metadata <- function(chunk) {
  ids_str <- paste(chunk, collapse = " ")
  q <- paste0(sparql_prefix, "
    select distinct ?cmp (sample(?label) as ?cmp_label) where {
      values ?cmp { ", ids_str, " }
      {
        ?cmp (^sen:hasCompound|(sen:maps_to+/^sen:hasCompound))/rdf:type sen:use .
      } UNION {
        values ?outcome { ", paste(outcomes, collapse = " "), " }
        ?assayResult a sen:assayResult .
        ?assayResult sen:measureOutcome ?outcome .
        ?cmp ^sen:hasCompound ?assayResult .
        ?assayResult sen:targetGene|sen:targetProtein ?target .
        ?target ^sen:hasTarget/rdf:type sen:use .
      }
      ?cmp sen:lcLabel ?label
    }
    group by ?cmp
  ")
  SPARQL(endpoint, q, ns = prefix, extra = query_options, format = 'json')$results
}

get_cmp_activities <- function(cmp_id) {
  q <- paste0(sparql_prefix, "
    select distinct ?cmp ?isoSmiles ?inchi ?act ?act_label (group_concat(distinct ?path; separator=\"|\") as ?path)
           (count(distinct ?tgt) as ?target_count_for_this_compound_activity)
           (group_concat(distinct ?tgt; separator=\"|\") as ?target) ?assay_source where {
      select ?assay_source ?cmp ?isoSmiles ?inchi ?act ?act_label ?path ?tgt ?primary_symbol where {
        ?use a sen:use .
        values ?cmp { ", cmp_id, " } .
        {
          ?cmp (^sen:hasCompound|(sen:maps_to+/^sen:hasCompound)) ?use .
        } UNION {
          values ?outcome { ", paste(outcomes, collapse = " "), " }
          ?assayResult a sen:assayResult .
          ?assayResult sen:measureOutcome ?outcome .
          ?cmp ^sen:hasCompound ?assayResult .
          ?assayResult (sen:targetGene/sen:maps_to*)|(sen:targetProtein/sen:maps_to*) ?tgt_entity .
          ?use sen:hasTarget/sen:maps_to* ?tgt_entity .
          ?tgt_entity sen:primary_symbol ?tgt .
          ?tgt_entity sen:tox21 \"false\" .
          ?assayResult sen:hasSource ?src .
          ?src rdfs:label ?assay_source .
        }
        ?use sen:hasActivity ?act .
        ?act sen:lcLabel ?act_label .
        OPTIONAL { ?cmp sen:IsomericSMILES ?isoSmiles . }
        OPTIONAL { ?cmp sen:InChI ?inchi . }
        BIND(IF(bound(?tgt), \"indirect\", \"direct\") as ?path)
      } order by ?path ?cmp ?act
    } group by ?assay_source ?cmp ?isoSmiles ?inchi ?act ?act_label
  ")
  SPARQL(endpoint, q, ns = prefix, extra = query_options, format = 'json')$results
}

# -------------------- Execute Queries in Parallel --------------------
all_metadata <- bind_rows(mclapply(chunks, get_cmp_metadata, mc.cores = detectCores()))
all_activities <- bind_rows(mclapply(all_metadata$cmp, get_cmp_activities, mc.cores = detectCores()))

# -------------------- Final Output --------------------
final_df <- fix_sparql_ids(all_activities) %>%
  left_join(all_metadata, by = "cmp")

if (!dir.exists(dirname(parsed_args$out))) dir.create(dirname(parsed_args$out), recursive = TRUE)
write_csv(final_df, parsed_args$out)
message("✅ Output written to ", parsed_args$out)
