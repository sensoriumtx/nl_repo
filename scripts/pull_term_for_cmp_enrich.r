# Unified Act/Use Label Enrichment Script
suppressMessages(library(tidyverse))
suppressMessages(library(stringi))
suppressMessages(library(parallel))
suppressMessages(source("scripts/common.r"))
suppressMessages(source("scripts/SPARQL.R"))

options(scipen = 100, digits = 4)

parseArgs <- function(args) {
  out <- list()
  i <- 1
  while (i <= length(args)) {
    arg <- args[i]
    val <- if (i + 1 <= length(args)) args[i + 1] else NA
    if (is.na(val)) stop(paste("Missing value for", arg))
    if (arg == "--endpoint") {
      out$endpoint <- if (val == "prod") endpoint_prod else endpoint_dev
    } else if (arg == "--terms") {
      out$terms_raw <- val
    } else if (arg == "--in_file") {
      out$in_file <- val
    } else if (arg == "--filter_out_file") {
      out$filter_out_file <- val
    } else if (arg == "--compound") {
      out$ids <- val
    } else if (arg == "--out") {
      out$outFile <- val
    } else if (arg == "--use_targets") {
      out$useTargets <- TRUE
      i <- i - 1  # flag argument has no value
    }
    i <- i + 2
  }
  return(out)
}

cli_args <- parseArgs(commandArgs(TRUE))

endpoint <- cli_args$endpoint
terms_raw <- cli_args$terms_raw
in_file <- cli_args$in_file
filter_out_file <- cli_args$filter_out_file
outFile <- cli_args$outFile
ids <- cli_args$ids
useTargets <- ifelse(is.null(cli_args$useTargets), FALSE, TRUE)

if (is.null(outFile)) {
  stop("--out must be specified")
}

# Create output directory if needed
if (!dir.exists(dirname(outFile))) dir.create(dirname(outFile), recursive = TRUE)

# Step 1: Gather terms

gather_terms <- function() {
  if (!is.null(terms_raw)) {
    return(trimws(unlist(str_split(terms_raw, "\\|"))))
  }
  if (!is.null(in_file)) {
    df <- read_csv(in_file, show_col_types = FALSE)
    for (colname in c("term", "Term", "terms", "Terms")) {
      if (colname %in% colnames(df)) {
        return(trimws(df[[colname]]))
      }
    }
    stop("No valid term column found in input file")
  }
  stop("Either --terms or --in_file must be specified")
}

terms_input <- gather_terms()

# Step 2: Infer term types

infer_term_type <- function(term) {
  q_act <- paste(sparql_prefix, sprintf('SELECT ?label WHERE { ?x a sen:activity; sen:lcLabel "%s" } LIMIT 1', term))
  q_use <- paste(sparql_prefix, sprintf('SELECT ?label WHERE { ?x a sen:use; rdfs:label "%s" } LIMIT 1', term))
  is_act <- nrow(SPARQL(endpoint, q_act, ns=prefix, extra=query_options, format='json')$results) > 0
  is_use <- nrow(SPARQL(endpoint, q_use, ns=prefix, extra=query_options, format='json')$results) > 0
  if (is_act) return("act")
  if (is_use) return("use")
  return(NA)
}

term_map <- tibble(term = terms_input) %>%
  mutate(term_type = unlist(mclapply(term, infer_term_type, mc.cores = 10))) %>%
  filter(!is.na(term_type))

unmatched_terms <- setdiff(terms_input, term_map$term)

# Step 3: Output warning if unmatched
if (length(unmatched_terms) > 0) {
  warning_msg <- paste0("\u26A0\uFE0F Unmatched terms: ", paste(unmatched_terms, collapse = "|"))
  message(warning_msg)
}

# Step 4: Enrichment counts (placeholder - insert core logic here)
# For each term type, build the query logic and join results

message("Resolving compounds for enrichment analysis...")

get_compounds_from_terms <- function(term, type) {
  if (type == "act") {
    q <- paste0(sparql_prefix, sprintf('
    SELECT DISTINCT ?cmp ?cmp_label ?label WHERE {
      ?use sen:hasActivity ?act .
      ?act sen:lcLabel "%s" .
      ?cmp (^sen:hasCompound|(sen:maps_to+/^sen:hasCompound)) ?use .
      ?cmp rdfs:label ?cmp_label
    }', term))
  } else {
    q <- paste0(sparql_prefix, sprintf('
    SELECT DISTINCT ?cmp ?cmp_label ?label WHERE {
      ?use a sen:use ; rdfs:label "%s" .
      ?cmp (^sen:hasCompound|(sen:maps_to+/^sen:hasCompound)) ?use .
      ?cmp rdfs:label ?cmp_label
    }', term))
  }
  df <- SPARQL(endpoint, q, ns=prefix, extra=query_options, format='json')$results
  df$term <- term
  df$term_type <- type
  return(df)
}

compound_term_links <- bind_rows(mclapply(1:nrow(term_map), function(i) {
  get_compounds_from_terms(term_map$term[i], term_map$term_type[i])
}, mc.cores = 10))

message("Calculating enrichment statistics...")

counts <- compound_term_links %>%
  group_by(term, term_type, cmp) %>%
  summarise(count = n(), .groups = 'drop')

total_cmp_count_q <- paste0(sparql_prefix, 'SELECT (COUNT(DISTINCT ?cmp) as ?count) WHERE { ?cmp a sen:compound }')
total_cmp <- as.integer(SPARQL(endpoint, total_cmp_count_q, ns=prefix, extra=query_options)$results$count[1])

enrichment_stats <- counts %>%
  group_by(term, term_type) %>%
  summarise(
    n = n(),
    pval = phyper(n() - 1, n(), total_cmp - n(), total_cmp, lower.tail = FALSE)
  ) %>%
  ungroup %>%
  mutate(padj = p.adjust(pval, method = "fdr"))

# Step 5: Output results
message(paste0("Writing enrichment results to: ", outFile))
write_csv(enrichment_stats, outFile)
message("Done.")
