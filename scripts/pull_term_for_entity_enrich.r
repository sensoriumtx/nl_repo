# Unified Act/Use Label Enrichment Script for Compounds and Plants
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
    } else if (arg == "--raw_out") {
      out$rawOut <- val
    } else if (arg == "--use_targets") {
      out$useTargets <- TRUE
      i <- i - 1
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
rawOut <- cli_args$rawOut
ids <- cli_args$ids
useTargets <- ifelse(is.null(cli_args$useTargets), FALSE, TRUE)

if (is.null(outFile)) stop("--out must be specified")
if (!dir.exists(dirname(outFile))) dir.create(dirname(outFile), recursive = TRUE)

# Step 1: Gather terms
gather_terms <- function() {
  if (!is.null(terms_raw)) return(trimws(unlist(str_split(terms_raw, "\\|"))))
  if (!is.null(in_file)) {
    df <- read_csv(in_file, show_col_types = FALSE)
    for (colname in c("term", "Term", "terms", "Terms")) {
      if (colname %in% colnames(df)) return(trimws(df[[colname]]))
    }
    stop("No valid term column found in input file")
  }
  stop("Either --terms or --in_file must be specified")
}

terms_input <- gather_terms()

# Step 2: Normalize and resolve ?use URI directly
normalize_term <- function(x) str_to_lower(stri_trans_general(trimws(x), "Latin-ASCII"))

resolve_term_uris <- function(term) {
  term_norm <- normalize_term(term)
  q <- paste(sparql_prefix, sprintf('SELECT ?use ?label WHERE { ?use a sen:use; rdfs:label ?label . FILTER(LCASE(STR(?label)) = "%s") }', term_norm))
  res <- tryCatch(SPARQL(endpoint, q, ns=prefix, extra=query_options, format='json')$results, error = function(e) NULL)
  if (!is.null(res) && is.data.frame(res) && nrow(res) > 0 && !is.null(res$use[1])) return(res$use[1])
  return(NA_character_)
}

term_map <- tibble(term = terms_input) %>%
  mutate(use_uri = suppressWarnings(map_chr(term, resolve_term_uris)))

valid_terms <- term_map %>% filter(!is.na(use_uri))
unmatched_terms <- setdiff(terms_input, valid_terms$term)
if (length(unmatched_terms) > 0) message(paste0("\u26A0\uFE0F Unmatched terms: ", paste(unmatched_terms, collapse = "|")))

# Step 3: Enrich from each resolved ?use
if (nrow(valid_terms) == 0) stop("No valid compound or plant mappings found. Exiting.")

pull_enrichment <- function(uri) {
  q <- paste(sparql_prefix, sprintf('SELECT DISTINCT ?cmp ?cmp_label ?pln ?pln_label ?act_label ?use ?use_label WHERE {
    BIND(<%s> AS ?use)
    ?use rdfs:label ?use_label .
    OPTIONAL {
      ?cmp ^sen:hasCompound|(^sen:hasCompound/(sen:targetProtein|sen:targetGene)/^sen:hasTarget) ?use .
      ?cmp (sen:lcLabel|sen:lcAltLabel)|(sen:maps_to/(sen:lcLabel|sen:lcAltLabel)) ?cmp_label .
    }
    OPTIONAL {
      ?use sen:hasActivity/rdfs:label ?act_label .
    }
    OPTIONAL {
      ?use sen:hasTaxon ?pln .
      ?pln rdfs:label ?pln_label .
    }
  }', uri))

  res <- tryCatch(SPARQL(endpoint, q, ns=prefix, extra=query_options, format='json')$results, error = function(e) NULL)

  if (!is.null(res) && is.data.frame(res) && nrow(res) > 0) {
    res$use_uri <- uri
    return(as_tibble(res))
  } else {
    return(tibble())
  }
}

results <- valid_terms$use_uri %>% map(pull_enrichment) %>% bind_rows()
if (!is.null(results) && nrow(results) > 0) {
  write_csv(results, outFile)
} else {
  message("\u26A0\uFE0F No enrichment results returned")
}
