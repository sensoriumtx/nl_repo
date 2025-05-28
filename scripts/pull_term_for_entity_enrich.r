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

# Step 2: Infer term types with normalization
normalize_term <- function(x) str_to_lower(stri_trans_general(trimws(x), "Latin-ASCII"))

infer_term_type <- function(term) {
  term_norm <- normalize_term(term)
  q_act <- paste(sparql_prefix, sprintf('SELECT ?label WHERE { ?x a sen:activity; sen:lcLabel ?label . FILTER(LCASE(STR(?label)) = "%s") } LIMIT 1', term_norm))
  q_use <- paste(sparql_prefix, sprintf('SELECT ?label WHERE { ?x a sen:use; rdfs:label ?label . FILTER(LCASE(STR(?label)) = "%s") } LIMIT 1', term_norm))
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
if (length(unmatched_terms) > 0) message(paste0("\u26A0\uFE0F Unmatched terms: ", paste(unmatched_terms, collapse = "|")))

# Step 3: Get cmp and pln mappings
get_compounds_from_terms <- function(term, type) {
  q <- if (type == "act") paste0(sparql_prefix, sprintf('SELECT DISTINCT ?cmp ?cmp_label WHERE { ?use sen:hasActivity ?act . ?act sen:lcLabel ?l . FILTER(LCASE(STR(?l)) = "', normalize_term(term), '") ?cmp (^sen:hasCompound|(sen:maps_to+/^sen:hasCompound)) ?use . ?cmp rdfs:label ?cmp_label }')) else paste0(sparql_prefix, sprintf('SELECT DISTINCT ?cmp ?cmp_label WHERE { ?use a sen:use ; rdfs:label ?l . FILTER(LCASE(STR(?l)) = "', normalize_term(term), '") ?cmp (^sen:hasCompound|(sen:maps_to+/^sen:hasCompound)) ?use . ?cmp rdfs:label ?cmp_label }'))
  df <- SPARQL(endpoint, q, ns=prefix, extra=query_options, format='json')$results
  df$term <- term
  df$term_type <- type
  return(df)
}

get_plants_from_terms <- function(term, type) {
  q <- if (type == "act") paste0(sparql_prefix, sprintf('SELECT DISTINCT ?pln ?pln_label WHERE { ?use sen:hasActivity ?act . ?act sen:lcLabel ?l . FILTER(LCASE(STR(?l)) = "', normalize_term(term), '") ?pln ^sen:hasTaxon ?use . ?pln rdfs:label ?pln_label }')) else paste0(sparql_prefix, sprintf('SELECT DISTINCT ?pln ?pln_label WHERE { ?use a sen:use ; rdfs:label ?l . FILTER(LCASE(STR(?l)) = "', normalize_term(term), '") ?pln ^sen:hasTaxon ?use . ?pln rdfs:label ?pln_label }'))
  df <- SPARQL(endpoint, q, ns=prefix, extra=query_options, format='json')$results
  df$term <- term
  df$term_type <- type
  return(df)
}

compound_term_links <- if (nrow(term_map) > 0) {
  bind_rows(mclapply(1:nrow(term_map), function(i) {
    get_compounds_from_terms(term_map$term[i], term_map$term_type[i])
  }, mc.cores = 10))
} else {
  tibble()
}

plant_term_links <- if (nrow(term_map) > 0) {
  bind_rows(mclapply(1:nrow(term_map), function(i) {
    get_plants_from_terms(term_map$term[i], term_map$term_type[i])
  }, mc.cores = 10))
} else {
  tibble()
}

if (!is.null(rawOut)) {
  write_csv(compound_term_links, paste0(rawOut, "_compound_term_links.csv"))
  write_csv(plant_term_links, paste0(rawOut, "_plant_term_links.csv"))
}

if (nrow(compound_term_links) == 0 && nrow(plant_term_links) == 0) {
  message("No valid compound or plant mappings found. Exiting.")
  quit("no", 0)
}

# Step 4: Calculate enrichment stats for cmp and pln
message("Calculating enrichment statistics...")

cmp_counts <- compound_term_links %>% group_by(term, term_type, cmp = cmp) %>% summarise(count = n(), .groups = 'drop')
pln_counts <- plant_term_links %>% group_by(term, term_type, pln = pln) %>% summarise(count = n(), .groups = 'drop')

total_cmp <- as.integer(SPARQL(endpoint, paste0(sparql_prefix, 'SELECT (COUNT(DISTINCT ?cmp) as ?count) WHERE { ?cmp a sen:compound }'), ns=prefix, extra=query_options)$results$count[1])
total_pln <- as.integer(SPARQL(endpoint, paste0(sparql_prefix, 'SELECT (COUNT(DISTINCT ?pln) as ?count) WHERE { ?pln a sen:taxon }'), ns=prefix, extra=query_options)$results$count[1])

cmp_enrich <- cmp_counts %>%
  group_by(term, term_type) %>%
  summarise(n = n(), pval = phyper(n() - 1, n(), total_cmp - n(), total_cmp, lower.tail = FALSE)) %>%
  ungroup %>%
  mutate(padj = p.adjust(pval, method = "fdr"), entity_type = "compound")

pln_enrich <- pln_counts %>%
  group_by(term, term_type) %>%
  summarise(n = n(), pval = phyper(n() - 1, n(), total_pln - n(), total_pln, lower.tail = FALSE)) %>%
  ungroup %>%
  mutate(padj = p.adjust(pval, method = "fdr"), entity_type = "plant")

final <- bind_rows(cmp_enrich, pln_enrich)

# Step 5: Output results
message(paste0("Writing enrichment results to: ", outFile))
write_csv(final, outFile)
message("Done.")
