#!/usr/bin/env Rscript

suppressMessages(library(tidyverse))
suppressMessages(library(parallel))
suppressMessages(source("scripts/SPARQL.R"))
suppressMessages(source("scripts/common.r"))

options(scipen = 100, digits = 4)
args <- commandArgs(TRUE)

# ------------------ Argument Defaults ------------------
ids <- NULL
endpoint <- NULL
outdir <- NULL
loop <- TRUE

# ------------------ Argument Parsing ------------------
while (loop) {
  if (args[1] == "--endpoint") {
    if (args[2] == "dev") {
      endpoint <- endpoint_dev
    }
  }

  if (args[1] == "--plants") {
    ids <- args[2]
  }

  if (args[1] == "--outdir") {
    outdir <- args[2]
  }

  if (length(args) > 1) {
    args <- args[2:length(args)]
  } else {
    loop <- FALSE
  }
}

if (is.null(ids) || is.null(outdir)) {
  stop("Both --plants and --outdir must be provided.")
}

# Construct output path
outFile <- file.path(outdir, "step2_cmp_for_plants.csv")
dir.create(dirname(outFile), recursive = TRUE, showWarnings = FALSE)

# ------------------ SPARQL: Get PLN Metadata ------------------
message("Verifying provided plant URIs...")

pln_uris <- trimws(unlist(str_split(ids, "\\|")))
quoted_uris <- paste(pln_uris, collapse = " ")

plant_query <- paste(sparql_prefix, paste0(
  "select distinct ?pln ?pln_label ?taxid ?rel where {
    values ?pln { ", quoted_uris, " }
    ?pln rdf:type sen:taxon .
    OPTIONAL { ?pln rdfs:label ?pln_label }
    OPTIONAL {
      values ?rel { rdfs:subClassOf sen:has_taxid }
      ?pln ?rel ?taxid .
      ?taxid ncbit:has_rank ?rank
    }
  }"
))

df.id <- SPARQL(endpoint, plant_query, ns = prefix, extra = query_options, format = "json")$results
if (nrow(df.id) == 0) {
  warning("No plants found. Exiting with empty result.")
  write_csv(tibble(), outFile)
  quit("no", 0)
}


# ------------------ SPARQL: Get Compounds for Plants ------------------
message("Fetching compounds associated with plants...")

getCompoundsForPlant <- function(SENPLN) {
  q <- paste(sparql_prefix, paste0(
    "SELECT distinct ?pln ?cmp
      (group_concat(distinct(?cmp_label); separator=\"|\") as ?cmp_labels)
      (group_concat(distinct(?part_name); separator=\"|\") as ?plant_parts)
      (group_concat(distinct(?reference); separator=\"|\") as ?references)
      (group_concat(distinct(?pcx_src_label); separator=\"|\") as ?plant_compound_source)
    WHERE {
      values ?pln { ", SENPLN, " } .
      ?pln ^sen:inPlant ?pcx .
      ?pcx sen:hasCompound ?cmp .
      ?cmp rdfs:label ?cmp_label .
      ?pcx sen:hasSource ?pcx_src .
      ?pcx_src rdfs:label ?pcx_src_label .
      OPTIONAL {
        ?pcx sen:hasReference ?ref .
        ?ref rdfs:label ?reference
      }
      OPTIONAL {
        ?prt ^sen:inComponent ?pcx .
        ?prt rdfs:label ?part_name
      }
    }
    group by ?pln ?cmp
    order by ?part_name"
  ))
  SPARQL(endpoint, q, ns = prefix, extra = query_options, format = "json")$results
}

compound_results <- mclapply(df.id$pln, getCompoundsForPlant, mc.cores = 10)
combined <- suppressWarnings(fix_sparql_ids(bind_rows(compound_results)))

# ------------------ Safe Join ------------------
if ("pln" %in% names(combined) && "pln" %in% names(df.id)) {
  message("Joining plant metadata with compound results")
  combined <- combined %>%
    left_join(df.id, by = "pln")
} else {
  warning("Join skipped: 'pln' column missing in one or both data frames.")
}

# ------------------ Save Output ------------------
if (nrow(combined) == 0) {
  message("No compound mappings found.")
} else {
  message(paste("Total rows retrieved:", nrow(combined)))
}
write_csv(combined, outFile)
