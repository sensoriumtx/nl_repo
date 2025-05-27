suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(parallel))
suppressMessages(source("scripts/SPARQL.R"))
suppressMessages(source("scripts/common.r"))

options(scipen = 100, digits = 4)
args <- commandArgs(TRUE)

ids = NULL
outFile = NULL
loop = TRUE

# Argument Parsing
while (loop) {
  if (args[1] == "--endpoint") {
    if (args[2] == "dev") {
      endpoint = endpoint_dev
    }
  }

  if (args[1] == "--plants") {
    ids = args[2]
  }

  if (args[1] == "--out") {
    outFile = args[2]
  }

  if (length(args) > 1) {
    args <- args[2:length(args)]
  } else {
    loop = FALSE
  }
}

# Validate Arguments
if (is.null(ids) | is.null(outFile)) {
  stop("Error: Missing required arguments (--plants or --out).")
}

# Make output directory if it doesn't exist
if (!dir.exists(dirname(outFile))) dir.create(dirname(outFile), recursive = TRUE)

# Split Plant Labels and Deduplicate
plant_labels <- unique(trimws(unlist(str_split(ids, r"(\|)"))))
cat("Total Plant Labels:", length(plant_labels), "\n")

# Set Parallel Options
num_cores <- detectCores() - 1
cat("Using", num_cores, "parallel workers.\n")

# Chunking Plant Labels for Parallel Querying (Batch Size 100)
chunk_size <- 100
plant_chunks <- split(plant_labels, ceiling(seq_along(plant_labels) / chunk_size))
cat("Total Chunks:", length(plant_chunks), "\n")

# Optimized SPARQL Query Function (Vectorized)
getCompoundsForPlants <- function(plants) {
  cat("\n[Chunk] Querying Compounds for Plants:", length(plants), "labels\n")
  
  query <- paste0(sparql_prefix, "
    SELECT DISTINCT ?pln ?cmp ?cmp_label ?part_name ?reference ?pcx_src_label WHERE {
      VALUES ?pln { ", paste0('<', plants, '>', collapse = " "), " }
      ?pln ^sen:inPlant ?pcx .
      ?pcx sen:hasCompound ?cmp .
      ?cmp rdfs:label ?cmp_label .
      OPTIONAL { 
        ?pcx sen:hasSource ?pcx_src .
        ?pcx_src rdfs:label ?pcx_src_label .
      }
      OPTIONAL { 
        ?prt ^sen:inComponent ?pcx .
        ?prt rdfs:label ?part_name .
      }
      OPTIONAL { 
        ?pcx sen:hasReference ?ref .
        ?ref rdfs:label ?reference .
      }
    }
  ")
  
  # Execute SPARQL Query
  res <- SPARQL(endpoint, query, ns = prefix, extra = query_options, format = "json")$results
  
  if (nrow(res) == 0) {
    return(NULL)
  } else {
    return(res)
  }
}

# Parallel Execution
cl <- makeCluster(num_cores)
clusterExport(cl, list("getCompoundsForPlants", "endpoint", "sparql_prefix", "prefix", "query_options"))
results_list <- parLapply(cl, plant_chunks, getCompoundsForPlants)
stopCluster(cl)

# Combine Results Efficiently
cat("\nCombining Results...\n")
results_combined <- rbindlist(results_list, fill = TRUE)

# Write Results Efficiently
fwrite(results_combined, outFile)
cat("\n[Success] Results written to:", outFile, "\n")
