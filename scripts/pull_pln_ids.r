suppressMessages(library(tidyverse))
suppressMessages(library(stringi))
suppressMessages(library(parallel))
suppressMessages(source("scripts/common.r"))
suppressMessages(source("scripts/SPARQL.R"))

options(scipen = 100, digits = 4)
args <- commandArgs(TRUE)

# ids = "Psilocin|Mesembrine"
outFile = NULL
loop = TRUE

while (loop) {
	if (args[1] == "--endpoint") {
		if (args[2] == "dev") {
      		endpoint = endpoint_dev
    	}
	}

	if (args[1] == "--out") {
		outFile = args[2]
	}

	if (length(args)>1) {
		args<-args[2:length(args)]
	} else {
		loop = FALSE
	}
}

if (
    is.null(outFile)
) {
    q("no", 1, FALSE)
}

# make out folder if doesn't exist
if (!dir.exists(dirname(outFile))) dir.create(dirname(outFile), recursive = TRUE)

message(paste0("endpoint: ", endpoint))

# split cmps
q = paste(sparql_prefix, paste0(
"select distinct ?pln ?pln_label ?taxid ?rel where {
    ?pln ^sen:hasTaxon ?use .
    ?use a sen:use .
    ?pln rdfs:label ?pln_label
    OPTIONAL {
        values ?rel { rdfs:subClassOf sen:has_taxid }
        ?pln ?rel ?taxid .
        ?taxid ncbit:has_rank ?rank
    }
}
"))
df.cmp_ids <- fix_sparql_ids(SPARQL(endpoint, q, ns=prefix, extra=query_options, format='json')$results)

write_csv(df.cmp_ids, outFile)