library(tidyverse)
source("scripts/SPARQL.R")
source("scripts/common.r")

options(scipen = 100, digits = 4)
args <- commandArgs(TRUE)

ids = NULL
senIds = NULL
outFile = NULL
compoundSmiles = NULL
cmpActFile = NULL
loop = TRUE
cmpIdCol = "cmp"

while (loop) {
	if (args[1] == "--endpoint") {
		if (args[2] == "dev") {
      		endpoint = endpoint_dev
    	}
	}

	if (args[1] == "--pubchem_ids") {
		ids = args[2]
	}

	if (args[1] == "--compound_names") {
		compoundSmiles = args[2]
	}

	if (args[1] == "--cmp_id_column") {
		cmpIdCol = args[2]
	}

	if (args[1] == "--compound_activity_file") {
		cmpActFile = args[2]
	}

	if (args[1] == "--sen_ids") {
		senIds = args[2]
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
    is.null(ids) &
    is.null(outFile)
) {
    q("no", 1, FALSE)
}

# make out folder if doesn't exist
if (!dir.exists(dirname(outFile))) dir.create(dirname(outFile), recursive = TRUE)

message("fetching all plants")
q = paste(sparql_prefix, paste0(
"select distinct ?pln ?pln_label ?taxid ?rel where {
    ?pln a sen:taxon .
    ?pln rdfs:label ?pln_label
    OPTIONAL {
        ?pln ^sen:hasTaxon ?use .
        ?use a sen:use .
    }
    OPTIONAL {
        values ?rel { rdfs:subClassOf sen:has_taxid }
        ?pln ?rel ?taxid .
        ?taxid ncbit:has_rank ?rank
    }
}
"))
suppressWarnings(
df.id <- SPARQL(endpoint, q, ns=prefix, extra=query_options, format='json')$results %>%
  group_by(pln) %>%
  summarize(
    pln_label = pln_label[which.max(nchar(pln_label))], .groups = 'drop',
    taxid,
    rel
  ) %>%
  unique %>%
  pivot_wider(names_from = "rel", values_from = "taxid") %>%
  mutate(across(starts_with("sen:") | starts_with("rdfs:") | starts_with("NA"), ~map_chr(., toString))) %>%
  mutate(
    taxid = case_when(
        `sen:has_taxid` == "" & `rdfs:subClassOf` == "" ~ NA_character_,
        `sen:has_taxid` != "" & `rdfs:subClassOf` == "" ~ `sen:has_taxid`,
        `sen:has_taxid` == "" & `rdfs:subClassOf` != "" ~ `rdfs:subClassOf`
    ),
    rel = case_when(
        `sen:has_taxid` == "" & `rdfs:subClassOf` == "" ~ NA_character_,
        `sen:has_taxid` != "" & `rdfs:subClassOf` == "" ~ "sen:has_taxid",
        `sen:has_taxid` == "" & `rdfs:subClassOf` != "" ~ "rdfs:subClassOf"
    )
  ) %>%
  dplyr::select(-c(`sen:has_taxid`, `rdfs:subClassOf`, `NA`))
)


# split pubmed IDs
if (!is.null(ids)) {
    v.id <- unlist(str_split(ids, ","))

    values <- paste0(v.id, collapse = " ")

    q = paste(sparql_prefix, paste0(
    "SELECT distinct ?chem_name ?CID ?plant_name ?taxid 
        (group_concat(distinct(?part_name); separator=\",\") as ?plant_parts)
        (group_concat(distinct(?reference); separator=\",\") as ?references)
    WHERE { 
        VALUES ?pcid { ",values," }
        ?cmp sen:PUBCHEM_CID ?pcid .
        ?cmp rdf:type sen:compound .
        ?cmp rdfs:label ?chem_name .
        ?plant ^sen:inPlant ?pcx .
        ?pcx sen:hasCompound ?cmp .
        ?plant rdfs:label ?plant_name .
        OPTIONAL { 
            ?pcx sen:hasReference ?ref . 
            ?ref rdfs:label ?reference
        }
        OPTIONAL { 
            ?prt ^sen:inComponent ?pcx .
            ?prt rdfs:label ?part_name
        }
        OPTIONAL { ?plant rdfs:subClassOf|sen:has_taxid ?taxon }
        bind( strafter(str(?pcid),str(pcid:)) as ?CID) .
        bind( strafter(str(?taxon),str(ncbi_taxon:)) as ?taxid) .
    } 
    group by ?chem_name ?CID ?plant_name ?taxid
    order by ?chem_name ?plant_name ?part_name"))

    res <- SPARQL(endpoint, q, ns=prefix, extra=query_options, format="json")$results

    write_csv(res, outFile)

}

# compoundSmiles = "alpha-Amyrin|Baicalein|Baicalin|beta-Amyrin|beta-Sitosterol-d-glucoside|Dihydrokavain|Dihydromethysticin|Honokiol|Isoliquiritigenin|L-tetrahydropalmatine|Magnolol|Methysticin|Tenuifolin|Valerenic-acid|Vincamine|Yangonin|alpha-Pinene"
if (!is.null(cmpActFile)) {
    message(paste0("reading compounds from: ", cmpActFile))
    if (grepl("tsv$", cmpActFile)) {
        df.cmpAct <- read_tsv(cmpActFile, locale = locale(encoding = "UTF-8"))
    } else {
        df.cmpAct <- read_csv(cmpActFile, locale = locale(encoding = "UTF-8"))
    }
    compoundIds = unique(df.cmpAct[,cmpIdCol]) %>% unlist
}

if (!is.null(compoundIds)) {
    # compoundSmiles = unlist(str_split(compoundSmiles, "\\|"))
    # cmpSmiless = paste0(compoundSmiles, collapse = "\" \"")

    getPlantsForCompound <- function(cmpID) {
        q = paste(sparql_prefix, paste0(
        "SELECT distinct ?cmp ?cmp_label ?cmp_smiles ?pln 
        (group_concat(distinct(?part_name); separator=\"|\") as ?plant_parts)
        (group_concat(distinct(?reference); separator=\"|\") as ?references)
        (group_concat(distinct(?location_label); separator=\"|\") as ?locations)
        (group_concat(distinct(?pcx_src_label); separator=\"|\") as ?plant_compound_source)
        (group_concat(distinct(?use_src_label); separator=\"|\") as ?plant_location_source)
        WHERE { 
            # cmp in plant
            values ?cmp { ", cmpID, " }
            ?cmp sen:lcLabel ?cmp_label .
            ?pln ^sen:inPlant ?pcx .
            ?pcx sen:hasCompound ?cmp .
            ?pcx sen:hasSource ?pcx_src .
            ?pcx_src rdfs:label ?pcx_src_label .
            # ?pln rdfs:label ?plant_name .
            OPTIONAL {
                ?cmp sen:StdIsomericSMILES ?cmp_smiles .
            }
            OPTIONAL { 
                ?pcx sen:hasReference ?ref . 
                ?ref rdfs:label ?reference
            }
            OPTIONAL { 
                ?prt ^sen:inComponent ?pcx .
                ?prt rdfs:label ?part_name
            }
            #OPTIONAL { ?pln rdfs:subClassOf|sen:has_taxid ?taxon }
            #bind( strafter(str(?taxon),str(ncbi_taxon:)) as ?taxid) .
            # plant in location
            OPTIONAL {
                ?pln ^sen:hasTaxon ?use .
                ?use sen:hasLocation ?location .
                ?location rdfs:label ?location_label .
                ?use sen:hasSource ?use_src .
                ?use_src rdfs:label ?use_src_label .
            }
        } 
        group by ?cmp ?cmp_label ?cmp_smiles ?pln
        order by ?cmp_label ?part_name"))

        SPARQL(endpoint, q, ns=prefix, extra=query_options, format="json")$results
    }
    message("getting plants for compounds")
    suppressWarnings(
        plantsForCompounds <- fix_sparql_ids(bind_rows(mclapply(compoundIds, getPlantsForCompound, mc.cores = 10))) %>%
        left_join(df.id, by = "pln")
    )

    write_csv(plantsForCompounds, outFile)

}