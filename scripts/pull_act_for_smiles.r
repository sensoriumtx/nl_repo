suppressMessages(library(tidyverse))
suppressMessages(library(stringi))
suppressMessages(library(parallel))
suppressMessages(source("scripts/common.r"))
suppressMessages(source("scripts/SPARQL.R"))

options(scipen = 100, digits = 4)
args <- commandArgs(TRUE)

# ids = "Psilocin|Mesembrine"
smiles = NULL
acts = NULL
outFile = NULL
loop = TRUE

while (loop) {
  if (args[1] == "--endpoint") {
    if (args[2] == "dev") {
      endpoint = endpoint_dev
    }
  }
  
  if (args[1] == "--smiles") {
    smiles = args[2]
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
  is.null(smiles) |
  is.null(outFile)
) {
  q("no", 1, FALSE)
}

# make out folder if doesn't exist
if (!dir.exists(dirname(outFile))) dir.create(dirname(outFile), recursive = TRUE)



# split cmps
v.smiles <- unlist(str_split(smiles, "\\|"))

# v.smiles <- c(
#     r"(COc1cc2c(cc1O)C[C@@H]1NCCc3cc(OC)c(O)c-2c31)",
#     r"(O=C(Cc1ccccc1)Oc1c(OC(=O)Cc2ccccc2)c(-c2ccc(O)cc2)c(O)c(O)c1-c1ccc(O)cc1)"
# )

outcomes = c(
  "\"Agonist\"", "\"Inhibitor\"", "\"Antagonist\"", "\"Activator\"", 
  "\"Inhibition\"", "\"Blocker\"", "\"Channel blocker\"", 
  "\"Blocker (channel blocker)\"", "\"Full agonist\"", "\"Partial agonist\"", 
  "\"Activation\"", "\"Inhibitor (gating inhibitor)\"", "\"Stimulator\"",
  "\"Inducer\"", "\"Positive\"", "\"Pore blocker\"", "\"Suppressor\"",
  "\"Potentiation\"", "\"Gating inhibitor\"", "\"Opener\"", "\"Active\"", "\"Unspecified\""
)

message("fetching compounds")
q = paste(sparql_prefix, paste0(
  "select distinct ?cmp (sample(?label) as ?cmp_label) ?cmp_smiles WHERE {
    values ?cmp_smiles { \"", paste0(v.smiles, collapse = "\" \""), "\" }
    ?cmp rdf:type sen:compound .
    ?cmp sen:StdIsomericSMILES ?cmp_smiles .
    ?cmp sen:lcLabel ?label .
    {
        ?cmp (^sen:hasCompound|(sen:maps_to+/^sen:hasCompound))/rdf:type sen:use .
    } UNION {
        values ?outcome { ", paste0(outcomes, collapse = " "), " }
        ?assayResult a sen:assayResult .
        ?assayResult sen:measureOutcome ?outcome .
        ?cmp ^sen:hasCompound ?assayResult .
        ?assayResult sen:targetGene|sen:targetProtein ?target .
        ?target ^sen:hasTarget/rdf:type sen:use .
    }
}
group by ?cmp ?cmp_smiles
"))
df.id <- SPARQL(endpoint, q, ns=prefix, extra=query_options, format='json')$results

getCmpActFromCmp <- function(SENCMP) {
  q = paste(sparql_prefix, paste0(
    "select distinct ?cmp ?isoSmiles ?inchi ?act ?act_label (group_concat(distinct ?path; separator=\"|\") as ?path) (count(distinct ?tgt) as ?target_count_for_this_compound_activity) (group_concat(distinct ?tgt; separator=\"|\") as ?target) ?assay_source where {
        select ?assay_source ?cmp ?isoSmiles ?inchi ?act ?act_label ?path ?tgt ?primary_symbol where {
        ?use a sen:use .
        values ?cmp { ", SENCMP, " } .
        {
            ?cmp (^sen:hasCompound|(sen:maps_to+/^sen:hasCompound)) ?use .
        } UNION {
            values ?outcome { ", paste0(outcomes, collapse = " "), " }
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
        OPTIONAL {
            ?cmp sen:IsomericSMILES ?isoSmiles .
        }
        OPTIONAL {
            ?cmp sen:InChI ?inchi
        }
        BIND(IF(bound(?tgt), \"indirect\", \"direct\") as ?path)
        } order by ?path ?cmp ?act
    } group by ?assay_source ?cmp ?isoSmiles ?inchi ?act ?act_label
    "))
  res <- SPARQL(endpoint, q, ns=prefix, extra=query_options, format='json')$results
  res
}
xrefResults <- fix_sparql_ids(bind_rows(mclapply(df.id$cmp, getCmpActFromCmp, mc.cores = 10))) %>%
  left_join(df.id, by = "cmp")
write_csv(xrefResults, outFile)