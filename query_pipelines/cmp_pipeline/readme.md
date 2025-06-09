# Semantic Association Pipelines (`cmp_pipeline.r`, `smiles_pipeline.r`, `wildcard_cmp.r`)

## Overview

These pipelines support semantic discovery of plant-compound-activity relationships by starting from different identifiers:

* **`cmp_pipeline.r`** – Starts with known compound IDs (`cmp`)
* **`smiles_pipeline.r`** – Starts with chemical structures (SMILES)
* **`wildcard_cmp.r`** – Starts with wildcard queries (partial strings) to match `cmp` records

Each script extracts structured biological relationships from a Sensorium knowledge graph, performs enrichment steps, and outputs tabular results that link compounds to plants and activities.

---

## General Requirements

* R ≥ 4.1
* Packages: `tidyverse`, `lubridate`, `parallel`
* All scripts must be run with `Rscript` from the command line
* External dependency scripts must exist in the `scripts/` folder:

  * `pull_acts_for_specific_cmp_ids.r`
  * `pull_plant_for_compound_ids.r`
  * `pull_acts_for_specific_pln.r`
  * (For SMILES only) `pull_act_for_smiles.r`

---

## `cmp_pipeline.r`

### Description

Starts with one or more compound identifiers (e.g., `sen:SENCMP000000000006`) and performs:

1. Retrieval of activities linked to the compounds
2. Mapping of plants containing those compounds
3. Enrichment of activities linked to those plants

### Arguments

| Argument     | Type   | Description                                  |
| ------------ | ------ | -------------------------------------------- |
| `--cmp`      | string | One or more pipe-delimited `cmp` IDs         |
| `--endpoint` | string | Target SPARQL endpoint (`dev`, `prod`, etc.) |
| `--outdir`   | string | Output directory path                        |

### Use Cases

#### Single `cmp` string

```bash
Rscript nl_repo/query_pipelines/cmp_pipeline/cmp_pipeline.r \
  --endpoint dev \
  --cmp "sen:SENCMP000000000006" \
  --outdir /sensorium-research-kb/dev/data/query_output/testing/for_nick/20250604_cmp_pipeline_single_cmp_test
```

#### Multiple `cmp` strings

```bash
Rscript nl_repo/query_pipelines/cmp_pipeline/cmp_pipeline.r \
  --endpoint dev \
  --cmp "sen:SENCMP000000000006|sen:SENCMP000000000033|sen:SENCMP000000000046" \
  --outdir /sensorium-research-kb/dev/data/query_output/testing/for_nick/20250604_cmp_pipeline_multi_cmp_test
```

---

## `smiles_pipeline.r`

### Description

Starts from one or more SMILES strings and performs:

1. SMILES → Compound resolution + direct activity mapping
2. Compound → Plant resolution
3. Plant → Activity enrichment

### Arguments

| Argument     | Type   | Description                                  |
| ------------ | ------ | -------------------------------------------- |
| `--smiles`   | string | Pipe-delimited SMILES strings                |
| `--endpoint` | string | Target SPARQL endpoint (`dev`, `prod`, etc.) |
| `--outdir`   | string | Output directory path                        |

### Use Cases

#### Single SMILES string

```bash
Rscript nl_repo/query_pipelines/cmp_pipeline/smiles_pipeline.r \
  --endpoint dev \
  --smiles "COc1ccc2c(c1)O[C@H]1c3ccc(O)cc3OC[C@@H]21" \
  --outdir /sensorium-research-kb/dev/data/query_output/testing/for_nick/20250604_smile_pipeline_test
```

#### Multiple SMILES strings

```bash
Rscript nl_repo/query_pipelines/cmp_pipeline/smiles_pipeline.r \
  --endpoint dev \
  --smiles "SMILES1||SMILES2" \
  --outdir /sensorium-research-kb/dev/data/query_output/testing/for_nick/20250604_smile_pipeline_test
```

---

## `wildcard_cmp.r`

### Description

Starts with a partial string (`search`) and performs:

1. Case-insensitive fuzzy search across all string fields in a master compound file
2. Extracts matched `cmp` identifiers
3. Proceeds with compound → plant → activity mapping (like `cmp_pipeline.r`)

### Arguments

| Argument     | Type   | Description                                                       |
| ------------ | ------ | ----------------------------------------------------------------- |
| `--search`   | string | Pipe-delimited partial string(s) to search across fields          |
| `--in_file`  | string | Master CSV file containing compound metadata (must include `cmp`) |
| `--endpoint` | string | Target SPARQL endpoint (`dev`, `prod`, etc.)                      |
| `--outdir`   | string | Output directory path                                             |

### Use Cases

#### Single Wildcard

```bash
Rscript nl_repo/query_pipelines/cmp_pipeline/wildcard_cmp.r \
  --endpoint dev \
  --search "1S/C21H21NO6/..." \
  --in_file /sensorium-research-kb/dev/data/query_output/activity/20250603_master_cmp_w_np.csv \
  --outdir /sensorium-research-kb/dev/data/query_output/testing/for_nick/20250604_wildcard_single_test
```

#### Multiple Wildcards

```bash
Rscript nl_repo/query_pipelines/cmp_pipeline/bulk_wildcard.r \
  --endpoint dev \
  --search "kavain|1S/C15H10O5/..." \
  --in_file /sensorium-research-kb/dev/data/query_output/activity/20250603_master_cmp_w_np.csv \
  --outdir /sensorium-research-kb/dev/data/query_output/testing/for_nick/20250604_wildcard_multi_test
```

#### Wildcard With Blank Compound Activities

Even if no compound activities are found, the pipeline continues to query plant and plant-associated activities.

```bash
Rscript nl_repo/query_pipelines/cmp_pipeline/bulk_wildcard.r \
  --endpoint dev \
  --search "1s/c10h16/c..." \
  --in_file /sensorium-research-kb/dev/data/query_output/activity/20250603_master_cmp_w_np.csv \
  --outdir /sensorium-research-kb/dev/data/query_output/testing/for_nick/20250604_wildcard_single_test_only_cmpact_blanks
```

---

## Output Files

All pipelines create a log and CSV outputs such as:

* `step0_acts_for_cmp.csv` / `step0_acts_for_smiles.csv`
* `step1_plants_for_cmp.csv`
* `step2_acts_for_pln.csv`
* (For wildcard: `step0_wildcard_matched_rows.csv`)

These outputs represent progressive mappings across:

* **Compound → Activity**
* **Compound → Plant**
* **Plant → Activity**

---

## Notes

* All scripts fail gracefully with logging when partial or blank data is encountered.
* Output folder must be specified (`--outdir`) or the pipeline will terminate early.
* Scripts are modular and reusable in batch or workflow automation.

---

Let me know if you'd like me to export this to separate README files or integrate it into a shared documentation structure.
