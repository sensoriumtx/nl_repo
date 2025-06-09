# Plant-Based Semantic Association Pipelines (`pln_pipeline.r`, `batch_pln_pipeline.r`, `pln_wildcard.r`)

## Overview

These R-based semantic pipelines begin with one or more **plant names** (PLN values) and resolve their associated:

1. **Compounds found in the plant**
2. **Activities associated with those compounds**
3. **Activities directly associated with the plant**

Each script supports a different input or performance modality:

* `pln_pipeline.r`: Lightweight, single-shot pipeline from a user-defined plant string.
* `batch_pln_pipeline.r`: High-performance, chunked execution for large batches.
* `pln_wildcard.r`: Partial string search across plant metadata followed by semantic enrichment.

---

## Shared Requirements

* R ≥ 4.1
* Packages: `tidyverse`, `parallel`, `lubridate`
* Required helper scripts (in `scripts/` directory):

  * `pull_cmp_for_pln.r`
  * `pull_acts_for_specific_cmp_ids.r`
  * `pull_acts_for_specific_pln.r`

---

## `pln_pipeline.r`

### Description

Starts with a direct list of plant names or labels. For each:

* Resolves compounds found in the plant
* Pulls activities associated with those compounds
* Pulls activities directly associated with the plant

### Arguments

| Argument     | Type   | Description                                       |              |
| ------------ | ------ | ------------------------------------------------- | ------------ |
| `--plants`   | string | Pipe-delimited plant names (e.g., \`Panax ginseng | Aloe vera\`) |
| `--endpoint` | string | SPARQL endpoint (`dev`, `prod`)                   |              |
| `--outdir`   | string | Output directory path                             |              |

### Example Use

```bash
Rscript nl_repo/query_pipelines/pln_pipeline/pln_pipeline.r \
  --endpoint dev \
  --plants "Panax ginseng|Withania somnifera" \
  --outdir /sensorium-research-kb/dev/data/query_output/testing/for_nick/20250604_pln_pipeline_test
```

---

## `batch_pln_pipeline.r`

### Description

Optimized for large-scale plant input. This pipeline:

* Chunks plant names into batches of 300
* Queries compounds, plant-level activities, and compound-level activities in parallel
* Merges all chunked results into final output files

### Arguments

| Argument     | Type   | Description                                         |
| ------------ | ------ | --------------------------------------------------- |
| `--plants`   | string | Pipe-delimited plant names or labels                |
| `--endpoint` | string | SPARQL endpoint (`dev`, `prod`, etc.)               |
| `--outdir`   | string | Output directory for all merged and chunked outputs |

### Example Use

```bash
Rscript nl_repo/query_pipelines/pln_pipeline/batch_pln_pipeline.r \
  --endpoint dev \
  --plants "Panax ginseng|Withania somnifera|Aloe vera|Curcuma longa" \
  --outdir /sensorium-research-kb/dev/data/query_output/testing/for_nick/20250604_batch_pln_pipeline_test
```

---

## `pln_wildcard.r`

### Description

Starts from **partial string queries** (e.g., "ginseng", "ashwagandha") and:

1. Performs a fuzzy case-insensitive search over a plant metadata CSV
2. Extracts plant labels from matches
3. Resolves compounds and all associated activities

### Arguments

| Argument     | Type   | Description                                                     |           |
| ------------ | ------ | --------------------------------------------------------------- | --------- |
| `--search`   | string | Pipe-delimited search terms for partial match (e.g., \`"ginseng | ashwa"\`) |
| `--in_file`  | string | Master plant metadata CSV (must include `pln_label` column)     |           |
| `--endpoint` | string | SPARQL endpoint (`dev`, `prod`, etc.)                           |           |
| `--outdir`   | string | Output directory path                                           |           |

### Example Use

#### Standard wildcard query

```bash
Rscript nl_repo/query_pipelines/pln_pipeline/pln_wildcard.r \
  --endpoint dev \
  --search "ginseng|ashwagandha" \
  --in_file /sensorium-research-kb/dev/data/query_output/master_pln_file.csv \
  --outdir /sensorium-research-kb/dev/data/query_output/testing/for_nick/20250604_wildcard_pln_test
```

#### Graceful handling of unmatched inputs

Even if the plant or compound associations are missing, the script continues and writes valid (but empty) output files to preserve pipeline continuity.

---

## Output Files

Each script creates a consistent set of outputs in the specified `--outdir`:

| File                                | Description                                                |
| ----------------------------------- | ---------------------------------------------------------- |
| `step1_compounds_for_pln.csv`       | Plant → Compound relationships                             |
| `step2_acts_for_cmp.csv`            | Compound → Activity mappings                               |
| `step3_acts_for_pln.csv`            | Plant → Activity mappings                                  |
| `step0_wildcard_matched_plants.csv` | (Only in `pln_wildcard.r`) Match results from fuzzy search |

---

## Summary

| Script                 | Input Type                | Use Case                             |
| ---------------------- | ------------------------- | ------------------------------------ |
| `pln_pipeline.r`       | Direct plant string       | Lightweight direct queries           |
| `batch_pln_pipeline.r` | Large-scale plant strings | High-performance parallel execution  |
| `pln_wildcard.r`       | Fuzzy search              | Partial matches over metadata fields |

---

Please log any issues and bugs found when running this script and email bug description to nick.laskowski@sensorium.bio