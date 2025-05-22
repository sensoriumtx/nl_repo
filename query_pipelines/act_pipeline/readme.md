````
---

```markdown
# Compound Classification Pipeline

This R-based pipeline automates the retrieval, enrichment, and classification of plant-derived compounds based on therapeutic activities. It uses SPARQL-based queries against a knowledge graph, merges and deduplicates the results, and classifies compounds using structural information.

---

## Directory Structure

```

.
├── scripts/
│   ├── pull\_act\_for\_plant\_enrich.r
│   ├── pull\_cmp\_for\_pln.r
│   ├── pull\_act\_for\_cmp\_enrich.r
│   └── pull\_class\_for\_cmp.r
├── ethnobotany_pipeline.r         
├── readme.md
└── output/
├── step1\_plants.csv
├── step2\_cmp.csv
├── step3\_cmp\_act.csv
├── final\_merged\_output.csv
├── deliverable.csv
├── final\_deliverable.csv
└── pipeline\_log.txt

````

---
````
## Pipeline Steps

### **Step 1: Fetch Plants for Activity**
Queries plants directly associated with specified activity terms (`act_pln`).

### **Step 2: Fetch Compounds for Plants**
Finds all compounds from those plants using `pln_label` as input.

### **Step 3: Fetch Activities for Compounds**
Searches for compounds directly associated with the same activities (`act_cmp`). Not all compounds will match.

### **Step 4: Merge**
Combines plant-compound and compound-activity data. Also appends unmatched but active compounds with no associated plant.

### **Step 5: Final Deliverable**
Aggregates all plant and compound relationships per compound. Optionally merges scoring data.

### **Step 6: Compound Classification**
Runs a classifier on the compound list and appends predictions to the final deliverable.
````
---
````
## Usage

```bash
Rscript ethnobotany_pipeline.r \
  --endpoint dev \
  --acts "obesity|inflammation|neuroprotection" \
  --filter_out_act /path/to/filter.csv \
  --scoring /path/to/scoring.csv \
  --outdir /path/to/output/
```

You can also provide a CSV of terms with a `terms` column instead of `--acts`:

```bash
Rscript ethnobotany_pipeline.r \
  --endpoint dev \
  --act_file /path/to/activities.csv \
  --filter_out_act /path/to/filter.csv \
  --scoring /path/to/scoring.csv \
  --outdir /path/to/output/
```

> Note: You must provide either `--acts` or `--act_file`, but not both. The script will exit if neither is present.

---

## Output Files

| File                      | Description                                    |
| ------------------------- | ---------------------------------------------- |
| `step1_plants.csv`        | Plants matched to the queried activities       |
| `step2_cmp.csv`           | Compounds matched to those plants              |
| `step3_cmp_act.csv`       | Compounds matched to the queried activities    |
| `final_merged_output.csv` | Combined raw dataframe of compounds and plants |
| `deliverable.csv`         | Aggregated summary by compound with scoring    |
| `final_deliverable.csv`   | Final summary with structural classification   |
| `pipeline_log.txt`        | Full execution log                             |
````
---
````
## Notes

* All returned plants are confirmed to be associated with the specified activity terms.
* Some compounds may not map to any plant (and vice versa) but are retained if they match the activity query.
* The pipeline gracefully handles partial data and unmapped entries.

---

## Requirements

* R ≥ 4.0.0
* R packages: `tidyverse`, `lubridate`, `parallel`
* Scripts must be accessible in the `scripts/` folder
* A containerized version exists with all dependencies pre-installed

---

## Design Principles

* **Traceability**: Intermediate files are preserved and labeled.
* **Scalability**: Uses parallel processing (10 threads).
* **Validation**: Ensures all inputs are aligned and labeled properly.

---

## Contact

For issues, questions, or suggestions, please contact:
`nick.laskowski@sensorium.bio`

```
````