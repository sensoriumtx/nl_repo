# `act_pipeline.r`

**Plant-Compound-Activity Semantic Enrichment Pipeline**

## Overview

This R script (`act_pipeline_no_class_query.r`) performs a multi-stage semantic enrichment pipeline that starts from input activities (e.g., “obesity”, “anxiety”) and identifies:

1. **Plants** associated with those activities
2. **Compounds** found in those plants
3. **Activities** associated with those compounds

It then merges all outputs into a final structured dataset, optionally integrating compound-level scoring data.

The pipeline enables end-to-end exploration of bioactive plant compounds through knowledge graph queries, including filtering, enrichment, and merging.

---

## Features

* Accepts activity terms via direct string input or CSV file
* Queries knowledge graph endpoints in parallel
* Filters negative or irrelevant activities (optional)
* Produces 3 core outputs:

  * Plants ↔ Activities
  * Plants ↔ Compounds
  * Compounds ↔ Activities
* Final merged file with optional scoring overlay
* Parallel processing for high performance

---

## Arguments

| Argument           | Type   | Description                                                               |             |
| ------------------ | ------ | ------------------------------------------------------------------------- | ----------- |
| `--endpoint`       | string | Target SPARQL endpoint (e.g., `dev`, `prod`)                              |             |
| `--acts`           | string | Pipe-delimited activity terms (e.g., \`"obesity                           | anxiety"\`) |
| `--act_file`       | string | CSV file containing one column named `term` or `terms` with activities    |             |
| `--filter_out_act` | string | Path to a CSV of activity terms to exclude (e.g., flagged negative terms) |             |
| `--outdir`         | string | Output directory path where results will be saved                         |             |
| `--scoring`        | string | (Optional) CSV file with compound-level scores to merge into final output |             |

**Note:** Either `--acts` or `--act_file` is required.

---

## Output Files

| File                                    | Description                                         |
| --------------------------------------- | --------------------------------------------------- |
| `step1_plants.csv`                      | Plants associated with activity terms               |
| `step2_cmp.csv`                         | Compounds found in those plants                     |
| `step3_chunks/step3_cmp_enrichment.csv` | Activities associated with those compounds          |
| `final_merged_output.csv`               | Merged mapping of plants, compounds, and activities |
| `final_deliverable.csv`                 | Aggregated, scored compound-plant-activity mapping  |
| `pipeline_log.txt`                      | Log file for execution steps                        |

---

## Use Cases

### 1. **Single Activity Term**

```bash
Rscript nl_repo/query_pipelines/act_pipeline/act_pipeline_no_class_query.r \
  --endpoint dev \
  --acts "obesity" \
  --filter_out_act /sensorium-research-kb/dev/data/filter/master_act_snomed_mapping_flag_negative_properties.csv \
  --scoring /sensorium-research-kb/dev/data/query_output/activity/20250603_master_cmp_w_np.csv \
  --outdir /sensorium-research-kb/dev/data/query_output/testing/for_nick/20250605_act_pipeline_single_string_test
```

### 2. **Multiple Activity Terms**

```bash
Rscript nl_repo/query_pipelines/act_pipeline/act_pipeline_no_class_query.r \
  --endpoint dev \
  --acts "obesity|anxiety" \
  --filter_out_act /sensorium-research-kb/dev/data/filter/master_act_snomed_mapping_flag_negative_properties.csv \
  --scoring /sensorium-research-kb/dev/data/query_output/activity/20250603_master_cmp_w_np.csv \
  --outdir /sensorium-research-kb/dev/data/query_output/testing/for_nick/20250605_act_pipeline_multi_string_test
```

### 3. **Activity Terms from File**

```bash
Rscript nl_repo/query_pipelines/act_pipeline/act_pipeline_no_class_query.r \
  --endpoint dev \
  --act_file /sensorium-research-kb/dev/data/query_output/testing/for_nick/act_files/20250605_act_file.csv \
  --filter_out_act /sensorium-research-kb/dev/data/filter/master_vct_snomed_mapping_flag_negative_properties.csv \
  --scoring /sensorium-research-kb/dev/data/query_output/activity/20250603_master_cmp_w_np.csv \
  --outdir /sensorium-research-kb/dev/data/query_output/testing/for_nick/20250605_act_pipeline_act_file_test1
```

---

## Performance Notes

* This pipeline is **bandwidth and computation heavy**.
* Enrichment queries span multiple relationships in the knowledge graph.
* Expect longer runtimes for large sets of input activity terms.

---

## Developer Notes

* Parallelization is controlled via `parallel::makeCluster(workers)`
* Subscripts used:

  * `scripts/pull_act_for_plant_enrich.r`
  * `scripts/pull_cmp_for_pln.r`
  * `scripts/pull_act_for_cmp_enrich.r`

Ensure these helper scripts exist and are accessible from the script’s root directory.

---

## Contact

For issues or extensions, please contact the Sensorium Knowledge Graph team.

---

Let me know if you’d like a `README.md` file export or want this added to a documentation repo.
