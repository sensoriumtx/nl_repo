# Use Case Examples

## 1. cmp_pipeline.r use cases

### Single cmp string

```
Rscript nl_repo/query_pipelines/cmp_pipeline/cmp_pipeline.r --endpoint dev --cmp "sen:SENCMP000000000006" --outdir /sensorium-research-kb/dev/data/query_output/testing/for_nick/20250604_cmp_pipeline_single_cmp_test
```
### Multiple cmp string

```
Rscript nl_repo/query_pipelines/cmp_pipeline/cmp_pipeline.r --endpoint dev --cmp "sen:SENCMP000000000006|sen:SENCMP000000000033|sen:SENCMP000000000046" --outdir /sensorium-research-kb/dev/data/query_output/testing/for_nick/20250604_cmp_pipeline_multi_cmp_test
```

## 2. smiles_pipeline.r use cases

### Single smiles string
```
Rscript nl_repo/query_pipelines/cmp_pipeline/smiles_pipeline.r --endpoint dev --smiles "COc1ccc2c(c1)O[C@H]1c3ccc(O)cc3OC[C@@H]21" --outdir /sensorium-research-kb/dev/data/query_output/testing/for_nick/20250604_smile_pipeline_test
```

### Multiple smiles strings
```
Rscript nl_repo/query_pipelines/cmp_pipeline/smiles_pipeline.r --endpoint dev --smiles "COc1ccc2c(c1)O[C@H]1c3ccc(O)cc3OC[C@@H]21||COc1ccc2c(c1OC)C(=O)O[C@@H]2[C@H]1c2c(cc3c(c2OC)OCO3)CCN1C" --outdir /sensorium-research-kb/dev/data/query_output/testing/for_nick/20250604_smile_pipeline_test
```

## 3. wildcard_cmp.r use cases

### Single wildcard

```
Rscript nl_repo/query_pipelines/cmp_pipeline/wildcard_cmp.r --endpoint dev --search "1S/C21H21NO6/c1-22-7-6-11-8-15-16(27-10-26-15)9-13(11)18(22)19-12-4-5-14(24-2)20(25-3)17(12)21(23)28-19/h4-5,8-9,18-19H,6-7,10H2,1-3H3/t18-,19-/m1/s1" --in_file /sensorium-research-kb/dev/data/query_output/activity/20250603_master_cmp_w_np.csv --outdir /sensorium-research-kb/dev/data/query_output/testing/for_nick/20250604_wildcard_single_test  
```
### Multiple wildcard

```
Rscript nl_repo/query_pipelines/cmp_pipeline/bulk_wildcard.r --endpoint dev --search "kavain|1S/C15H10O5/c1-6-2-8-12(10(17)3-6)15(20)13-9(14(8)19)4-7(16)5-11(13)18/h2-5,16-18H,1H3" --in_file /sensorium-research-kb/dev/data/query_output/activity/20250603_master_cmp_w_np.csv --outdir /sensorium-research-kb/dev/data/query_output/testing/for_nick/20250604_wildcard_multi_test 
```

### Wildcard_ with blanks
    if there is blank outputs for compound activity associations the script will error and continue to find associations with plants and plant activities

```
Rscript nl_repo/query_pipelines/cmp_pipeline/bulk_wildcard.r --endpoint dev --search "1s/c10h16/c1-7-6-8-4-5-9(7)10(8,2)3/h8-9h,1,4-6h2,2-3h3/t8-,9+/m0/s1 " --in_file /sensorium-research-kb/dev/data/query_output/activity/20250603_master_cmp_w_np.csv --outdir /sensorium-research-kb/dev/data/query_output/testing/for_nick/20250604_wildcard_single_test_only_cmpact_blanks
```

## 4. pln_pipeline.r use cases

### pln single string

```
Rscript nl_repo/query_pipelines/pln_pipeline/pln_pipeline.r --endpoint dev --plants "Papaver somniferum" --outdir /sensorium-research-kb/dev/data/query_output/testing/for_nick/20250605_pln_pipeline_1_string_test
```

### pln multi string

```
Rscript nl_repo/query_pipelines/pln_pipeline/pln_pipeline.r --endpoint dev --plants "Papaver somniferum|galium aparine" --outdir /sensorium-research-kb/dev/data/query_output/testing/for_nick/20250605_pln_pipeline_2_string_test
```

### batch_pln_pipeline.r

    This script will intake a massive list of plant labels and chunk them accross each step to fall below dropout ceilings. (best if needed for >100 specific pln_labels)
```
Rscript nl_repo/query_pipelines/pln_pipeline/batch_pln_pipeline.r --endpoint dev --plants "galium aparine|ginko baloba|piper methysticum|Papaver somniferum" --outdir /sensorium-research-kb/dev/data/query_output/testing/for_nick/20250606_batch_pln_test1
```

### Single Wildcard

```
Rscript nl_repo/query_pipelines/pln_pipeline/pln_wildcard.r --endpoint dev --search "galium aparine" --in_file /sensorium-research-kb/dev/data/query_output/activity/20250606_master_pln_dev.csv --outdir /sensorium-research-kb/dev/data/query_output/testing/for_nick/20250606_pln_wildcard_string_input
```

### Multiple Wildcard
        if there is blank outputs for compound activity associations the script will error and continue to find associations with plants and plant activities

```
Rscript nl_repo/query_pipelines/pln_pipeline/pln_wildcard.r --endpoint dev --search "galium aparine|piper methysticum" --in_file /sensorium-research-kb/dev/data/query_output/activity/20250606_master_pln_dev.csv --outdir /sensorium-research-kb/dev/data/query_output/testing/for_nick/20250606_pln_wildcard_string_input
```

## 5. act_pipeline.r

##### This query is bandwidth heavy and will take a bit of time to compile. The enrichment queries pull the full graph in about 5 different ways and there are 2 enrichment steps.

### act single string use case

```
Rscript nl_repo/query_pipelines/act_pipeline/act_pipeline_no_class_query.r --endpoint dev --acts "obesity" --filter_out_act /sensorium-research-kb/dev/data/filter/master_act_snomed_mapping_flag_negative_properties.csv --scoring /sensorium-research-kb/dev/data/query_output/activity/20250603_master_cmp_w_np.csv --outdir /sensorium-research-kb/dev/data/query_output/testing/for_nick/20250605_act_pipeline_single_string_test
```

### act multi string use case

```
Rscript nl_repo/query_pipelines/act_pipeline/act_pipeline_no_class_query.r --endpoint dev --acts "obesity|anxiety" --filter_out_act /sensorium-research-kb/dev/data/filter/master_act_snomed_mapping_flag_negative_properties.csv --scoring /sensorium-research-kb/dev/data/query_output/activity/20250603_master_cmp_w_np.csv --outdir /sensorium-research-kb/dev/data/query_output/testing/for_nick/20250605_act_pipeline_multi_string_test
```

### act 'act_file' use case

```
Rscript nl_repo/query_pipelines/act_pipeline/act_pipeline_no_class_query.r --endpoint dev _file /sensorium-research-kb/dev/data/query_output/testing/for_nick/act_files/20250605_act_file.csv --filter_out_act /sensorium-research-kb/dev/data/filter/master_vct_snomed_mapping_flag_negative_properties.csv --scoring /sensorium-research-kb/dev/data/query_output/activity/20250603_master_cmp_w_np.csv --outdir /sensorium-research-kb/dev/data/query_output/testing/for_nick/20250605_act_pipeline_act_file_test1
```



