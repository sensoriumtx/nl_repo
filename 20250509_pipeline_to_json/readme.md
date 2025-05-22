Usage

--primary # file contains all cmp/pln associations
--acts # file contains all activities for both cmp and pln can be the same file as primary
--scoring # pulls all scoring metrics mapped to cmp within a file, can be the master file or other
--assays # boassay fetch mapped to cmp
--out # .json format output



Example:

python json_pipeline.py \
--primary 20250509_primary.csv \
--acts 20250415_act_pln_cmp_act_merge.csv \
--scoring 20250509_primary_w_scoring.csv \
--assays 20250509_bioassay_od.csv \
--out 20250509_od_35_cmp.json

Starting JSON builder...
Loading primary file: 20250509_primary.csv
Primary file loaded with 4944 rows.
Initialized JSON for 35 compounds.
Loading activities file: 20250415_act_pln_cmp_act_merge.csv
Activities file loaded with 168535 rows.
Finished mapping activities.
Loading scoring file: 20250509_primary_w_scoring.csv
Scoring file loaded with 4944 rows.
Finished mapping scoring properties.
Loading assays file: 20250509_bioassay_od.csv
Assays file loaded with 2893 rows.
Finished mapping assays.
Saving JSON...