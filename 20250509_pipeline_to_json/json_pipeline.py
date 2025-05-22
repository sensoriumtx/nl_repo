import pandas as pd
import json
import argparse
import os
import traceback

# Enhanced JSON Builder: Full Diagnostic Debugging

def build_dynamic_json(primary_file, acts_file, output_json, scoring_file=None, assays_file=None, log_file='json_builder_log.txt'):
    try:
        print("Starting JSON builder...")
        with open(log_file, 'w') as log:
            log.write("Starting JSON builder...\n")

            # Step 1: Load Primary File (Anchor for cmp)
            print(f"Loading primary file: {primary_file}")
            if not os.path.exists(primary_file):
                print(f"Error: Primary file {primary_file} does not exist.")
                return

            try:
                df_primary = pd.read_csv(primary_file).fillna("")
                print(f"Primary file loaded with {len(df_primary)} rows.")
            except Exception as e:
                print(f"Error loading primary file: {e}")
                return

            if 'cmp' not in df_primary.columns:
                print("Error: Primary file must contain a 'cmp' column.")
                return

            # Initialize JSON structure
            df_primary['cmp'] = df_primary['cmp'].str.strip().str.lower()
            cmp_list = df_primary['cmp'].unique().tolist()
            print(f"Initialized JSON for {len(cmp_list)} compounds.")

            json_output = {cmp: {"cmp": cmp, "cmp_labels": [], "properties": {}, "activities": {"compound_activities": {}, "plant_associations": {}}, "assays": []} for cmp in cmp_list}

            # Step 2: Map CMP Labels and PLN Labels
            for _, row in df_primary.iterrows():
                cmp = row['cmp'].strip().lower()
                if cmp in json_output:
                    cmp_label = str(row.get('cmp_label', '')).strip()
                    cmp_labels = str(row.get('cmp_labels', '')).split('|') if row.get('cmp_labels', '') else []
                    json_output[cmp]['cmp_labels'] = list(set(json_output[cmp]['cmp_labels'] + [cmp_label] + cmp_labels))

            # Step 3: Map Activities First
            print(f"Loading activities file: {acts_file}")
            df_activities = pd.read_csv(acts_file).fillna("")
            print(f"Activities file loaded with {len(df_activities)} rows.")

            for _, row in df_activities.iterrows():
                cmp = str(row.get('cmp', '')).strip().lower()
                if cmp in json_output:
                    pln = str(row.get('pln', '')).strip()
                    pln_label = str(row.get('pln_label', '')).strip()
                    act_pln = str(row.get('act', '')).strip()
                    act_label_pln = str(row.get('act_label_x', '')).strip()
                    act_cmp = str(row.get('act.x', '')).strip()
                    act_label_cmp = str(row.get('act_label_y', '')).strip()

                    # Map Compound Activities
                    if act_cmp and act_label_cmp:
                        json_output[cmp]['activities']["compound_activities"].setdefault(act_cmp, []).append(act_label_cmp)
                        json_output[cmp]['activities']["compound_activities"][act_cmp] = list(set(json_output[cmp]['activities']["compound_activities"][act_cmp]))

                    # Map Plant Activities
                    if pln and pln_label:
                        json_output[cmp]['activities']["plant_associations"].setdefault(pln, {"pln_label": pln_label, "act_pln": {}})
                        if act_pln and act_label_pln:
                            json_output[cmp]['activities']["plant_associations"][pln]["act_pln"].setdefault(act_pln, []).append(act_label_pln)
                            json_output[cmp]['activities']["plant_associations"][pln]["act_pln"][act_pln] = list(set(json_output[cmp]['activities']["plant_associations"][pln]["act_pln"][act_pln]))

            print("Finished mapping activities.")

            # Step 4: Map Scoring
            # Step 4: Map Scoring (Minimized Fields)
            if scoring_file:
                print(f"Loading scoring file: {scoring_file}")
                df_scoring = pd.read_csv(scoring_file).fillna("")
                print(f"Scoring file loaded with {len(df_scoring)} rows.")

                scoring_columns = [
                    "cmp","label","inchi","isoSmiles","pcid","casid",
                    "primary_SMILES","sensBBB","SCScore","ro5_violations",
                    "is_PAINS","is_macrocycle","MolWt","HeavyAtomMolWt","TPSA",
                    "HeavyAtomCount","NumAliphaticCarbocycles","NumAliphaticHeterocycles",
                    "NumAliphaticRings","NumAromaticCarbocycles","NumAromaticHeterocycles",
                    "NumAromaticRings","NumHAcceptors","NumHDonors","NumHeteroatoms","NumRotatableBonds",
                    "NumSaturatedCarbocycles","NumSaturatedHeterocycles","NumSaturatedRings","RingCount","MolLogP",
                    "QED","Priority","PriorityRationale"
                ]

                for _, row in df_scoring.iterrows():
                    cmp = str(row.get('cmp', '')).strip().lower()
                    if cmp in json_output:
                        for col in scoring_columns:
                            value = str(row.get(col, '')).strip()
                            if value:
                                json_output[cmp]['properties'][col] = value

            print("Finished mapping scoring properties.")

            # Step 5: Map Assays
            if assays_file:
                print(f"Loading assays file: {assays_file}")
                df_assays = pd.read_csv(assays_file).fillna("")
                print(f"Assays file loaded with {len(df_assays)} rows.")

                for _, row in df_assays.iterrows():
                    cmp = str(row.get('cmp', '')).strip().lower()
                    if cmp in json_output:
                        assay_data = {col: str(row.get(col, '')).strip() for col in row.index if col != 'cmp'}
                        json_output[cmp]['assays'].append(assay_data)

            print("Finished mapping assays.")

            # Save JSON
            print("Saving JSON...")
            temp_json = output_json.replace('.json', '_TEMP.json')
            with open(temp_json, 'w') as json_file:
                json.dump(json_output, json_file, indent=4)

            print(f"Temporary JSON saved to {temp_json}.")

    except Exception as e:
        print(f"Critical Error encountered: {e}")
        traceback.print_exc()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Dynamic JSON Builder - Full Diagnostic Debugging")
    parser.add_argument('--primary', type=str, required=True, help='Primary CMP file')
    parser.add_argument('--acts', type=str, required=True, help='Activities file')
    parser.add_argument('--scoring', type=str, required=False, help='Scoring file')
    parser.add_argument('--assays', type=str, required=False, help='Assays file')
    parser.add_argument('--out', type=str, required=True, help='Output JSON file')

    args = parser.parse_args()

    build_dynamic_json(args.primary, args.acts, args.out, args.scoring, args.assays)
