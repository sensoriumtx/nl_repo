import pandas as pd
from chembl_webresource_client.new_client import new_client
from tqdm import tqdm
import argparse
import os

def get_targets_for_chembl_id(chembl_id):
    try:
        activities = new_client.activity.filter(molecule_chembl_id=chembl_id)
        targets = []
        for act in activities:
            if act.get('target_chembl_id'):
                targets.append({
                    'chembl_id': chembl_id,
                    'assay_chembl_id': act.get('assay_chembl_id'),
                    'target_chembl_id': act.get('target_chembl_id'),
                    'target_type': act.get('target_type'),
                    'activity_type': act.get('standard_type'),
                    'activity_value': act.get('standard_value'),
                    'activity_units': act.get('standard_units')
                })
        return targets
    except Exception as e:
        print(f"[WARN] Failed to fetch activities for {chembl_id}: {e}")
        return []

def get_gene_symbol_for_target(target_id):
    try:
        target = new_client.target.filter(target_chembl_id=target_id).only("target_components")[0]
        components = target.get("target_components", [])
        for comp in components:
            for xref in comp.get("target_component_xrefs", []):
                if xref.get("xref_src_db") == "UniProt":
                    return comp.get("accession")
        return None
    except Exception as e:
        print(f"[WARN] Failed to fetch target for {target_id}: {e}")
        return None

def enrich_targets_with_genes(target_data):
    enriched = []
    for entry in tqdm(target_data):
        gene_symbol = get_gene_symbol_for_target(entry["target_chembl_id"])
        entry["gene_symbol"] = gene_symbol
        enriched.append(entry)
    return enriched

def run_pipeline(input_csv, output_csv):
    df = pd.read_csv(input_csv)
    if 'chembl_id' not in df.columns:
        raise ValueError("Input CSV must contain a column named 'chembl_id'.")

    all_results = []
    for chembl_id in tqdm(df['chembl_id'].dropna().unique()):
        target_data = get_targets_for_chembl_id(chembl_id)
        all_results.extend(target_data)

    enriched_results = enrich_targets_with_genes(all_results)
    out_df = pd.DataFrame(enriched_results)
    out_df.to_csv(output_csv, index=False)
    print(f"[DONE] Output saved to {output_csv}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Map ChEMBL compound IDs to gene networks.")
    parser.add_argument('--in', required=True, help="Path to input CSV with column 'chembl_id'")
    parser.add_argument('--out', required=True, help="Path to output CSV")

    args = parser.parse_args()
    run_pipeline(args.in, args.out)
