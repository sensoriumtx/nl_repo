import pandas as pd
from chembl_webresource_client.new_client import new_client
import logging
from typing import Optional, List, Dict
import sys
from rdkit import Chem
from tenacity import retry, stop_after_attempt, wait_fixed
from concurrent.futures import ThreadPoolExecutor
import threading

# Logging setup
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(threadName)s - %(message)s',
    handlers=[
        logging.FileHandler('chembl_pipeline.log'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

def validate_smiles(smiles: str) -> Optional[str]:
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return Chem.MolToSmiles(mol, canonical=True)
        else:
            logger.warning(f"Invalid SMILES: {smiles}")
            return None
    except Exception as e:
        logger.warning(f"SMILES validation error ({smiles}): {str(e)}")
        return None

def clean_smiles_list(smiles_str: str) -> List[str]:
    return [s.strip() for s in smiles_str.split(';') if s.strip()]

def read_csv(file_path: str) -> pd.DataFrame:
    try:
        df = pd.read_csv(file_path)
        if 'primary_SMILES' not in df.columns or 'cmp' not in df.columns:
            raise ValueError("CSV must contain 'primary_SMILES' and 'cmp' columns")

        logger.info(f"Loaded {len(df)} rows from CSV")

        df['primary_SMILES_list'] = df['primary_SMILES'].apply(
            lambda x: clean_smiles_list(x) if isinstance(x, str) else []
        )

        df['primary_SMILES_canonical'] = df['primary_SMILES_list'].apply(
            lambda lst: [s for s in (validate_smiles(sm) for sm in lst) if s]
        )

        if 'InChI' in df.columns:
            df['InChI_list'] = df['InChI'].apply(lambda x: x.split(';') if isinstance(x, str) else [])
        else:
            df['InChI_list'] = [[] for _ in range(len(df))]

        logger.info(f"Total validated SMILES: {df['primary_SMILES_canonical'].apply(len).sum()}")
        return df

    except Exception as e:
        logger.error(f"Error loading CSV: {str(e)}")
        raise

@retry(stop=stop_after_attempt(3), wait=wait_fixed(2))
def get_chembl_id_from_smiles(smiles: str, molecule_client) -> Optional[Dict]:
    try:
        result = molecule_client.filter(molecule_structures__canonical_smiles=smiles).only(['molecule_chembl_id'])
        if result:
            return {'chembl_id': result[0]['molecule_chembl_id'], 'identifier': smiles, 'identifier_type': 'SMILES'}
        else:
            logger.warning(f"No ChEMBL match for SMILES: {smiles}")
            return None
    except Exception as e:
        logger.error(f"ChEMBL query error for SMILES {smiles}: {str(e)}")
        return None

@retry(stop=stop_after_attempt(3), wait=wait_fixed(2))
def get_target_details(target_id: str, target_client) -> Dict:
    try:
        target = target_client.filter(target_chembl_id=target_id).only(
            ['target_chembl_id', 'pref_name', 'organism', 'target_type']
        )
        if target:
            return {
                'target_chembl_id': target[0]['target_chembl_id'],
                'target_name': target[0]['pref_name'] or 'Unknown',
                'target_organism': target[0]['organism'] or 'Unknown',
                'target_type': target[0]['target_type'] or 'Unknown'
            }
        return {'target_chembl_id': target_id, 'target_name': 'Unknown', 'target_organism': 'Unknown', 'target_type': 'Unknown'}
    except Exception as e:
        logger.error(f"Target lookup error for {target_id}: {str(e)}")
        return {'target_chembl_id': target_id, 'target_name': 'Unknown', 'target_organism': 'Unknown', 'target_type': 'Unknown'}

def get_assay_data(chembl_id: str, activity_client, target_client) -> List[Dict]:
    try:
        activities = activity_client.filter(molecule_chembl_id=chembl_id).only(
            ['assay_chembl_id', 'assay_type', 'assay_description', 'standard_type', 
             'standard_value', 'standard_units', 'target_chembl_id']
        )
        assay_data = []
        for act in activities:
            target_info = get_target_details(act['target_chembl_id'], target_client) if act['target_chembl_id'] else {}

            if target_info.get('target_organism') not in ['Homo sapiens', 'Human']:
                continue
            if target_info.get('target_type') != 'SINGLE PROTEIN':
                continue

            assay_data.append({
                'assay_chembl_id': act.get('assay_chembl_id'),
                'assay_type': act.get('assay_type', 'Unknown'),
                'assay_description': act.get('assay_description', 'N/A'),
                'standard_type': act.get('standard_type', 'N/A'),
                'standard_value': act.get('standard_value'),
                'standard_units': act.get('standard_units', 'N/A'),
                'target_chembl_id': target_info.get('target_chembl_id'),
                'target_name': target_info.get('target_name'),
                'target_organism': target_info.get('target_organism'),
                'target_type': target_info.get('target_type')
            })
        return assay_data
    except Exception as e:
        logger.error(f"Error retrieving assays for {chembl_id}: {str(e)}")
        return []

def process_single_cmp(row: pd.Series) -> tuple[str, Dict]:
    cmp = row['cmp']
    smiles_list = row['primary_SMILES_canonical']
    logger.info(f"Processing {cmp}: {len(smiles_list)} SMILES")

    molecule_client = new_client.molecule
    activity_client = new_client.activity
    target_client = new_client.target

    chembl_ids = []
    for smiles in smiles_list:
        res = get_chembl_id_from_smiles(smiles, molecule_client)
        if res:
            chembl_ids.append(res)

    assays = []
    for chembl_info in chembl_ids:
        chembl_id = chembl_info['chembl_id']
        bioactivity = get_assay_data(chembl_id, activity_client, target_client)
        for a in bioactivity:
            a.update({'identifier': chembl_info['identifier'], 'identifier_type': chembl_info['identifier_type']})
        assays.extend(bioactivity)

    return cmp, {'chembl_ids': chembl_ids, 'assays': assays}

def process_smiles_list(df: pd.DataFrame, max_workers: int) -> Dict[str, Dict]:
    results = {}
    with ThreadPoolExecutor(max_workers=max_workers, thread_name_prefix="ChEMBLWorker") as executor:
        future_map = {executor.submit(process_single_cmp, row): row['cmp'] for _, row in df.iterrows()}
        for future in future_map:
            cmp = future_map[future]
            try:
                cmp, result = future.result()
                results[cmp] = result
            except Exception as e:
                logger.error(f"Error for cmp {cmp}: {str(e)}")
                results[cmp] = {'chembl_ids': [], 'assays': []}
    return results

def explode_assays(df: pd.DataFrame, results: Dict[str, Dict]) -> pd.DataFrame:
    """
    Create a DataFrame where each assay is in its own row, mapping back to compound, SMILES, and ChEMBL ID.
    """
    exploded_rows = []

    for _, row in df.iterrows():
        cmp = row['cmp']
        assays = results.get(cmp, {}).get('assays', [])
        chembl_ids = results.get(cmp, {}).get('chembl_ids', [])

        if assays:
            for assay in assays:
                matched_smiles = assay.get('identifier')
                matched_chembl = next(
                    (cid['chembl_id'] for cid in chembl_ids if cid['identifier'] == matched_smiles),
                    None
                )

                new_row = row.drop(['primary_SMILES_list', 'primary_SMILES_canonical', 'InChI_list']).to_dict()
                new_row.update({
                    'cmp': cmp,
                    'ChEMBL_ID': matched_chembl,
                    'identifier': matched_smiles,
                    'identifier_type': assay.get('identifier_type'),
                    'target_chembl_id': assay.get('target_chembl_id'),
                    'target_name': assay.get('target_name'),
                    'target_organism': assay.get('target_organism'),
                    'target_type': assay.get('target_type'),
                    'assay_chembl_id': assay.get('assay_chembl_id'),
                    'assay_type': assay.get('assay_type'),
                    'assay_description': assay.get('assay_description'),
                    'standard_type': assay.get('standard_type'),
                    'standard_value': assay.get('standard_value'),
                    'standard_units': assay.get('standard_units')
                })
                exploded_rows.append(new_row)
        else:
            # No assay data found, still record for completeness
            new_row = row.drop(['primary_SMILES_list', 'primary_SMILES_canonical', 'InChI_list']).to_dict()
            new_row.update({
                'cmp': cmp,
                'ChEMBL_ID': None,
                'identifier': None,
                'identifier_type': None,
                'target_chembl_id': None,
                'target_name': None,
                'target_organism': None,
                'target_type': None,
                'assay_chembl_id': None,
                'assay_type': None,
                'assay_description': None,
                'standard_type': None,
                'standard_value': None,
                'standard_units': None
            })
            exploded_rows.append(new_row)

    return pd.DataFrame(exploded_rows)


def main(input_path: str, output_path: str, max_workers: int):
    try:
        df = read_csv(input_path)
        results = process_smiles_list(df, max_workers)
        output_df = explode_assays(df, results)
        output_df.to_csv(output_path, index=False)
        logger.info(f"Results written to {output_path}")
        logger.info(f"{output_df['ChEMBL_ID'].notna().sum()} ChEMBL IDs, "
                    f"{output_df['target_chembl_id'].notna().sum()} targets, "
                    f"{output_df['assay_chembl_id'].notna().sum()} assays found.")
    except Exception as e:
        logger.error(f"Pipeline failed: {str(e)}")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Extract ChEMBL data from SMILES")
    parser.add_argument("--input", required=True, help="Path to input CSV")
    parser.add_argument("--output", default="chembl_pipeline_output.csv", help="Output CSV path")
    parser.add_argument("--workers", type=int, default=10, help="Parallel worker count")
    args = parser.parse_args()
    main(args.input, args.output, args.workers)
