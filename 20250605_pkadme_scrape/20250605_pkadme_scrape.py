import pandas as pd
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, AllChem
from rdkit import DataStructs
import requests
import json
import time
from typing import Dict, List
import logging
from tdc.single_pred import ADME

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def read_smiles_csv(file_path: str) -> pd.DataFrame:
    """Read CSV file containing isomeric SMILES."""
    try:
        df = pd.read_csv(file_path)
        if 'isomeric_smiles' not in df.columns:
            raise ValueError("CSV must contain 'isomeric_smiles' column")
        return df
    except Exception as e:
        logger.error(f"Error reading CSV file: {e}")
        raise

def validate_smiles(smiles: str) -> bool:
    """Validate SMILES string using RDKit."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        return mol is not None
    except:
        return False

def canonicalize_smiles(smiles: str) -> str:
    """Convert SMILES to canonical form."""
    try:
        return Chem.CanonSmiles(smiles)
    except:
        return ""

def get_pubchem_data(smiles: str) -> Dict:
    """Fetch physicochemical and ADME-related properties from PubChem."""
    try:
        compounds = pcp.get_compounds(smiles, 'smiles')
        if not compounds:
            logger.warning(f"No PubChem data found for SMILES: {smiles}")
            return {}
       
        compound = compounds[0]
        props = pcp.get_properties([
            'MolecularWeight', 'XLogP', 'HBondDonorCount', 'HBondAcceptorCount',
            'TPSA', 'RotatableBondCount'
        ], 'smiles', smiles)
       
        return {
            'pubchem_cid': compound.cid,
            'molecular_weight': compound.molecular_weight,
            'logp': props[0].get('XLogP') if props else None,
            'h_bond_donors': props[0].get('HBondDonorCount') if props else None,
            'h_bond_acceptors': props[0].get('HBondAcceptorCount') if props else None,
            'tpsa': props[0].get('TPSA') if props else None,
            'rotatable_bonds': props[0].get('RotatableBondCount') if props else None
        }
    except Exception as e:
        logger.error(f"Error fetching PubChem data for SMILES {smiles}: {e}")
        return {}

def get_rdkit_data(smiles: str) -> Dict:
    """Compute ADME-related properties using RDKit (proxy for SwissADME)."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return {}
       
        return {
            'logp_rdkit': Descriptors.MolLogP(mol),
            'bioavailability_score': Descriptors.FractionCSP3(mol),
            'molecular_weight_rdkit': Descriptors.MolWt(mol),
            'lipinski_violations': sum([
                Descriptors.MolWt(mol) > 500,
                Descriptors.MolLogP(mol) > 5,
                Descriptors.NumHDonors(mol) > 5,
                Descriptors.NumHAcceptors(mol) > 10
            ]),
            'synthetic_accessibility': rdMolDescriptors.CalcSAscore(mol)
        }
    except Exception as e:
        logger.error(f"Error computing RDKit data for SMILES {smiles}: {e}")
        return {}

def get_chembl_data(smiles: str, target_chembl_id: str = None) -> Dict:
    """Fetch ADME and bioactivity data from ChEMBL."""
    try:
        canonical_smiles = canonicalize_smiles(smiles)
        url = f"https://www.ebi.ac.uk/chembl/api/data/activity.json?canonical_smiles={canonical_smiles}"
        if target_chembl_id:
            url += f"&target_chembl_id={target_chembl_id}"
        response = requests.get(url, timeout=10)
        if response.status_code != 200:
            logger.warning(f"No ChEMBL data for SMILES: {smiles}")
            return {}
       
        data = response.json()
        activities = data.get('activities', [])
        if not activities:
            return {}
       
        return {
            'chembl_id': activities[0].get('molecule_chembl_id'),
            'ic50': activities[0].get('standard_value') if activities[0].get('standard_type') == 'IC50' else None,
            'oral_bioavailability': activities[0].get('ro3_pass') == 'Y',
            'target_name': activities[0].get('target_pref_name')
        }
    except Exception as e:
        logger.error(f"Error fetching ChEMBL data for SMILES {smiles}: {e}")
        return {}

def get_pharmabench_data(smiles: str, pharmabench_df: pd.DataFrame) -> Dict:
    """Fetch ADMET data from PharmaBench dataset."""
    try:
        canonical_smiles = canonicalize_smiles(smiles)
        if not canonical_smiles:
            return {}
        match = pharmabench_df[pharmabench_df['smiles'] == canonical_smiles]
        if match.empty:
            return {}
        return {
            'caco2_permeability': match['caco2_permeability'].iloc[0] if 'caco2_permeability' in match.columns else None,
            'clearance_hepatocyte': match['clearance_hepatocyte'].iloc[0] if 'clearance_hepatocyte' in match.columns else None,
            'hia': match['hia'].iloc[0] if 'hia' in match.columns else None
        }
    except Exception as e:
        logger.error(f"Error fetching PharmaBench data for SMILES {smiles}: {e}")
        return {}

def get_tdc_data(smiles: str, tdc_df: pd.DataFrame) -> Dict:
    """Fetch ADMET data from Therapeutics Data Commons."""
    try:
        canonical_smiles = canonicalize_smiles(smiles)
        if not canonical_smiles:
            return {}
        match = tdc_df[tdc_df['Drug'] == canonical_smiles]
        if match.empty:
            return {}
        return {
            'solubility_aqsol': match['Y'].iloc[0] if 'Y' in match.columns else None
        }
    except Exception as e:
        logger.error(f"Error fetching TDC data for SMILES {smiles}: {e}")
        return {}

def get_drugbank_data(smiles: str, drugbank_df: pd.DataFrame) -> Dict:
    """Fetch PK data from DrugBank open data subset."""
    try:
        canonical_smiles = canonicalize_smiles(smiles)
        if not canonical_smiles:
            return {}
        match = drugbank_df[drugbank_df['smiles'] == canonical_smiles]
        if match.empty:
            return {}
        return {
            'drugbank_id': match['drugbank_id'].iloc[0] if 'drugbank_id' in match.columns else None,
            'bioavailability': match['bioavailability'].iloc[0] if 'bioavailability' in match.columns else None,
            'half_life': match['half_life'].iloc[0] if 'half_life' in match.columns else None
        }
    except Exception as e:
        logger.error(f"Error fetching DrugBank data for SMILES {smiles}: {e}")
        return {}

def get_bindingdb_data(smiles: str, bindingdb_df: pd.DataFrame) -> Dict:
    """Fetch binding affinity data from BindingDB."""
    try:
        canonical_smiles = canonicalize_smiles(smiles)
        if not canonical_smiles:
            return {}
        match = bindingdb_df[bindingdb_df['smiles'] == canonical_smiles]
        if match.empty:
            return {}
        return {
            'bindingdb_id': match['bindingdb_id'].iloc[0] if 'bindingdb_id' in match.columns else None,
            'ic50_target': match['ic50'].iloc[0] if 'ic50' in match.columns else None
        }
    except Exception as e:
        logger.error(f"Error fetching BindingDB data for SMILES {smiles}: {e}")
        return {}

def get_tox21_data(smiles: str, tox21_df: pd.DataFrame) -> Dict:
    """Fetch toxicity data from Tox21 dataset."""
    try:
        canonical_smiles = canonicalize_smiles(smiles)
        if not canonical_smiles:
            return {}
        match = tox21_df[tox21_df['smiles'] == canonical_smiles]
        if match.empty:
            return {}
        return {
            'tox21_herg': match['herg_activity'].iloc[0] if 'herg_activity' in match.columns else None,
            'tox21_ames': match['ames_activity'].iloc[0] if 'ames_activity' in match.columns else None
        }
    except Exception as e:
        logger.error(f"Error fetching Tox21 data for SMILES {smiles}: {e}")
        return {}

def get_chemical_diversity(smiles: str, reference_smiles: List[str]) -> Dict:
    """Compute chemical diversity (Tanimoto similarity)."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2)
        similarities = []
        for ref_smiles in reference_smiles:
            ref_mol = Chem.MolFromSmiles(ref_smiles)
            if ref_mol:
                ref_fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(ref_mol, 2)
                similarities.append(DataStructs.TanimotoSimilarity(fp, ref_fp))
        return {
            'tanimoto_similarity_mean': sum(similarities) / len(similarities) if similarities else None
        }
    except Exception as e:
        logger.error(f"Error computing chemical diversity for SMILES {smiles}: {e}")
        return {}

def process_compound(smiles: str, pharmabench_df: pd.DataFrame, tdc_df: pd.DataFrame, drugbank_df: pd.DataFrame,
                    bindingdb_df: pd.DataFrame, tox21_df: pd.DataFrame, reference_smiles: List[str],
                    target_chembl_id: str = None) -> Dict:
    """Process a single SMILES string and fetch PK/ADME and additional data."""
    if not validate_smiles(smiles):
        logger.warning(f"Invalid SMILES: {smiles}")
        return {'isomeric_smiles': smiles, 'valid': False}
   
    result = {'isomeric_smiles': smiles, 'valid': True}
    result.update(get_pubchem_data(smiles))
    result.update(get_rdkit_data(smiles))
    result.update(get_chembl_data(smiles, target_chembl_id))
    result.update(get_pharmabench_data(smiles, pharmabench_df))
    result.update(get_tdc_data(smiles, tdc_df))
    result.update(get_drugbank_data(smiles, drugbank_df))
    result.update(get_bindingdb_data(smiles, bindingdb_df))
    result.update(get_tox21_data(smiles, tox21_df))
    result.update(get_chemical_diversity(smiles, reference_smiles))
    return result

def main(file_path: str, output_path: str, pharmabench_path: str, drugbank_path: str,
         bindingdb_path: str, tox21_path: str, reference_smiles_path: str, target_chembl_id: str = None):
    """Main pipeline to process SMILES and extract PK/ADME and additional data."""
    # Read input CSV
    df = read_smiles_csv(file_path)
    logger.info(f"Loaded {len(df)} compounds from {file_path}")
   
    # Load datasets
    try:
        pharmabench_df = pd.read_csv(pharmabench_path)
        logger.info(f"Loaded PharmaBench dataset from {pharmabench_path}")
    except Exception as e:
        logger.error(f"Error loading PharmaBench dataset: {e}")
        pharmabench_df = pd.DataFrame()
   
    try:
        drugbank_df = pd.read_csv(drugbank_path)
        logger.info(f"Loaded DrugBank dataset from {drugbank_path}")
    except Exception as e:
        logger.error(f"Error loading DrugBank dataset: {e}")
        drugbank_df = pd.DataFrame()
   
    try:
        bindingdb_df = pd.read_csv(bindingdb_path)
        logger.info(f"Loaded BindingDB dataset from {bindingdb_path}")
    except Exception as e:
        logger.error(f"Error loading BindingDB dataset: {e}")
        bindingdb_df = pd.DataFrame()
   
    try:
        tox21_df = pd.read_csv(tox21_path)
        logger.info(f"Loaded Tox21 dataset from {tox21_path}")
    except Exception as e:
        logger.error(f"Error loading Tox21 dataset: {e}")
        tox21_df = pd.DataFrame()
   
    try:
        tdc_data = ADME(name='Solubility_AqSolDB')
        tdc_df = tdc_data.get_data()
        logger.info("Loaded TDC Solubility_AqSolDB dataset")
    except Exception as e:
        logger.error(f"Error loading TDC dataset: {e}")
        tdc_df = pd.DataFrame()
   
    try:
        ref_smiles_df = pd.read_csv(reference_smiles_path)
        reference_smiles = ref_smiles_df['smiles'].tolist()
        logger.info(f"Loaded {len(reference_smiles)} reference SMILES")
    except Exception as e:
        logger.error(f"Error loading reference SMILES: {e}")
        reference_smiles = []
   
    # Process each SMILES
    results = []
    for idx, row in df.iterrows():
        smiles = row['isomeric_smiles']
        logger.info(f"Processing compound {idx + 1}/{len(df)}: {smiles}")
        result = process_compound(smiles, pharmabench_df, tdc_df, drugbank_df, bindingdb_df, tox21_df,
                                 reference_smiles, target_chembl_id)
        results.append(result)
        time.sleep(0.5)  # Avoid API rate limits
   
    # Convert results to DataFrame and save
    output_df = pd.DataFrame(results)
    output_df.to_csv(output_path, index=False)
    logger.info(f"Saved results to {output_path}")

if __name__ == "__main__":
    input_file = "smiles.csv"
    output_file = "pkadme_extended_data.csv"
    pharmabench_file = "pharmabench_admet.csv"
    drugbank_file = "drugbank_open.csv"
    bindingdb_file = "bindingdb_data.csv"
    tox21_file = "tox21_data.csv"
    reference_smiles_file = "reference_smiles.csv"
    target_chembl_id = None  # Set to specific target (e.g., 'CHEMBL4282' for EGFR)
    main(input_file, output_file, pharmabench_file, drugbank_file, bindingdb_file, tox21_file,
         reference_smiles_file, target_chembl_id)