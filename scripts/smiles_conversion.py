#!/usr/bin/env python3

import argparse
import pandas as pd
from rdkit import Chem
from rdkit.Chem import MolToSmiles


def is_isomeric(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid SMILES"
        canonical = Chem.MolToSmiles(mol, isomericSmiles=False, canonical=True)
        isomeric = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
        return canonical != isomeric, "Valid"
    except Exception as e:
        return False, f"Error: {str(e)}"


def standardize_smiles(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None, "Invalid SMILES"
        iso = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
        return iso, "Converted"
    except Exception as e:
        return None, f"Error: {str(e)}"


def process_csv(input_path, output_path, column):
    df = pd.read_csv(input_path)
    if column not in df.columns:
        raise ValueError(f"Column '{column}' not found in CSV.")

    df["isomeric_SMILES"] = None
    df["is_already_isomeric"] = False
    df["status"] = None

    for i, row in df.iterrows():
        smiles = row[column]
        if pd.isna(smiles) or not isinstance(smiles, str):
            df.at[i, "status"] = "Missing or invalid SMILES"
            continue

        is_isom, msg = is_isomeric(smiles)
        df.at[i, "is_already_isomeric"] = is_isom
        df.at[i, "status"] = msg

        if is_isom:
            df.at[i, "isomeric_SMILES"] = smiles
        else:
            std, status = standardize_smiles(smiles)
            df.at[i, "isomeric_SMILES"] = std
            df.at[i, "status"] = status

    df.to_csv(output_path, index=False)
    print(f"✅ Output written to {output_path}")


def process_inline(smiles_str):
    smiles_list = smiles_str.strip().split("|")
    results = []

    for s in smiles_list:
        is_isom, msg = is_isomeric(s)
        if is_isom:
            iso = s
            status = msg
        else:
            iso, status = standardize_smiles(s)

        results.append({
            "input_SMILES": s,
            "isomeric_SMILES": iso,
            "is_already_isomeric": is_isom,
            "status": status
        })

    df = pd.DataFrame(results)
    print(df.to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Standardize SMILES to isomeric format")
    parser.add_argument("--in", dest="in_file", help="Path to input CSV file")
    parser.add_argument("--out", dest="out_file", help="Path to output CSV file")
    parser.add_argument("--smiles_column", default="primary_SMILES", help="Column name in CSV")
    parser.add_argument("--smiles", help="Pipe-separated SMILES string for inline use")
    args = parser.parse_args()

    if args.smiles:
        process_inline(args.smiles)
    elif args.in_file and args.out_file:
        process_csv(args.in_file, args.out_file, args.smiles_column)
    else:
        parser.print_help()
