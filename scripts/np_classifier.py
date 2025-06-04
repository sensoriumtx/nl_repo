import os
import glob
import time
import argparse
import pandas as pd
import requests
from concurrent.futures import ThreadPoolExecutor, as_completed

CLASSIFIER_URL = "https://npclassifier.gnps2.org/classify"

def classify_smiles(smiles):
    try:
        if not smiles or pd.isna(smiles) or smiles.lower().strip() == 'nan':
            return "", "", "", smiles
        r = requests.get(CLASSIFIER_URL, params={"smiles": smiles}, timeout=10)
        r.raise_for_status()
        result = r.json()
        return (
            ', '.join(result.get('pathway_results', [])),
            ', '.join(result.get('superclass_results', [])),
            ', '.join(result.get('class_results', [])),
            smiles
        )
    except Exception:
        with open("failed_smiles.txt", "a") as f:
            f.write(smiles + "\n")
        return "", "", "", smiles

def process_chunk(chunk_df, chunk_index, outdir, max_workers):
    smiles_list = chunk_df['primary_SMILES'].fillna('').astype(str).tolist()

    results = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_smiles = {executor.submit(classify_smiles, s): s for s in smiles_list}
        for future in as_completed(future_to_smiles):
            results.append(future.result())

    pathways, superclasses, classes, smiles_out = zip(*results)
    chunk_df['primary_SMILES'] = smiles_out
    chunk_df['Pathway Results'] = pathways
    chunk_df['Superclass Results'] = superclasses
    chunk_df['Class Results'] = classes

    outfile = os.path.join(outdir, f"output_chunk_{chunk_index}.csv")
    chunk_df.to_csv(outfile, index=False)
    print(f"✅ Saved chunk {chunk_index}: {outfile}")

def merge_chunks(outdir, final_outfile):
    files = sorted(glob.glob(os.path.join(outdir, 'output_chunk_*.csv')))
    all_chunks = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    all_chunks.to_csv(final_outfile, index=False)
    print(f"🧬 Merged output written to: {final_outfile}")

def run_pipeline(input_file, outdir="classified_chunks", chunk_size=10000, max_workers=64):
    os.makedirs(outdir, exist_ok=True)
    df = pd.read_csv(input_file)
    total_chunks = (len(df) + chunk_size - 1) // chunk_size
    print(f"🚀 Classifying {len(df)} SMILES in {total_chunks} chunks using {max_workers} threads...")

    start_time = time.time()
    for i in range(total_chunks):
        start = i * chunk_size
        end = min((i + 1) * chunk_size, len(df))
        chunk_df = df.iloc[start:end].copy()
        t0 = time.time()
        process_chunk(chunk_df, i + 1, outdir, max_workers)
        print(f"⏱️ Chunk {i + 1} completed in {time.time() - t0:.2f} seconds")

    print(f"\n✅ All chunks done in {time.time() - start_time:.2f} seconds")
    merge_chunks(outdir, os.path.join(outdir, "final_npclassifier_output.csv"))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="High-speed NPClassifier SMILES classification pipeline")
    parser.add_argument("--in_file", required=True, help="Input CSV file with 'primary_SMILES' column")
    parser.add_argument("--out_dir", default="classified_chunks", help="Directory to save output chunks")
    parser.add_argument("--chunk_size", type=int, default=10000, help="Number of SMILES per chunk")
    parser.add_argument("--workers", type=int, default=64, help="Number of parallel threads")
    args = parser.parse_args()

    run_pipeline(
        input_file=args.in_file,
        outdir=args.out_dir,
        chunk_size=args.chunk_size,
        max_workers=args.workers
    )
