# lineage_lit_spider.py

import os
import re
import requests
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from scholarly import scholarly
from pyvis.network import Network
from crossref.restful import Works
from PyPDF2 import PdfReader

# Create folder for downloads
os.makedirs('downloads', exist_ok=True)

# Crossref API
works = Works()

# Extract text from PDF
def extract_text_from_pdf(pdf_path):
    reader = PdfReader(pdf_path)
    pages = []
    for page in reader.pages:
        pages.append(page.extract_text())
    return pages

# Extract entities (simple noun phrases for demo)
def extract_entities(text):
    entities = re.findall(r'\b([A-Z][a-z]+(?: [A-Z][a-z]+)*)\b', text)
    return list(set(entities))

# Extract references (DOIs)
def extract_references(text):
    dois = re.findall(r'10\.\d{4,9}/[-._;()/:A-Z0-9]+', text, re.IGNORECASE)
    return list(set(dois))

# Download paper from Unpaywall
def download_pdf(doi):
    url = f"https://api.unpaywall.org/v2/{doi}?email=your_email@example.com"
    r = requests.get(url).json()
    if 'best_oa_location' in r and r['best_oa_location']:
        pdf_url = r['best_oa_location']['url_for_pdf']
        if pdf_url:
            response = requests.get(pdf_url)
            if response.status_code == 200:
                path = f'downloads/{doi.replace("/", "_")}.pdf'
                with open(path, 'wb') as f:
                    f.write(response.content)
                return path
    return None

# Recursive spidering
def spider_paper(doi, layer, source_doi, all_data, visited_dois):
    if doi in visited_dois:
        return
    visited_dois.add(doi)
    metadata = works.doi(doi)
    title = metadata['title'][0] if metadata else 'Unknown'
    year = metadata['issued']['date-parts'][0][0] if metadata else 'Unknown'

    pdf_path = download_pdf(doi)
    evidence = []
    entities = []

    if pdf_path:
        pages = extract_text_from_pdf(pdf_path)
        for i, page in enumerate(pages):
            evidence.append({'Layer': layer, 'DOI': doi, 'Title': title, 'Year': year, 
                             'Source_DOI': source_doi, 'Evidence_Statement': page, 
                             'Entity': ', '.join(extract_entities(page)), 'Page_Number': i+1})

        references = []
        for page in pages[-3:]:
            references.extend(extract_references(page))
        references = list(set(references))

        for ref in references:
            spider_paper(ref, layer + 1, doi, all_data, visited_dois)

    if evidence:
        all_data.extend(evidence)

# Main Execution
start_doi = "10.3389/fphar.2023.1240295"  # Example starting DOI
all_data = []
visited_dois = set()

spider_paper(start_doi, layer=1, source_doi="ROOT", all_data=all_data, visited_dois=visited_dois)

# Save CSV
output_df = pd.DataFrame(all_data)
output_df.to_csv('literature_lineage_output.csv', index=False)

# Build Network Graph
G = nx.DiGraph()

for _, row in output_df.iterrows():
    G.add_node(row['DOI'], title=row['Title'], year=row['Year'])
    G.add_node(row['Entity'], type='entity')
    G.add_edge(row['DOI'], row['Entity'], label='MENTIONS')
    if row['Source_DOI'] != 'ROOT':
        G.add_edge(row['Source_DOI'], row['DOI'], label='CITES')

nt = Network(notebook=False, directed=True)
nt.from_nx(G)
nt.show('literature_lineage_network.html')
