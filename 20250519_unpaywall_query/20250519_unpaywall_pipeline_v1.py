import requests
import argparse
import os
import re
import urllib.parse

# Static email for Unpaywall API access
UNPAYWALL_EMAIL = "nick.laskowski@sensorium.bio"

def query_unpaywall(quantity, query_terms, output_dir="downloads"):
    """
    Query the Unpaywall API to find papers matching specified terms and download PDFs.

    Args:
        quantity (int): Number of papers to retrieve.
        query_terms (str): Search terms (e.g., "galphimine AND obesity").
        output_dir (str): Directory to save downloaded PDFs.

    Returns:
        list: List of dictionaries containing paper titles and file paths.
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Encode query terms for URL
    encoded_query = urllib.parse.quote(query_terms)
    url = f"https://api.unpaywall.org/v2/search?query={encoded_query}&email={UNPAYWALL_EMAIL}"

    try:
        response = requests.get(url)
        response.raise_for_status()
        data = response.json()

        results = []
        papers_found = 0

        # Process results up to the requested quantity
        for item in data.get('results', []):
            if papers_found >= quantity:
                break

            title = item.get('response', {}).get('title', 'untitled')
            pdf_url = None
            for location in item.get('response', {}).get('oa_locations', []):
                if location.get('url_for_pdf'):
                    pdf_url = location['url_for_pdf']
                    break

            if not pdf_url:
                continue

            safe_title = re.sub(r'[^\w\s-]', '', title).replace(' ', '_').strip()
            if not safe_title:
                safe_title = f"paper_{papers_found + 1}"

            file_path = os.path.join(output_dir, f"{safe_title}.pdf")

            try:
                pdf_response = requests.get(pdf_url, stream=True)
                pdf_response.raise_for_status()
                with open(file_path, 'wb') as f:
                    for chunk in pdf_response.iter_content(chunk_size=8192):
                        f.write(chunk)

                results.append({'title': title, 'file_path': file_path})
                papers_found += 1
                print(f"Downloaded: {title} -> {file_path}")

            except requests.RequestException as e:
                print(f"Failed to download PDF for {title}: {e}")

        if not results:
            print("No matching papers with accessible PDFs found.")
        return results

    except requests.RequestException as e:
        print(f"Error querying Unpaywall API: {e}")
        return []
    except ValueError as e:
        print(f"Error parsing API response: {e}")
        return []

def main():
    parser = argparse.ArgumentParser(description="Query Unpaywall API and download papers.")
    parser.add_argument('--quantity', type=int, required=True, help="Number of papers to retrieve")
    parser.add_argument('--query', type=str, required=True, help="Search terms (e.g., 'compound AND activity')")
    parser.add_argument('--outdir', type=str, default="downloads", help="Directory to save downloaded PDFs")

    args = parser.parse_args()

    results = query_unpaywall(args.quantity, args.query, args.outdir)

    for result in results:
        print(f"Paper: {result['title']}, Saved to: {result['file_path']}")

if __name__ == "__main__":
    main()
