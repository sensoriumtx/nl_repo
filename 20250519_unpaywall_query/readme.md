
````markdown
# Unpaywall PDF Retriever

This Python script queries the [Unpaywall API](https://unpaywall.org/products/api) to search for open-access scientific papers using custom search terms and downloads available PDF files. It is designed for researchers who need to programmatically retrieve and organize open-access literature.

## Features

- Queries Unpaywall using arbitrary search terms
- Automatically downloads available open-access PDFs
- Saves files using sanitized titles in a specified directory
- Simple command-line interface

## Requirements

- Python 3.6 or higher
- Required Python packages are listed in `requirements.txt`

Install them using:

```bash
pip install -r requirements.txt
````

## Directory Structure

```
.
├── unpaywall_pdf_retriever.py       # Main script
├── requirements.txt                 # Required packages
├── downloads/                       # Default directory for saved PDFs
└── README.md                        # Project documentation
```

## Usage

```bash
python unpaywall_pdf_retriever.py --quantity 5 --query "galphimine AND epilepsy" --outdir "papers"
```

### Arguments

| Argument     | Description                                                        |
| ------------ | ------------------------------------------------------------------ |
| `--quantity` | Number of papers to retrieve and attempt to download               |
| `--query`    | Search terms (e.g., `"kava AND anxiety"`)                          |
| `--outdir`   | Output directory for saving downloaded PDFs (default: `downloads`) |

## API Access

This script uses a registered email address to query the Unpaywall API:

```
nick.laskowski@sensorium.bio
```

If you plan to use this for broader or repeated use, replace the `UNPAYWALL_EMAIL` variable in the script with your own registered email.

## Output Example

```
Downloaded: Structure of Galphimine B -> papers/Structure_of_Galphimine_B.pdf
Downloaded: Conformational Study of Galphimines -> papers/Conformational_Study_of_Galphimines.pdf
```

## Limitations

* Only retrieves articles that are open-access and have an available PDF URL
* Article titles are sanitized for safe file naming
* Results depend on the coverage and accuracy of the Unpaywall index

## License

This project is distributed under the MIT License.

## Contact

**Nick Laskowski**
Email: [nick.laskowski@sensorium.bio](mailto:nick.laskowski@sensorium.bio)

```

```
