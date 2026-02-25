# retarats-pubmed-pipeline

A PubMed ingestion + classification pipeline that reads a Google Sheet config (molecules + search rules), queries PubMed, and writes results into Google Sheets dashboards.

## What the script reads (required)
This script reads a Google Spreadsheet named:
* `CONFIG_SHEET_NAME` (default: `Moleculessearch`)

That Google Sheet must contain **two tabs** with these exact names:
* `MOLECULES`
* `SEARCH_RULES`

The repo includes the exact CSV exports of those two tabs:

* `MOLECULES.csv`
* `SEARCH_RULES.csv`

You should import/paste these into your own Google Sheet to get started.

## Quick start
### 1) Create your config Google Sheet

1. Go to Google Sheets and create a new spreadsheet
2. Name it: `Moleculessearch` (or any name you want — just match `CONFIG_SHEET_NAME`)
3. Create two tabs:

   * `MOLECULES`
   * `SEARCH_RULES`

### 2) Import the CSVs into the tabs

* Open `MOLECULES.csv` from this repo, copy all rows, paste into the `MOLECULES` tab starting at cell A1
* Open `SEARCH_RULES.csv` from this repo, copy all rows, paste into the `SEARCH_RULES` tab starting at cell A1

(You can also use File → Import → Upload in Google Sheets and import each CSV into the correct tab.)

### 3) Set environment variables (“bring your own keys”)

Copy `.env.example` to `.env` and fill in:

* `NCBI_EMAIL` (required)
* `NCBI_API_KEY` (recommended)
* `CONFIG_SHEET_NAME` (default: `Moleculessearch`)

### 4) Install dependencies

```bash
pip install -r requirements.txt
```

### 5) Authenticate Google access

This script uses Google Application Default Credentials (ADC).

Common approaches:

* `gcloud auth application-default login`
* OR set `GOOGLE_APPLICATION_CREDENTIALS=/path/to/service_account.json`

### 6) Run

```bash
python retarats.py
```

## Output

The pipeline creates Google Sheets in your Drive (folder: `DRIVE_FOLDER_PATH`, default `My Drive/Retarats`) including:

* `..._PEPTIDE_DATA` / `..._PEPTIDE_TABS`
* `..._SMALL_MOLECULE_DATA` / `..._SMALL_MOLECULE_TABS`
* `..._MIXTURE_DATA` / `..._MIXTURE_TABS`
  with `PAPERS_MASTER`, `STATS`, and `QUALITY_ALERTS` tabs
