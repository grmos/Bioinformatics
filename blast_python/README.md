# Protein Sequence Analysis

## Overview

This project provides Python scripts for extracting protein sequences from a text file containing protein names. The main functionality involves querying the Uniprot database using the `requests` library, creating a FASTA file with protein sequences, building a BLAST database using the NCBI BLAST tools, and performing a grid search with BLASTP.

### Prerequisites

Before running the scripts, make sure to install the required dependencies:

- Python 3.x
- `requests` library
- `tqdm` library
- NCBI BLAST (version 2.15.0+)

```bash
pip install requests tqdm

Download and install NCBI BLAST from [https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

## Steps

### Step 1: Create a FASTA file with protein sequences

Run the following command to generate a FASTA file from a text file containing protein names:

```bash
python Data_generator.py -rf ../proteins.txt -wf ../proteins.fasta

## Step 2: Create the BLAST Database

Run the following command to create a BLAST database from the generated FASTA file:

```bash
makeblastdb -in ../proteins.fasta -dbtype prot -out ../database

## Step 3: Run Grid Search using BLASTP

Execute the following command to perform a grid search with BLASTP:

```bash
python grid_search.py -qf ../Q13148.fasta -df database/databasename -of ../results/ -bf ../ncbi-blast-2.15.0+/bin/blastp -e 1.5

**Adjust the paths and parameters according to your specific setup.**

## Acknowledgments

- NCBI BLAST: [https://blast.ncbi.nlm.nih.gov/](https://blast.ncbi.nlm.nih.gov/)
- Uniprot: [https://www.uniprot.org/](https://www.uniprot.org/)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.