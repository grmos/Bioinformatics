# Protein Sequence Analysis

## Overview

This project provides scripts for extracting protein sequences from a text file, creating a FASTA file, building a BLAST database, and performing a grid search using BLASTP. The scripts utilize the Uniprot database for protein information and the NCBI BLAST tools for sequence analysis.

### Prerequisites

Before running the scripts, make sure to install the required dependencies:

# Enviroment
- Python 3.x
- Install libraries: `pip3 install -r requirements.txt`
- NCBI BLAST (version 2.15.0+)

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
