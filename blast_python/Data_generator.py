import argparse
import requests
from tqdm import tqdm 

ap = argparse.ArgumentParser()
ap.add_argument("-rf", "--read_filename", type=str,  required=True,help="filename of the proteins IDs file that is read (.txt)")
ap.add_argument("-wf", "--write_filename", type=str,  required=True,help="filename of the proteins IDs file that is written (.fasta)")

args = vars(ap.parse_args())

def read_protein_ids_from_file(file_path):
    protein_ids = []
    with open(file_path, 'r') as file:
        for line in file:
            protein_ids.append(line.strip())
    return protein_ids

def fetch_protein_sequence(uniprot_id):
    url = f'https://www.uniprot.org/uniprot/{uniprot_id}.fasta'
    response = requests.get(url)
    
    if response.ok:
        return response.text
    else:
        print(f"Failed to fetch sequence for {uniprot_id}. Status code: {response.status_code}")
        return None



file_path =args ["read_filename"] 
uniprot_ids = read_protein_ids_from_file(file_path)

sequences=[]
for uniprot_id in tqdm(uniprot_ids, desc="Fetching sequences", unit=" protein"):
    sequence = fetch_protein_sequence(uniprot_id)
    
    if sequence:
        sequences.append(f"\n{sequence}")

# Save all sequences to a single file
with open(args ["write_filename"] , "w") as file:
    file.write("\n".join(sequences))



