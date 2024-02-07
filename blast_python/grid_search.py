import subprocess
import os
import argparse


ap = argparse.ArgumentParser()
ap.add_argument("-qf", "--query_file", type=str, required=True,help="filepath of the query protein file  (.fasta)")
ap.add_argument("-df", "--database_file", type=str, required=True, help="filepath of the database ")
ap.add_argument("-of", "--output_folder", type=str, required=True, help="path of the output folder in which results.txt will be located ")
ap.add_argument("-bf", "--blastp_folder", type=str, required=True, help="path of the blastp file (../ncbi-blast-2.15.0+/bin/blastp) ")
ap.add_argument("-e", "--e_value", type=float, default=0.01, help="E-value threshold for BLAST search")

args = vars(ap.parse_args())


def run_blast(query_file, db_path, out_file, word_size, threshold, matrix,gapopen,gapextend,window_size):
    command = [
        args["blastp_folder"],
        '-query', query_file,
        '-db', db_path,
        '-word_size', str(word_size),
        '-out', out_file,
        '-evalue', str(args["e_value"]),
        '-threshold', str(threshold),
        '-matrix', matrix,
        '-gapopen',str(gapopen),
        '-gapextend',str(gapextend),
        '-window_size',str(window_size)
    ]
    subprocess.run(command, check=True)

def parse_blast_results(result_file,matrix ):
    hits = []
    ids= []
    e_values = []
    with open(result_file, 'r') as f:
        for line in f:
            if line.startswith('sp|'):
                # Assume hit if E-value < 1
                e_value_line = line.strip()
                id=e_value_line.split()[0]
                e_value = float(e_value_line.split()[-1])
                hits.append(e_value_line)
                ids.append(id)
                e_values.append(e_value)
    return hits,ids,e_values 

def main():
    query_file = args["query_file"]
    db_path = args["database_file"]
    output_folder = args["output_folder"]

    # Define parameter grids
    word_sizes =  [2,3,4,5,6,7]
    thresholds = [i for i in range(1, 30,2)]
    matrices = ['BLOSUM45','BLOSUM62','PAM250','BLOSUM80',
'BLOSUM50','BLOSUM90','PAM30','PAM70']
    w_s=[40,60,80]
    gap_p={'BLOSUM45':[(13, 3),(14, 2),(16, 1)],
'BLOSUM62':[(11, 2),(7, 2),(13, 1),(9, 1)],'PAM250':[(15, 3),(17, 2),(14, 2),(17, 1)],'BLOSUM80':[(25, 2),(9, 2),(11, 1)],'BLOSUM50':[(13, 3),(16, 2),(18, 1)],'BLOSUM90':[(9, 2),(6, 2),(9, 1)],
'PAM30':[(7, 2),(8, 1),(15, 3)],'PAM70':[(8, 2),(10, 1),(12, 3)]}
    
    hits_dict = {}
    for word_size in word_sizes:
        for threshold in thresholds:
            for matrix in matrices:
                for gapopen,gapextend in gap_p[matrix]:
                    for window_size in w_s:
                        output_file = os.path.join(output_folder, f'results_w{word_size}_t{threshold}_{matrix}_{gapopen}_{gapextend}_{window_size}.txt')

                # Run BLASTp
                        run_blast(query_file, db_path, output_file, word_size, threshold, matrix,gapopen,gapextend,window_size)

                # Parse results
                        hits,ids,e_values  = parse_blast_results(output_file,matrix )
                        if hits:
                            sw=1
                            for (iter,(id,e_v)) in enumerate(zip(ids,e_values)):
                                if(id not in hits_dict):
                                    hits_dict[id]=1000 #set flag
                                if(hits_dict[id]>e_v):
                                    sw=0
                                    hits_dict[id]=e_v
                                    if(iter<1):
                                        print(f"Hits for w={word_size}, t={threshold}, matrix={matrix},, gapopen={gapopen},, gapextend={gapextend},, window_size={window_size}:")
                                        for hit in hits:
                                            print(hit)
                                if(sw and os.path.exists(output_file)):
                                    os.remove(output_file )
                        elif( os.path.exists(output_file)):
                            os.remove(output_file )

if __name__ == "__main__":
    main()
