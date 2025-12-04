"""
Output file writing and organization utilities.
"""

import os
import shutil
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils.MeltingTemp import Tm_NN


def write_primers_to_file(amino_acid_mutation, primer_trimers_list_sorted, output_dir=None):
    """
    Write primer results to a text file.
    
    Args:
        amino_acid_mutation (str): Mutation identifier (e.g., "K123E")
        primer_trimers_list_sorted (list): Sorted list of primer trimers
        output_dir (str, optional): Directory to write files to. If None, uses current directory.
    """
    file_name_str = amino_acid_mutation + "$" + ".txt"
    if output_dir:
        file_name_str = os.path.join(output_dir, file_name_str)
    
    with open(file_name_str, 'w') as primer_file:
        # Write header
        print("primer".ljust(20, ' '),
              "sequence".ljust(40, ' '),
              "length".ljust(10, ' '),
              "GC".ljust(10, ' '),
              "Tm".ljust(10, ' '),
              "Overlap_Tm".ljust(15, ' '),
              sep="",
              file=primer_file)

        # Write primer data
        for i in range(len(primer_trimers_list_sorted)):
            for j in range(2):
                print_primer_sequence_str = primer_trimers_list_sorted[i][j]
                
                if j == 0:
                    print("\n", file=primer_file)
                    S0_print_primer_line_str = amino_acid_mutation + "_overlap_5"
                else:
                    S0_print_primer_line_str = amino_acid_mutation + "_overlap_3"

                S1_print_primer_line_str = str(print_primer_sequence_str)
                S2_print_primer_line_str = str(len(S1_print_primer_line_str))
                S3_print_primer_line_str = str(round(gc_fraction(S1_print_primer_line_str) * 100, 2))
                S4_print_primer_line_str = str(round(Tm_NN(S1_print_primer_line_str), 2))
                S5_print_primer_line_str = str(primer_trimers_list_sorted[i][2])

                print(S0_print_primer_line_str.ljust(20, ' '),
                      S1_print_primer_line_str.ljust(40, ' '),
                      S2_print_primer_line_str.ljust(10, ' '),
                      S3_print_primer_line_str.ljust(10, ' '),
                      S4_print_primer_line_str.ljust(10, ' '),
                      S5_print_primer_line_str.ljust(15, ' '),
                      sep="",
                      file=primer_file)


def organize_primer_files(gene_file, current_path, output_dir=None):
    """
    Organize primer output files into a folder.
    
    Args:
        gene_file (str): Name of the gene file (used to create folder name)
        current_path (str): Current working directory
        output_dir (str, optional): Custom output directory. If None, creates folder in current_path
    """
    if output_dir:
        # If custom output directory specified, just ensure it exists
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        # Move files to output directory
        for file in os.listdir(current_path):
            if file.endswith('$.txt'):
                shutil.move(os.path.join(current_path, file), os.path.join(output_dir, file))
        return
    
    # Default behavior: create subfolder in current_path
    folder_name = gene_file.split(".")[0] + "_primers"
    dirs = os.listdir(current_path)

    if folder_name not in dirs:
        os.system('mkdir temp_foldername')
        os.system('mv *$.txt temp_foldername')
        os.rename("temp_foldername", folder_name)
    else:
        shutil.rmtree(folder_name)
        os.system('mkdir temp_foldername')
        os.system('mv *$.txt temp_foldername')
        os.rename("temp_foldername", folder_name)
