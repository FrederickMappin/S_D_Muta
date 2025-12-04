#usage:
#python3 main.py DNAsequence.txt mutaa.txt species [output_directory]
#
#species can be one of the following Ec, Yeast, Insect,
#Ce, Dm, Human, Mouse, Rat, Pig, Pp, At, Streptomyces,
#Zm, Tobacco, Sc, Cg.

"""
This script searches possible mutagenic primers for overlap extension PCR.
The resulting overlap primers are sorted according to the theoretical Tm value of
the "overlap" sequence during the overlapping stage. Also, the length, GC content
and Tm value of each primer are provided to help designing the gene primers and
set annealing temperature for the PCR run. We hope these listed primers are helpful
for preparing mutagenic primers with overlap extension PCR.
"""

import sys
import os

# Add parent directory to path to import from src
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from dna_utils import read_and_filter_dna_sequence, print_dna_sequence, validate_gene_length
from codon_utils import read_codon_preference_file, get_species_codon_preference, validate_species, get_best_codon
from mutation_utils import read_amino_acid_mutations, parse_amino_acid_mutation, validate_mutation_position, print_mutation_info
from primer_utils import generate_tetrad_primers, sort_primers, print_primer_summary
from output_utils import write_primers_to_file, organize_primer_files

# Get current path
current_path = os.getcwd()
print(f"Current path is {current_path}.")
print("\n")

# Get output directory if provided
output_dir = None
if len(sys.argv) > 4:
    output_dir = sys.argv[4]
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    print(f"Output directory: {output_dir}")
    print("\n")

# Read and validate gene sequence
gene_sequence_str = read_and_filter_dna_sequence(sys.argv[1])
print_dna_sequence(gene_sequence_str)

if not validate_gene_length(gene_sequence_str):
    sys.exit()

# Read and process codon preference file
species = sys.argv[3]
codon_file_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'codon_preferences', 'CodonPreference.txt')
codon_file_lines = read_codon_preference_file(codon_file_path)
species_codon_preference_table_list = get_species_codon_preference(species, codon_file_lines)

if not validate_species(species_codon_preference_table_list):
    sys.exit()

print("\n")

# Process each amino acid mutation
amino_acid_mutations = read_amino_acid_mutations(sys.argv[2])

for amino_acid_mutation in amino_acid_mutations:
    # Parse mutation
    original_amino_acid, amino_acid_position, new_amino_acid = parse_amino_acid_mutation(amino_acid_mutation)
    
    # Validate mutation position
    if not validate_mutation_position(amino_acid_position, len(gene_sequence_str)):
        print(f"Amino acid mutation file contains a line with mutation position within the first or last 65.")
        print(amino_acid_mutation)
        print("Mutagenic primers can be directly designed without overlap extension PCR.\n")
        print("Similar message may occur if there is still other such positions.\n")
        print("Please delete them all from the txt file of amino acid mutation and run script again. Thank you!\n")
        sys.exit()
    
    # Print mutation info and get best codon
    print_mutation_info(original_amino_acid, amino_acid_position, new_amino_acid)
    mutagenic_codon_str = get_best_codon(new_amino_acid, species_codon_preference_table_list)
    mutagenic_codon_used_str = mutagenic_codon_str[0:3]
    
    print("\n")
    
    # Generate and sort primers
    primer_trimers_list = generate_tetrad_primers(gene_sequence_str, amino_acid_position, mutagenic_codon_used_str)
    primer_trimers_list_sorted = sort_primers(primer_trimers_list)
    
    # Print primer summary
    print_primer_summary(primer_trimers_list, primer_trimers_list_sorted)
    print("\n")
    
    # Write results to file
    write_primers_to_file(amino_acid_mutation, primer_trimers_list_sorted, output_dir)

# Organize primer files
gene_file = sys.argv[1]
organize_primer_files(gene_file, current_path, output_dir)
