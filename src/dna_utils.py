"""
DNA sequence reading and validation utilities.
"""

import sys


def read_and_filter_dna_sequence(filename):
    """
    Read a DNA sequence file and filter out non-alphabetic characters.
    Skips the first line if it starts with '>' (FASTA header).
    """
    with open(filename) as dna_file:
        lines = dna_file.readlines()

    # Skip FASTA header if present
    start_idx = 1 if lines and lines[0].lstrip().startswith('>') else 0

    dna_sequence_list = []
    for line in lines[start_idx:]:
        for char in line:
            if char.isalpha():
                dna_sequence_list.append(char.upper())

    return ''.join(dna_sequence_list)


def print_dna_sequence(sequence):
    """
    Print DNA sequence in 100 character lines.
    
    Args:
        sequence (str): DNA sequence to print
    """
    print("The gene sequence is")
    for j in range(0, len(sequence), 100):
        dna_string_100_per_line = sequence[j:j+100]
        print(dna_string_100_per_line)
    print(f"The lenght of gene sequence is {len(sequence)} bps.")
    print("\n")


def validate_gene_length(gene_sequence):
    """
    Validate that the gene sequence is long enough for overlap extension PCR.
    
    Args:
        gene_sequence (str): DNA sequence to validate
        
    Returns:
        bool: True if sequence is valid, False otherwise
    """
    if len(gene_sequence) < 399:
        print("The gene sequence is less than 399 and has not enough length for the script to run.\n")
        print("In this situation, the whole plasmid PCR mutagenesis might be utilized.\n")
        print("Another way is to include the DNA sequences upstream and downstream of gene and recalculate\n")
        print("the position of mutagenic amino acids. Thank you!\n")
        return False
    return True
