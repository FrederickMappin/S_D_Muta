"""
Primer generation utilities for overlap extension PCR.
"""

from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils.MeltingTemp import Tm_NN
from operator import itemgetter


def generate_tetrad_primers(gene_sequence, amino_acid_position, mutagenic_codon_used_str):
    """
    Generate tetrad primers for a specific amino acid mutation.
    
    Args:
        gene_sequence (str): Full gene sequence
        amino_acid_position (int): Position of mutation in amino acid sequence
        mutagenic_codon_used_str (str): The mutagenic codon to use
        
    Returns:
        list: List of primer trimers [overlap5_seq, overlap3_seq, Tm_overlap]
    """
    p = amino_acid_position
    primer_trimers_list = []
    
    for i in range(4):
        for j in range(4):
            for k in range(4):
                for l in range(4):
                    overlap_for_left_str = gene_sequence[3*p-12-i:3*p-3]
                    overlap_for_right_str = gene_sequence[3*p:3*p+15+j]
                    overlap5_Seq = (overlap_for_left_str +
                                    mutagenic_codon_used_str +
                                    overlap_for_right_str)

                    overlap_rev_left_str = gene_sequence[3*p-18-l:3*p-3]
                    overlap_rev_right_str = gene_sequence[3*p:3*p+9+k]
                    overlap3_Seq = Seq(overlap_rev_left_str + mutagenic_codon_used_str + overlap_rev_right_str).reverse_complement()

                    overlap_Seq = Seq(gene_sequence[3*p-12-i:3*p-3] + mutagenic_codon_used_str + gene_sequence[3*p:3*p+9+k])

                    Tm_overlap5 = round(Tm_NN(overlap5_Seq), 2)
                    Tm_overlap3 = round(Tm_NN(overlap3_Seq), 2)
                    Tm_overlap = round(Tm_NN(overlap_Seq), 2)

                    primer_trimer_list = []
                    primer_trimer_list.append(overlap5_Seq)
                    primer_trimer_list.append(overlap3_Seq)
                    primer_trimer_list.append(Tm_overlap)

                    primer_trimers_list.append(primer_trimer_list)

    return primer_trimers_list


def sort_primers(primer_trimers_list):
    """
    Sort primers by Tm value (overlap Tm) and then by sequence.
    
    Args:
        primer_trimers_list (list): List of unsorted primer trimers
        
    Returns:
        list: Sorted primer trimers list
    """
    return sorted(primer_trimers_list, key=itemgetter(2, 0))


def print_primer_summary(primer_trimers_list, primer_trimers_list_sorted):
    """
    Print summary information about original and sorted primers.
    
    Args:
        primer_trimers_list (list): Original primer list
        primer_trimers_list_sorted (list): Sorted primer list
    """
    print(f"The number of original primers is {len(primer_trimers_list)}.")
    print("The original primers are shown in the following")
    for i in range(len(primer_trimers_list)):
        for j in range(3):
            print(primer_trimers_list[i][j])
    
    print(f"\nThe number of sorted primers is {len(primer_trimers_list_sorted)}.")
    print("The sorted primers are shown in the following")
    for i in range(len(primer_trimers_list_sorted)):
        for j in range(3):
            print(primer_trimers_list_sorted[i][j])


def calculate_primer_properties(primer_sequence):
    """
    Calculate properties for a primer sequence.
    
    Args:
        primer_sequence (str or Seq): Primer sequence
        
    Returns:
        dict: Dictionary with length, GC content, and Tm value
    """
    seq_str = str(primer_sequence)
    return {
        'length': len(seq_str),
        'gc': round(gc_fraction(seq_str) * 100, 2),
        'tm': round(Tm_NN(seq_str), 2)
    }
