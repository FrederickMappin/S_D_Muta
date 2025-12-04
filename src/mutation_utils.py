"""
Amino acid mutation parsing and validation utilities.
"""


def read_amino_acid_mutations(filename):
    """
    Read amino acid mutations from a file.
    
    Args:
        filename (str): Path to amino acid mutations file
        
    Returns:
        list: List of amino acid mutation strings
    """
    with open(filename) as amino_acid_mutation_file:
        amino_acid_mutation_file_lines = amino_acid_mutation_file.readlines()
    
    return [line.strip() for line in amino_acid_mutation_file_lines]


def parse_amino_acid_mutation(amino_acid_mutation):
    """
    Parse a single amino acid mutation string.
    
    Format example: "K123E" means replace lysine (K) at position 123 with glutamic acid (E)
    
    Args:
        amino_acid_mutation (str): Mutation string
        
    Returns:
        tuple: (original_aa, position, new_aa)
    """
    amino_acid_mutation_line_length = len(amino_acid_mutation)
    original_amino_acid = amino_acid_mutation[0]
    new_amino_acid = amino_acid_mutation[amino_acid_mutation_line_length - 1]
    amino_acid_position = int(str(amino_acid_mutation[1:amino_acid_mutation_line_length - 1]))
    
    return original_amino_acid, amino_acid_position, new_amino_acid


def validate_mutation_position(amino_acid_position, gene_length):
    """
    Validate that mutation position is not too close to gene ends.
    
    Requirements:
    1) amino_acid_position < 66
    2) amino_acid_position > len(gene)/3-65
    3) len(gene) > 400
    
    Args:
        amino_acid_position (int): Position of mutation in amino acid sequence
        gene_length (int): Length of gene sequence in bases
        
    Returns:
        bool: True if position is valid for overlap extension PCR, False otherwise
    """
    if amino_acid_position < 66 or amino_acid_position > int(gene_length / 3 - 65):
        return False
    return True


def print_mutation_info(original_aa, position, new_aa):
    """
    Print information about the amino acid mutation.
    
    Args:
        original_aa (str): Original amino acid
        position (int): Position of mutation
        new_aa (str): New amino acid
    """
    print(f"The orginal {original_aa} in position {position} will be changed to {new_aa}.")
