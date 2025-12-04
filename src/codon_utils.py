"""
Codon preference table handling utilities.
"""


def read_codon_preference_file(filename):
    """
    Read the codon preference file and return all lines.
    
    Args:
        filename (str): Path to codon preference file
        
    Returns:
        list: All lines from the codon preference file
    """
    with open(filename) as codon_file:
        return codon_file.readlines()


def get_species_codon_preference(species, codon_file_lines):
    """
    Extract codon preference table for a specific species.
    
    Args:
        species (str): Species code (e.g., 'Ec', 'Yeast', 'Human', etc.)
        codon_file_lines (list): Lines from codon preference file
        
    Returns:
        list: Species codon preference table or empty list if not found
    """
    species_str = species + ","
    species_codon_preference_table_list = []
    
    for line in codon_file_lines:
        if species_str in line:
            species_line_index = codon_file_lines.index(line)
            species_codon_preference_table_list = codon_file_lines[species_line_index+2:species_line_index+66]
            
            print(f"The species selected is {codon_file_lines[species_line_index]}")
            print("The codon preference table is as follows")
            for codon_line in species_codon_preference_table_list:
                print(codon_line.strip())
            
            break
    
    return species_codon_preference_table_list


def validate_species(species_codon_preference_list):
    """
    Validate that species was found in codon preference file.
    
    Args:
        species_codon_preference_list (list): Species codon preference table
        
    Returns:
        bool: True if species found, False otherwise
    """
    if len(species_codon_preference_list) == 0:
        print("Species not found. Please rerun with appropriate species. Thank you!")
        return False
    return True


def get_best_codon(new_amino_acid, species_codon_preference_table_list):
    """
    Get the codon with the highest preference value for a specific amino acid.
    
    Args:
        new_amino_acid (str): Target amino acid
        species_codon_preference_table_list (list): Species codon preference table
        
    Returns:
        str: Best codon line (format: 'GGG\tG\t0.21')
    """
    potential_codons_list = []
    for codon_line in species_codon_preference_table_list:
        codon_line = codon_line.strip()
        if new_amino_acid == codon_line[4]:
            potential_codons_list.append(codon_line)
    
    mutagnenic_codon_str = potential_codons_list[0]
    if len(potential_codons_list) == 1:
        print("The potential codon is as following")
        print(mutagnenic_codon_str)
    else:
        print("The potential codons are as follows")
        for i in range(len(potential_codons_list)):
            print(potential_codons_list[i])
            if (i >= 1) and (float(potential_codons_list[i][6:10]) > float(potential_codons_list[0][6:10])):
                mutagnenic_codon_str = potential_codons_list[i]
    
    print(f"The codon used is {mutagnenic_codon_str}")
    return mutagnenic_codon_str
