# SDM Mutagenesis - Site-Directed Mutagenesis Primer Generator

A professional Python tool for designing mutagenic primers for site-directed mutagenesis (SDM) using overlap extension PCR. Based on work of [mds-mutagenesis](https://github.com/shiqiang-lin/sdm-mutagenesis) 

## Features

- **Automatic primer design** - Generates tetrad primer combinations for SDM
- **Species-specific codons** - Uses codon preference tables for optimal expression
- **Thermal analysis** - Calculates Tm values and GC content for each primer
- **Sorted output** - Primers ranked by overlap Tm for easy selection
- **Modular architecture** - Clean separation of concerns with reusable utilities

## Quick Start

### Installation

```bash
# Clone the repository
git clone <repo>
cd sdm-mutagenesis

# Create virtual environment
python3 -m venv .venv
source .venv/bin/activate

# Install dependencies
pip install -r requirements.txt
```

### Usage

Activate the virtual environment first:
```bash
source .venv/bin/activate
```

Basic usage:
```bash
python main.py data/examples/input_gene.txt data/examples/input_mutaa.txt Human
```

With custom output directory:
```bash
python main.py data/examples/input_gene.txt data/examples/input_mutaa.txt Human output/my_results
```

Or run directly without activating:
```bash
.venv/bin/python main.py data/examples/input_gene.txt data/examples/input_mutaa.txt Human
```

### Input Files

**Gene sequence file** (`input_gene.txt`):
- DNA sequence in FASTA or plain text format
- Minimum length: 399 bps

**Mutations file** (`input_mutaa.txt`):
- One mutation per line
- Format: `[Original_AA][Position][New_AA]`
- Example: `R147A`, `K123E`, `T862R`

### Supported Species

`Ec`, `Yeast`, `Insect`, `Ce`, `Dm`, `Human`, `Mouse`, `Rat`, `Pig`, `Pp`, `At`, `Streptomyces`, `Zm`, `Tobacco`, `Sc`, `Cg`

### Output

Each mutation generates a file with primer pairs showing:
- **Primer name** - Mutation and orientation (5' or 3')
- **Sequence** - The oligonucleotide sequence
- **Length** - Number of bases
- **GC** - GC content percentage
- **Tm** - Melting temperature
- **Overlap_Tm** - Overlap region Tm for annealing

## Project Structure

```
sdm-mutagenesis/
├── src/                    # Source code modules
│   ├── main.py            # Entry point
│   ├── dna_utils.py       # DNA sequence handling
│   ├── codon_utils.py     # Codon preference tables
│   ├── mutation_utils.py  # Mutation parsing
│   ├── primer_utils.py    # Primer generation
│   └── output_utils.py    # File output
├── data/
│   ├── codon_preferences/ # Codon tables by species
│   └── examples/          # Example input files
├── output/                # Generated primer files
├── tests/                 # Unit tests
├── docs/                  # Documentation
└── legacy/                # Previous versions
```

## Development

### Running Tests

```bash
python3 -m pytest tests/
```

### Adding New Species

1. Update `data/codon_preferences/CodonPreference.txt`
2. Add species code to supported list
3. Test with example sequences

## Requirements

- Python 3.7+
- BioPython (for sequence analysis)

## License

See LICENSE file for details

## Contributing

Contributions welcome! Please ensure:
- Code follows PEP 8 style guide
- Functions have docstrings
- Tests cover new functionality

## References

- Overlap Extension PCR primer design
- Codon usage optimization for expression
