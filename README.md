# SDM Mutagenesis - Site-Directed Mutagenesis Primer Generator

A professional Python tool for designing mutagenic primers for site-directed mutagenesis (SDM) using overlap extension PCR. Based on work of [mds-mutagenesis](https://github.com/shiqiang-lin/sdm-mutagenesis) 

## Features

- **Automatic primer design** - Generates tetrad primer combinations for SDM
- **Species-specific codons** - Uses codon preference tables for optimal expression
- **Thermal analysis** - Calculates Tm values and GC content for each primer
- **Sorted output** - Primers ranked by overlap Tm for easy selection

## Quick Start

### Installation

```
# Clone the repository
git clone <repo>
```

### Usage


Basic usage:
```
python main.py data/examples/input_gene.txt data/examples/input_mutaa.txt Human
```

With custom output directory:
```
python main.py data/examples/input_gene.txt data/examples/input_mutaa.txt Human output/my_results
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



## Requirements

- Python 3.7+
- BioPython (for sequence analysis)


