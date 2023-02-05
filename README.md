# Basefinder

## Installation

## Usage
```
USAGE:
    basefinder [OPTIONS] --alignment <ALIGNMENT> --tsv <TSV> --ingroup <INGROUP> --outgroup <OUTGROUP>

OPTIONS:
    -a, --alignment <ALIGNMENT>
            Sequence alignment in FASTA-format containing one sequence per species
    -h, --help
            Print help information
    -i, --ingroup <INGROUP>
            Ingroup
        --no-header
            Hide header
    -o, --outgroup <OUTGROUP>
            Outgroup
    -s, --show-utr
            hide non-cds rows (only active in verbose mode)
    -t, --tsv <TSV>
            A _headerless_ tab-separated file containing 8 columns extracted from a VCF-file:

            #CHROM    POS    ALT    REF    HGVS.c    AF    VARIANT    STRAND

            Example:

            elon_1_1    166    G    A    c.-131G>A    0.0625    synonymous_variant    +
    -V, --version
            Print version information
```
## Example input & output
