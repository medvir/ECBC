# ECBC
Script subtype assignment and primer (constant region) trimming for antibody sequences wiht Error Correcting BarCodes

## Experimental Setup
### ECBC IgG Reads
- R1: fwd 4N, leader, variable
- R2: rev subtype determination, constant, variable
- I1: ECBC (12)
- I2: sample index

### ECBC klMA Reads
- R1: fwd 4N, leader, variable
- R2: rev ECBC (12), spacer (4), constant, variable
- I1: LC index TAAGGCGAGAGC (12)
- I2: sample index

## Workflow

`IgG_klMA_ECBC.py`
- split into IgG vs. klMA
- for LC do klMA discrimination
- for HC determine IgG subtypes
- write ECBC in front of R1
- pandaseq R1 und R2

## Usage
`IgG_klMA_ECBC_V2.py NAME_L001_R1_001.FASTQ[.GZ]`

###Helper scritps
`seq_count.sh` count fasta and fastq sequences

### Primer
- `primer_fwd_H.fasta` heavy chain forward primers
- `primer_fwd_k.fasta` kappa forward primers
- `primer_fwd_l.fasta` lambda forward primers