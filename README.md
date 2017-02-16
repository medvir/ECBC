# ECBC
Subtype assignment and primer trimming for antibody sequences with Error Correcting BarCodes (ECBC)

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
1. split into IgG vs. klMA
2. for LC do klMA discrimination
3. for HC determine IgG subtypes
4. write ECBC in front of R1
5. pandaseq R1 und R2

## Usage
`IgG_klMA_ECBC_V2.py NAME_L001_R1_001.FASTQ[.GZ]`

###Helper Scripts
`seq_count.sh` counts fasta and fastq sequences

`batch.ch` to run multiple samples

`imgt_combine.sh` to combine two IMGT output files from the same sample
 
### Primer
- `primer_fwd_H.fasta` heavy chain forward primers
- `primer_fwd_k.fasta` kappa forward primers
- `primer_fwd_l.fasta` lambda forward primers