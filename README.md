# ECBC
Collection of scripts for error correcting antibody sequening

## IgG Reads
- R1: fwd 4N, leader, variable
- R2: rev subtype determination, constant, variable
- I1: EVBC 12
- I2: sample

## klMA Reads
- R1: fwd 4N, leader, variable
- R2: rev ECBC 12, spacer 4, constant, variable
- I1: LC index TAAGGCGAGAGC 12
- I2: sample

## Workflow

`pipeline.ch`
- script that calls all the later ones

`pandaseq.sh`
combine R1 and R2 to full variable region sequence
Usage: `pandaseq.sh R1.FASTQ(.GZ)/DIRECTORY [MINIMAL_OVERLAP (default = 1)]`

`IgG_klMA_ECBC.py`
- pandaseq R1 und R2split into IgG vs. klMA
- for HC get IgG subtypes
- for LC get ECBE, klMA discrimination

`indices`
- count indices (only indices present >10 ?)
- collapse indices —> consensus (check for variation)
