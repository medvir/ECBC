# ECBC
Collection of scripts for error correcting antibody sequening

## IgG
- R1: fwd 4N, leader, variable
- R2: rev subtype determination, constant, variable
- I1: EVBC 12
- I2: sample

## klMA
- R1: fwd 4N, leader, variable
- R2: rev ECBC 12, spacer 4, constant, variable
- I1: LC index TAAGGCGAGAGC 12
- I2: sample

## workflow

1. pandaseq R1 und R2
2. split into IgG vs. klMA
3a HC IgG subtypes
3b. LC get ECBE, klMA discrimination
4. count indices, only indices present >10 ?
5. collapse indices —> consensus (check for variation)
