# ECBC

## IgG
R1 fwd 4N, leader, variable
R2 rev subtype determination, constant, variable
I1 EVBC 12
I2 sample

## klMA
R1 fwd 4N, leader, variable
R2 rev ECBC 12, spacer 4, constant, variable
I1 LC index TAAGGCGAGAGC 12
I2 sample

## workflow

1. pandaseq R1 und R2


2. split into IgG vs. klMA

HC
- IgG subtypes

LC
- get ECBE
- klMA discrimination

- count indices, only indices present >10 ?
- collapse indices —> consensus (check for variation)
