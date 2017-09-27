# plow
A phylodynamics workflow for pathogen sequence data

Phlow can be called using the following command from a unix terminal:

INFILE='fileName' snakemake -s phlow.snakefile

The input sequence file should have the name "fileName.fas" and the sequences should be in FASTA format. 

The header line for the sequences should be in the following format:

\>seqName_sampleDate
