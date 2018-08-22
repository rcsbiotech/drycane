# Perl helper stuff for drycane processing

Basic Perl utilities for drycane RNASeq processing (differential expression and annotation mostly)

- `getGI2.pl`: convert old NCBI gene identifiers (gi) to new accession numbers (AN). Made to be slow (~1 million GIs/hour), given NCBI's API limit policy.
