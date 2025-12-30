# NCBI Motif Finder API

A simple Python project that **consumes the NCBI Entrez (E-utilities) Web API**
to retrieve nucleotide sequences and search for DNA motifs in both the forward
strand and the reverse complement, reporting motif counts, positions, and densities.

---

## Features

- Fetch DNA sequences via the NCBI Entrez Web API using accession IDs
- Search for motifs in forward and reverse complement strands
- Report motif positions (zero-based), counts, sequence length, and motif density
- Graceful handling of network errors or missing sequences
- Optional export of results to a CSV file
  
---

## Repository Structure

- `motif_finder.py`  
  Contains the core functions:
  - `fetch_sequences_from_ncbi(accessions)` → fetch DNA sequences from NCBI Entrez Web API
  - `find_motif_positions(seq, motif)` → find all occurrences of a motif in a sequence
  - `reverse_complement(seq)` → compute the reverse complement of a sequence
  - `find_motif_in_sequence_list(motif, accessions, save=False)` → combines sequence retrieval and motif analysis, returning a pandas DataFrame

- `example_run.ipynb`  
  A Jupyter notebook demonstrating:
  - How to use `motif_finder.py`
  - Running motif searches on example accession numbers
  - Displaying results as a pandas DataFrame
  - Saving the results to CSV

---
