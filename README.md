# DNA_Protein_Alignment

The following libraries are used: Biopython, NumPy, and sys. Sys should be included if Python is installed on your device.

To install NumPy, run the following command in your terminal:
```
pip install numpy
```

To install Biopython, run the following command in your terminal:
```
pip install biopython
```

To run the code, download the alignment.py file to your device. In your terminal, navigate to the folder containing the alignment.py file and run the following command:
```
 python alignment.py <protein_sequences>.fasta <dna_sequences>.fasta <l/g> [gap_open_penalty] [gap_extend_penalty]
```
- <protein_sequences>.fasta contains the sequences you want to compare in the following format: <br>
  \>Organism 1 <br>
  MNIT...<br>
  \>Organism 2 <br>
  MNIT...
- <dna_sequences>.fasta contains the sequences you want to compare in the same format as the protein sequences
- 
  
