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
- <dna_sequences>.fasta contains the sequences you want to compare in the same format as the protein sequences.
- <l/g> specifies which type of alignment: local (l) or global (g). If you enter an incorrect character, the program defaults to a global alignment.
- The last two arguments specify the gap open and gap extend penalty. <b>Both of these arguments are optional</b>. If you choose not to enter anything, the program defaults to a -2 penalty for both. You can either enter both or none of these parameters--the program does not accept one of these parameters.

If you enter the parameters incorrectly, the usage is printed out.

Other error checking of the parameters is not there. For errors that have to do with incorrect file inputs, the program will either give an incorrect output or crash.

For the file inputs: The protein fasta file should have the same number of sequences as the DNA fasta file. The first sequence in the protein file will be compared with the first sequence in the DNA file, the second sequences in both files will be compared, and so on.

The output of the program will be the alignment. The rows starting with 'g' mark the given protein sequence (which came from the protein fasta file) and the 't' rows mark the translated protein sequence (which came from the DNA fasta file).
