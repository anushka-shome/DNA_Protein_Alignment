from Bio import SeqIO
from Bio.Seq import Seq
import sys
import numpy as np

def read_file(file):
  seqs = {}
  for i in SeqIO.parse(file, "fasta"):
    seqs[i.id] = str(i.seq)
  return seqs

def gen_matrix(file):
  matrix = {}
  with open(file, 'r') as f:
    lines = f.readlines()

  aa = lines[6].strip().split()
  aa.remove("*")
  #print(aa)
  for line in lines[7:]:
    row = line.strip().split()
    row_aa = row[0]
    scores = row[1:]
    for col_aa, score in zip(aa, scores):
      matrix[(row_aa, col_aa)] = int(score)
  return matrix

def translate(dna):
  seq = Seq(dna)
  protein = seq.translate()
  return str(protein).replace("*", "W")

def local_alignment(seq1, seq2, matrix, go, gp):
  scores = np.zeros((len(seq1) + 1, len(seq2) + 1))
  gaps1 = np.zeros((len(seq1) + 1, len(seq2) + 1)) #seq1 gaps
  gaps2 = np.zeros((len(seq1) + 1, len(seq2) + 1)) #seq2 gaps

  '''
  for i in range(1, len(seq1)+1):
    scores[i][0] = go + (i-1) * gp
    gaps1[i][0] = scores[i][0]
  for i in range(1, len(seq2)+1):
    scores[0][i] = go + (i-1) * gp
    gaps2[0][i] = scores[0][i]
  '''

  max_score = 0
  max_i = 0
  max_j = 0

  for i in range(1, len(seq1)+1):
    for j in range(1, len(seq2)+1):
      match = scores[i-1][j-1] + matrix[(seq1[i-1], seq2[j-1])]
      delete = max(scores[i-1][j] + go, gaps1[i-1][j] + gp)
      gaps1[i][j] = delete
      insert = max(scores[i][j-1] + go, gaps2[i][j-1] + gp)
      gaps2[i][j] = insert
      scores[i][j] = max(match, delete, insert, 0)
      if scores[i][j] > max_score:
        max_score = scores[i][j]
        max_i = i
        max_j = j

  aseq1 = []
  aseq2 = []
  i = max_i
  j = max_j
  if i != len(seq1): #more letters in seq1 than seq2
    temp_i = len(seq1)
    while (temp_i != i):
      aseq1.append(seq1[temp_i - 1])
      aseq2.append("-")
      temp_i = temp_i - 1
  if j != len(seq2): #more letters in seq2 than seq1
    temp_j = len(seq2)
    while (temp_j != j):
      aseq1.append("-")
      aseq2.append(seq2[temp_j - 1])
      temp_j = temp_j - 1
 
  while i > 0 or j > 0:
    if i > 0 and (j == 0 or scores[i][j] == gaps1[i][j]): #gap in sequence 2
      aseq1.append(seq1[i-1])
      aseq2.append("-")
      i = i-1
    elif j > 0 and (i == 0 or scores[i][j] == gaps2[i][j]): #gap in sequence 1
      aseq1.append("-")
      aseq2.append(seq2[j-1])
      j = j-1
    elif i > 0 and j > 0 and scores[i][j] == scores[i-1][j-1] + matrix[(seq1[i-1], seq2[j-1])]:
      aseq1.append(seq1[i-1])
      aseq2.append(seq2[j-1])
      i = i-1
      j = j-1
    else:
      break
  '''
  while i > 0:
    aseq1.append(seq1[i - 1])
    aseq2.append("-")
    i = i - 1
  while j > 0:
    aseq1.append("-")
    aseq2.apped(seq2[j-1
  '''

  aseq1.reverse()
  aseq2.reverse()
  return aseq1, aseq2

def global_alignment(seq1, seq2, matrix, go, gp):
  scores = np.zeros((len(seq1) + 1, len(seq2) + 1))
  gaps1 = np.zeros((len(seq1) + 1, len(seq2) + 1))
  gaps2 = np.zeros((len(seq1) + 1, len(seq2) + 1))
  for i in range(1, len(seq1)+1):
    scores[i][0] = go + (i-1) * gp
    gaps1[i][0] = scores[i][0]
  for i in range(1, len(seq2)+1):
    scores[0][i] = go + (i-1) * gp
    gaps2[0][i] = scores[0][i]
  for i in range(1, len(seq1)+1):
    for j in range(1, len(seq2)+1):
      match = scores[i-1][j-1] + matrix[(seq1[i-1], seq2[j-1])]
      delete = max(scores[i-1][j] + go, gaps1[i-1][j] + gp)
      gaps1[i][j] = delete
      insert = max(scores[i][j-1] + go, gaps2[i][j-1] + gp)
      gaps2[i][j] = insert
      scores[i][j] = max(match, delete, insert)

  aseq1 = []
  aseq2 = []
  i = len(seq1)
  j = len(seq2)

  while i > 0 or j > 0:
    if i > 0 and (j == 0 or scores[i][j] == gaps1[i][j]): #gap in sequence 2
      aseq1.append(seq1[i-1])
      aseq2.append("-")
      i = i-1
    elif j > 0 and (i == 0 or scores[i][j] == gaps2[i][j]): #gap in sequence 1
      aseq1.append("-")
      aseq2.append(seq2[j-1])
      j = j-1
    else: #match/mismatch case
      aseq1.append(seq1[i-1])
      aseq2.append(seq2[j-1])
      i = i-1
      j = j-1

  aseq1.reverse()
  aseq2.reverse()
  return aseq1, aseq2


def output(dna, protein, protein_dna, aseq1, aseq2, p_name, d_name):
  print(f"{p_name} (protein) vs. {d_name} (DNA):")
  print("")
  print(f"DNA Sequence: {dna}")
  print("")
  print(f"Protein Sequence (given): {protein}")
  print("")
  print(f"Translated protein sequence (dna): {protein_dna}")
  print("")
  print(f"Alignment [{p_name} (g) vs. {d_name} (t)]:")
  for i in range(0, len(aseq1), 75):
    if aseq1[i] == " ":
      i = i+1 #skip the space if it is the beginning of the line
    top = aseq1[i:i+75]
    bottom = aseq2[i:i+75]
    middle = ""
    for j in range(0, len(top)):
      if top[j] == " ": #space
        middle += " "
      elif top[j] == bottom[j]: #match
        middle += "|"
      else:
        middle+= " "
    print(f"g:{aseq1[i:i+75]}")
    print(f"  {middle}")
    print(f"t:{aseq2[i:i+75]}")

  print("") #newline
  print("")

def main():
  go = -2
  gp = -2
  align = 1
  if len(sys.argv) != 4 and len(sys.argv) != 6:
    print("Usage: python alignment.py <protein_sequences>.fasta <dna_sequences>.fasta <l/g> [gap_open_penalty] [gap_extend_penalty]")
    return
  elif len(sys.argv) == 5:
    go = int(sys.argv[3])
    gp = int(sys.argv[4])
  if sys.argv == "l":
    align = 0 #defaults to global alignment
  protein = read_file(sys.argv[1])
  p_keys = list(protein.keys())
  #print(protein)
  dna = read_file(sys.argv[2])
  d_keys = list(dna.keys())
  #print(dna)
  matrix = gen_matrix("BLOSUM62")
  #print(matrix)
  #protein_dna = {}
  for i in range(0, min(len(p_keys), len(d_keys))):
    protein_dna = translate(dna[d_keys[i]])
    #print(protein_dna)
    #print(entry)
    if align == 1:
      [aseq1, aseq2] = global_alignment(protein[p_keys[i]], protein_dna, matrix, go, gp)
      aseq1 = " ".join(aseq1)
      aseq2 = " ".join(aseq2)
      output(dna[d_keys[i]], protein[p_keys[i]], protein_dna, aseq1, aseq2, p_keys[i], d_keys[i])
    else:
      [aseq1, aseq2] = local_alignment(protein[p_keys[i]], protein_dna, matrix, go, gp)
      aseq1 = " ".join(aseq1)
      aseq2 = " ".join(aseq2)
      output(dna[d_keys[i]], protein[p_keys[i]], protein_dna, aseq1, aseq2, p_keys[i], d_keys[i])

    #print(aseq1)
    #print(aseq2)
    #protein_dna[entry].replace("*", "W")
  #print(protein_dna)
  #print(protein)

if __name__ == "__main__":
  main()
