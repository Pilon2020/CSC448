"""

CSC448 - Bioinformatics
Smith-Waterman Alignment Algorithm
@Author: Oliver Pilon

"""

# Library Imports
import numpy as np      #NumPy Documentation Link: https://numpy.org/doc/2.2/numpy-user.pdf       
import scipy as sp      #SciPy Documentation Link: https://docs.scipy.org/doc/scipy/tutorial/index.html#user-guide
import sklearn as sk    #sklearn Documentation Link: https://scikit-learn.org/stable/user_guide.html
from enum import IntEnum

dna1 = "TACGGGCCCGCTAC" #String 1 of DNA for comparison  
dna2 = "TAGCCCTATCGGTCA" #String 2 of DNA for comparison 

#Convert to lists of individual characters
seq1 = list(dna1) 
seq2 = list(dna2)

row = len(seq1)+1
col = len(seq2)+1

# Scoring Criteria
class Score(IntEnum):
    Match = 2
    Mismatch = -2
    Gap = -1
    GapStart = -5

class Trace(IntEnum):
    Stop = 0
    Left = 1
    Up = 2
    Diag = 3

matrix = np.zeros(shape=(row, col), dtype=int)
tracing = np.zeros(shape=(row, col), dtype=int)
max_score = 0
max_index = (0,0)


for i in range(1,row):
    for j in range(1,col):
        matched = Score.Match if seq1[i-1] == seq2[j-1] else Score.Mismatch

        if tracing[i-1,j] == 2 or tracing[i,j-1] == 1:
            gaplen = gaplen + 1
        else:
            gaplen = 0
        
        print(gaplen)    
        GapScore = (Score.Gap * gaplen) + Score.GapStart
        print(f'[{i},{j}]: {GapScore}')
        matrix[i,j], tracing[i,j] = max (
            (matrix[i-1][j-1] + matched, 3),
            (matrix[i-1][j] + Score.Gap, 2),
            (matrix[i][j-1] + Score.Gap, 1),
            (0, 0)
        )
        
        if matrix[i, j] >= max_score:
            max_index = (i,j)
            max_score = matrix[i, j]
        
#Taken from https://github.com/slavianap/Smith-Waterman-Algorithm/blob/master/Script.py
#Modified to work with variables
# Initialising the variables for tracing
aligned_seq1 = ""
aligned_seq2 = ""   
(max_i, max_j) = max_index

while tracing[max_i, max_j] != Trace.Stop:
    if tracing[max_i, max_j] == Trace.Diag:
        aligned_seq1 = seq1[max_i - 1] + aligned_seq1
        aligned_seq2 = seq2[max_j - 1] + aligned_seq2
        max_i -= 1
        max_j -= 1
    elif tracing[max_i, max_j] == Trace.Up:
        aligned_seq1 = seq1[max_i - 1] + aligned_seq1
        aligned_seq2 = '-' + aligned_seq2
        max_i -= 1
    elif tracing[max_i, max_j] == Trace.Left:
        aligned_seq1 = '-' + aligned_seq1
        aligned_seq2 = seq2[max_j - 1] + aligned_seq2
        max_j -= 1

print(matrix)
print("Tracing Array")
print(tracing)
print("Aligned Sequence 1:", aligned_seq1)
print("Aligned Sequence 2:", aligned_seq2)
