"""

CSC448 - Bioinformatics
Project 1 - Sequences and Evolutionary Trees
@Author: Oliver Pilon

"""
#import libraries
import numpy as np      #NumPy Documentation Link: https://numpy.org/doc/2.2/numpy-user.pdf       
import scipy as sp      #SciPy Documentation Link: https://docs.scipy.org/doc/scipy/tutorial/index.html#user-guide
import sklearn as sk    #sklearn Documentation Link: https://scikit-learn.org/stable/user_guide.html
import blosum as bl
import requests

# Use BLOSUM62 from https://github.com/not-a-feature/blosum
matrix = bl.BLOSUM(62)

#Initiallizing Constants/Arrays
sequence_data = []
gap_penalty = -12 #From Project description

url = 'https://raw.githubusercontent.com/CalPoly-MLBio/CSC448_Spring2025/refs/heads/main/data/protein_sequences.txt'

response = requests.get(url)


## LOADING DATA SECTION
#Load in data set
for line in response.text.strip().split('\n'):
    parts = line.strip().split('\t')
    if len(parts) == 2:
        sequence_data.append((parts[0], parts[1]))
            
sequence_array = np.array(sequence_data, dtype=object)

## FUNCTIONS

#Smith Waterman Algorithm
def smith_waterman(seq1, seq2, blosum62, gap_penalty):
    row = len(seq1)+1
    col = len(seq2)+1

    matrix = np.zeros(shape=(row, col), dtype=int)
    tracing = np.zeros(shape=(row, col), dtype=int)

    max_score = 0
    max_index = (0, 0)
    
    for i in range(1,row):
        for j in range(1,col):
            matched = blosum62[seq1[i-1]][seq2[j-1]]

            matrix[i,j], tracing[i,j] = max (
                (matrix[i-1][j-1] + matched, 3),    # DIAG - Match
                (matrix[i-1][j] + gap_penalty, 2),  # UP
                (matrix[i][j-1] + gap_penalty, 1),  # LEFT
                (0, 0)                              # None
            )
            
            if matrix[i, j] >= max_score:
                max_index = (i,j)
                max_score = matrix[i, j]
            sw_score = matrix[i,j]    
    return sw_score

# Normalizing Function
def normalize(seq1, seq2):
    sw = smith_waterman(seq1, seq2, matrix, gap_penalty)
    sw_1 = smith_waterman(seq1, seq1, matrix, gap_penalty)
    sw_2 = smith_waterman(seq2, seq2, matrix, gap_penalty)
    max_sw = max(sw_1,sw_2)
    normalized = sw/max_sw
    return normalized

sim_scores = [[0 for _ in range(len(sequence_array))] for _ in range(len(sequence_array))]

for i in range(len(sequence_array)):
    for j in range(len(sequence_array)):
        seq1 = sequence_array[i][1] #String 1 of DNA for comparison  
        seq2 = sequence_array[j][1] #String 2 of DNA for comparison 
        normalized = normalize(seq1, seq2)
        sim_scores[i][j] = float(normalized)

# ## Small Scale Testing
# size = 4
# sim_scores = [[0 for _ in range(size)] for _ in range(size)]

# for i in range(size):
#     for j in range(size):
#         seq1 = sequence_array[i][1] #String 1 of DNA for comparison  
#         seq2 = sequence_array[j][1] #String 2 of DNA for comparison 
#         normalized = normalize(seq1, seq2)
#         sim_scores[i][j] = float(normalized)
        
print(sim_scores)