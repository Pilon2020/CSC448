"""
CSC448 - Bioinformatics
Project 1 - Sequences and Evolutionary Trees
@Author: Oliver Pilon
"""

# Import libraries
import numpy as np      # NumPy Documentation Link: https://numpy.org/doc/2.2/numpy-user.pdf       
import scipy as sp      # SciPy Documentation Link: https://docs.scipy.org/doc/scipy/tutorial/index.html#user-guide
import sklearn as sk    # sklearn Documentation Link: https://scikit-learn.org/stable/user_guide.html
import blosum as bl
import requests

from matplotlib import pyplot as plt
from tqdm import tqdm
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import dendrogram


# Use BLOSUM62 from https://github.com/not-a-feature/blosum
matrix = bl.BLOSUM(62)

# Initializing constants/arrays
sequence_data = []
gap_penalty = -12  # From project description

url = 'https://raw.githubusercontent.com/CalPoly-MLBio/CSC448_Spring2025/refs/heads/main/data/protein_sequences.txt'
response = requests.get(url)

## LOADING DATA SECTION
# Load in data set
for line in response.text.strip().split('\n'):
    parts = line.strip().split('\t')
    if len(parts) == 2:
        sequence_data.append((parts[0], parts[1]))
            
sequence_array = np.array(sequence_data, dtype=object)
bacteria_ids = sequence_array[:, 0]

 
## FUNCTIONS

# Smith-Waterman Algorithm
def smith_waterman(seq1, seq2, blosum62, gap_penalty):
    row = len(seq1) + 1
    col = len(seq2) + 1

    score_matrix = np.zeros((row, col), dtype=int)
    tracing = np.zeros((row, col), dtype=int)

    max_score = 0
    
    for i in range(1, row):
        for j in range(1, col):
            matched = blosum62[seq1[i-1]][seq2[j-1]]
            score, direction = max(
                (score_matrix[i-1][j-1] + matched, 3),  # DIAG - match
                (score_matrix[i-1][j] + gap_penalty, 2),   # UP
                (score_matrix[i][j-1] + gap_penalty, 1),   # LEFT
                (0, 0)                                    # None
            )
            score_matrix[i, j] = score
            tracing[i, j] = direction

            if score > max_score:
                max_score = score
    return max_score

# Normalizing Smithwaterman function
def normalize(seq1, seq2):
    sw = smith_waterman(seq1, seq2, matrix, gap_penalty)
    sw_1 = smith_waterman(seq1, seq1, matrix, gap_penalty)
    sw_2 = smith_waterman(seq2, seq2, matrix, gap_penalty)
    max_sw = max(sw_1, sw_2)
    normalized = sw / max_sw
    return normalized

#dendrogram plot
def plot_dendrogram(model, labels, **kwargs):
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack(
        [model.children_, model.distances_, counts]
    ).astype(float)

    # Pass labels to the dendrogram function using the `labels` argument.
    dendrogram(linkage_matrix, labels=labels, **kwargs)

# Compute similarity scores with a progress bar that updates per sequence comparison
test_len = 10
num_sequences = len(sequence_array) # test_len
dist_scores = [[0.0 for _ in range(num_sequences)] for _ in range(num_sequences)]
norm_scores = [[0.0 for _ in range(num_sequences)] for _ in range(num_sequences)]
total_iterations = num_sequences + ((num_sequences*(num_sequences-1))/2)

with tqdm(total=total_iterations, desc="Computing distance scores") as pbar:
    for i in range(num_sequences):
        for j in range(i, num_sequences):
            seq1 = sequence_array[i][1]
            seq2 = sequence_array[j][1]
            normalized = normalize(seq1, seq2)
            norm_scores[i][j] = normalized  
            norm_scores[j][i] = normalized
            dist = 1 - float(normalized)
            dist_scores[i][j] = dist
            dist_scores[j][i] = dist
            pbar.update(1)
        np.savetxt("nSW_scores_2.csv", np.array(norm_scores), delimiter=",")


# Save the normalized Smith-Waterman scores as a CSV file

print("Normalized Smith-Waterman scores saved")

clustering = AgglomerativeClustering(metric='precomputed',linkage='complete',compute_distances=True).fit(dist_scores)
#Clustering Agglomerative Guide : https://scikit-learn.org/stable/modules/generated/sklearn.cluster.AgglomerativeClustering.html

plt.title(f"Hierarchical Clustering Dendrogram Using {clustering.linkage}")
plot_dendrogram(clustering)
plt.tight_layout()
plt.show()