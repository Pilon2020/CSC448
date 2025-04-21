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
from sklearn.cluster import AgglomerativeClustering # Clustering Agglomerative Guide : https://scikit-learn.org/stable/modules/generated/sklearn.cluster.AgglomerativeClustering.html
from scipy.cluster.hierarchy import dendrogram
from sklearn.metrics import silhouette_score

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
def normalize(seq1, seq2, matrix, gap_penalty):
    sw = smith_waterman(seq1, seq2, matrix, gap_penalty)
    sw_1 = smith_waterman(seq1, seq1, matrix, gap_penalty)
    sw_2 = smith_waterman(seq2, seq2, matrix, gap_penalty)
    max_sw = max(sw_1, sw_2)
    normalized = sw / max_sw
    return normalized

#dendrogram plot
def plot_dendrogram(model, **kwargs):
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
    dendrogram(linkage_matrix, **kwargs)

# Compute similarity scores with a progress bar that updates per sequence comparison
test_len = 20
num_sequences = len(sequence_array) # test_len #  
dist_scores = np.zeros((num_sequences, num_sequences))
norm_scores = np.zeros((num_sequences, num_sequences))
np.fill_diagonal(norm_scores,1)
total_iterations = ((num_sequences*(num_sequences-1))/2)

with tqdm(total=total_iterations, desc="Computing distance scores") as pbar:
    for i in range(num_sequences):
        for j in range(i+1, num_sequences):
            seq1 = sequence_array[i][1]
            seq2 = sequence_array[j][1]
            normalized = normalize(seq1, seq2, matrix, gap_penalty)
            norm_scores[i][j] = normalized  
            norm_scores[j][i] = normalized
            dist = 1 - float(normalized)
            dist_scores[i][j] = dist
            dist_scores[j][i] = dist
            pbar.update(1)
        np.savetxt("nSW_scores_2.csv", np.array(norm_scores), delimiter=",")



# Silhouette score vs. number of clusters
linkage_type = 'complete'  # options: 'ward', 'complete', 'average', 'single'
cluster_range = range(2, min(15, num_sequences))

# Silhouette analysis
silhouette_scores = []
for k in cluster_range:
    model = AgglomerativeClustering(
        n_clusters=k,
        metric='precomputed',
        linkage=linkage_type,
        compute_distances=True
    ).fit(dist_scores)
    score = silhouette_score(dist_scores, model.labels_, metric='precomputed')
    silhouette_scores.append(score)

# Plot silhouette scores
plt.figure("Silhouette Scores")
plt.plot(list(cluster_range), silhouette_scores, marker='o')
plt.xlabel('Number of Clusters')
plt.ylabel('Silhouette Score')
plt.title('Silhouette Score vs Number of Clusters')
plt.savefig(f'SilhouetteScores_{linkage_type}.png',dpi='figure',format = 'png',bbox_inches='tight',pad_inches=0.2)
plt.tight_layout()
plt.show()

# Final clustering and dendrogram
best_k = cluster_range[np.argmax(silhouette_scores)]
final_model = AgglomerativeClustering(
    n_clusters=best_k,
    metric='precomputed',
    linkage=linkage_type,
    compute_distances=True
).fit(dist_scores)

# Dendrogram plot
def plot_dendrogram(model, **kwargs):
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        counts[i] = sum((c < n_samples and 1) or counts[c-n_samples] for c in merge)
    linkage_matrix = np.column_stack([model.children_, model.distances_, counts]).astype(float)
    dendrogram(linkage_matrix, **kwargs)

plt.figure("Dendrogram",figsize=(12,8), dpi=300)
plt.title(f"Hierarchical Clustering Dendrogram (k={best_k}, linkage={linkage_type})")
plot_dendrogram(final_model, labels=bacteria_ids[:num_sequences],leaf_rotation=90,leaf_font_size=8)
plt.savefig(f'Dendrogram_{linkage_type}.png',dpi='figure',format = 'png',bbox_inches='tight',pad_inches=0.2)
plt.tight_layout()
plt.show()
