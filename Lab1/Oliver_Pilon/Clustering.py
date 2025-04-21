import numpy as np
from matplotlib import pyplot as plt
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import dendrogram
from sklearn.metrics import silhouette_score

# 1. Load normalized SW scores
#    Assumes 'nSW_scores_2.csv' in working directory
norm_scores = np.loadtxt('nSW_scores_2.csv', delimiter=',')

# 2. Convert similarity to distance matrix
#    dist_scores[i,j] = 1 - normalized_score
dist_scores = 1.0 - norm_scores
linkage_type = 'single'
# 3. (Optional) Load sequence IDs for labeling dendrogram leaves
#    Expects a tab-delimited file 'protein_sequences.txt' with ID\tsequence
try:
    ids = []
    with open('protein_sequences.txt') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 1:
                ids.append(parts[0])
except FileNotFoundError:
    ids = None

# 4. Determine range of cluster counts to test
n_samples = dist_scores.shape[0]
cluster_range = range(2, min(15, n_samples))

# 5. Compute silhouette scores for each k
silhouette_scores = []
for k in cluster_range:
    model = AgglomerativeClustering(
        n_clusters=k,
        metric='precomputed',
        linkage=linkage_type,  # change linkage as desired
        compute_distances=True
    ).fit(dist_scores)
    score = silhouette_score(dist_scores, model.labels_, metric='precomputed')
    silhouette_scores.append(score)

# 6. Plot silhouette score vs. number of clusters
plt.figure(figsize=(15, 10))
plt.plot(list(cluster_range), silhouette_scores, marker='o')
plt.xlabel('Number of Clusters')
plt.ylabel('Silhouette Score')
plt.title('Silhouette Score vs. Number of Clusters')
plt.tight_layout()
plt.savefig(f'SilhouetteScores_{linkage_type}.png')
plt.show()

# 7. Identify best number of clusters
best_k = list(cluster_range)[np.argmax(silhouette_scores)]
print(f"Best number of clusters by silhouette score: {best_k}")

# 8. Perform final clustering with best_k
final_model = AgglomerativeClustering(
    n_clusters=best_k,
    metric='precomputed',
    linkage=linkage_type,
    compute_distances=True
).fit(dist_scores)

# 9. Build linkage matrix for dendrogram
counts = np.zeros(final_model.children_.shape[0])
n_samples = len(final_model.labels_)
for i, merge in enumerate(final_model.children_):
    current_count = 0
    for child_idx in merge:
        if child_idx < n_samples:
            current_count += 1
        else:
            current_count += counts[child_idx - n_samples]
    counts[i] = current_count

linkage_matrix = np.column_stack([
    final_model.children_,
    final_model.distances_,
    counts
]).astype(float)

# 10. Compute cutoff for exactly best_k clusters
distances = linkage_matrix[:, 2]
distances_sorted = np.sort(distances)
# the distance at which tree goes from k+1 to k clusters is distances_sorted[-(best_k-1)]
cutoff = distances_sorted[-(best_k - 1)]

# 11. Plot dendrogram with colored clusters
plt.figure(figsize=(15,10))
plt.title(f"Hierarchical Clustering Dendrogram (k={best_k}) Using {linkage_type}")

ddata = dendrogram(
    linkage_matrix,
    labels=ids,
    leaf_rotation=90,
    leaf_font_size=5,
    color_threshold=cutoff,
    p=3
)
plt.axhline(y=cutoff, linestyle='--', linewidth=0.5)
plt.tight_layout()
plt.savefig(f'Dendrogram_{linkage_type}.png', dpi=150)
plt.show()

import numpy as np
from scipy.cluster.hierarchy import to_tree

def get_splits(node, n_leaves):
    splits = set()
    def recurse(n):
        if n.is_leaf():
            return {n.id}
        left = recurse(n.get_left())
        right = recurse(n.get_right())
        split = frozenset(left if len(left) <= len(right) else right)
        if 0 < len(split) < n_leaves:
            splits.add(split)
        return left.union(right)
    recurse(node)
    return splits

def split_based_distance(Z1, Z2, n_leaves):
    tree1 = to_tree(Z1, rd=False)
    tree2 = to_tree(Z2, rd=False)
    splits1 = get_splits(tree1, n_leaves)
    splits2 = get_splits(tree2, n_leaves)
    return len(splits1.symmetric_difference(splits2))

def build_linkage_matrix(model):
    n = len(model.labels_)
    counts = np.zeros(model.children_.shape[0])
    for i, merge in enumerate(model.children_):
        counts[i] = sum((c < n and 1) or counts[c-n] for c in merge)
    return np.column_stack([model.children_, model.distances_, counts]).astype(float)

# --- Compare 'average' vs 'complete' and 'average' vs 'single' linkage trees ---

# Ensure best_k and dist_scores are defined from above

# Build models
avg_model = AgglomerativeClustering(
    n_clusters=best_k,
    metric='precomputed',
    linkage='average',
    compute_distances=True
).fit(dist_scores)

single_model = AgglomerativeClustering(
    n_clusters=best_k,
    metric='precomputed',
    linkage='single',
    compute_distances=True
).fit(dist_scores)

complete_model = AgglomerativeClustering(
    n_clusters=best_k,
    metric='precomputed',
    linkage='complete',
    compute_distances=True
).fit(dist_scores)

# Build linkage matrices
Z_avg     = build_linkage_matrix(avg_model)
Z_single  = build_linkage_matrix(single_model)
Z_complete= build_linkage_matrix(complete_model)

n_leaves = dist_scores.shape[0]

# Compute split-based distances
dist_avg_vs_complete = split_based_distance(Z_avg, Z_complete, n_leaves)
dist_avg_vs_single   = split_based_distance(Z_avg, Z_single, n_leaves)

print(f"Split-based distance between 'average' and 'complete': {dist_avg_vs_complete}")
print(f"Split-based distance between 'average' and 'single':   {dist_avg_vs_single}")
