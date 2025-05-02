"""

CSC448 - Bioinformatics
Project 2. Wine making yeast Part 1.1
@Author: Oliver Pilon
Date: May 1, 2025

"""

# Library Imports
import numpy as np      #NumPy Documentation Link: https://numpy.org/doc/2.2/numpy-user.pdf       
import scipy as sp      #SciPy Documentation Link: https://docs.scipy.org/doc/scipy/tutorial/index.html#user-guide
from scipy.cluster.hierarchy import linkage, leaves_list,fcluster
from sklearn.metrics import silhouette_score
import pandas as pd
import matplotlib.pyplot as plt

data230 = pd.read_csv('Lab2/data/230genes_log_expression.txt',sep='\t', usecols=[1,2,3,4,5,6,7,8,9])
my230 = pd.read_csv('Lab2/top_230_genes_log.txt',sep=',')

data230['Name'] = data230['Name'] + data230['Extension'].fillna('').apply(lambda x: f" {x}" if x else '')

# my230 = pd.read_csv('Lab2/data/230genes_log_expression.txt',sep='\t', usecols=[1,2,3,4,5,6,7,8,9])
# my230['Name'] = my230['Name'] + my230['Extension'].fillna('').apply(lambda x: f" {x}" if x else '')
# my230 = my230.drop(columns=['Extension'])

# Find which merged names are duplicated
duplicates = data230['Name'][data230['Name'].duplicated(keep=False)]
# print(duplicates)

# Drop Extension
data230 = data230.drop(columns=['Extension'])

# print(f"data230 rows: {data230.shape[0]}")
# print(f"my230 rows: {my230.shape[0]}")

# print(data230['Name'].nunique())
# print(my230['Name'].nunique())


overlap = set(data230['Name']).intersection(set(my230['Name']))

print(f"Number of overlapping values: {len(overlap)}")



# Jaccard Index

jaccard = len(overlap) / (data230['Name'].nunique() + my230['Name'].nunique() - len(overlap))
print(f"Jaccard Coefficient: {jaccard}")





#Plots
numeric_data230 = data230.iloc[:, 1:].to_numpy()
numeric_my230 = my230.iloc[:,1:].to_numpy()

fig, axes = plt.subplots(2, 2, figsize=(16, 14))

axes[0,0].imshow(numeric_data230, aspect="auto", vmin=0, vmax=1, cmap="plasma")
axes[0,0].set_title("Provided Data Pre-Cluster")
axes[0,0].set_xticks(np.arange(numeric_data230.shape[1]))
axes[0,0].set_yticks([])

n_samples = numeric_data230.shape[0]
cluster_range = range(2, min(15, n_samples))

silhouette_scores = []
c_data = linkage(numeric_data230, method="average")
for k in cluster_range:
    labels = fcluster(c_data, k, criterion='maxclust')
    score = silhouette_score(numeric_data230, labels)
    silhouette_scores.append(score)

best_k = list(cluster_range)[np.argmax(silhouette_scores)]
print(f"Best number of clusters for provided 230 by silhouette score: {best_k}")
print(f"Silhouette score: {max(silhouette_scores):.4f}")

ordered = leaves_list(c_data)
clustered = numeric_data230[ordered, :]
axes[0,1].imshow(clustered, aspect="auto", vmin=0, vmax=1, cmap="plasma")
axes[0,1].set_title("Provided Data Clustered")
axes[0,1].set_xticks(np.arange(clustered.shape[1]))
axes[0,1].set_yticks([])

#My 230
axes[1,0].imshow(numeric_my230, aspect="auto", vmin=0, vmax=1, cmap="plasma")
axes[1,0].set_title("My Data Pre-Cluster")
axes[1,0].set_xticks(np.arange(numeric_my230.shape[1]))
axes[1,0].set_yticks([])

n_samples = numeric_my230.shape[0]
cluster_range = range(2, min(15, n_samples))

silhouette_scores = []
c_data = linkage(numeric_my230, method="average")
for k in cluster_range:
    labels = fcluster(c_data, k, criterion='maxclust')
    score = silhouette_score(numeric_my230, labels)
    silhouette_scores.append(score)

best_k = list(cluster_range)[np.argmax(silhouette_scores)]
print(f"Best number of clusters for my 230 by silhouette score: {best_k}")
print(f"Silhouette score: {max(silhouette_scores):.4f}")

ordered = leaves_list(c_data)
clustered = numeric_my230[ordered, :]
axes[1,1].imshow(clustered, aspect="auto", vmin=0, vmax=1, cmap="plasma")
axes[1,1].set_title("My Data Clustered")
axes[1,1].set_xticks(np.arange(clustered.shape[1]))
axes[1,1].set_yticks([])

# Layout and display
plt.tight_layout()
plt.savefig('Lab2/230Compare.png', bbox_inches="tight")
plt.show()

extent = axes[1,0].get_window_extent().transformed(fig.dpi_scale_trans.inverted())
fig.savefig('Lab2/my230_precluster.png', bbox_inches=extent.expanded(1.025,1.105))

extent = axes[1,1].get_window_extent().transformed(fig.dpi_scale_trans.inverted())
fig.savefig('Lab2/my230_cluster.png', bbox_inches=extent.expanded(1.025,1.105))