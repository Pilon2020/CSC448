"""

CSC448 - Bioinformatics
Project 2. Wine making yeast Part 1
@Author: Oliver Pilon
Date: April 22, 2025

"""

# Library Imports
import numpy as np      #NumPy Documentation Link: https://numpy.org/doc/2.2/numpy-user.pdf       
import scipy as sp      #SciPy Documentation Link: https://docs.scipy.org/doc/scipy/tutorial/index.html#user-guide
from scipy.cluster.hierarchy import linkage, leaves_list,fcluster
from sklearn.metrics import silhouette_score
import pandas as pd
import matplotlib.pyplot as plt
from enum import IntEnum

df = pd.read_csv('Lab2/data/diauxic_raw_ratios.txt',sep='\t',usecols=[1,3,4,5,6,7,8,9])


names = df['Name'].tolist()
ratios = df.columns[1:].tolist()        # ['R1.Ratio', …, 'R7.Ratio']
data   = df.iloc[:, 1:].astype(float).to_numpy()


fig, axes = plt.subplots(1, 2, figsize=(16, 14))

axes[0].imshow(data, aspect="auto", vmin=0, vmax=1, cmap="plasma")
axes[0].set_title("Unnormalized Pre-Cluster")
axes[0].set_xticks(np.arange(data.shape[1]))
axes[0].set_yticks([])

n_samples = data.shape[0]
cluster_range = range(2, min(15, n_samples))

silhouette_scores = []
c_data = linkage(data, method="average")
for k in cluster_range:
    labels = fcluster(c_data, k, criterion='maxclust')
    score = silhouette_score(data, labels)
    silhouette_scores.append(score)

best_k = list(cluster_range)[np.argmax(silhouette_scores)]
# print(f"Best number of clusters by silhouette score: {best_k}")
# print(f"Silhouette score: {max(silhouette_scores):.4f}")

ordered = leaves_list(c_data)
clustered = data[ordered, :]
axes[1].imshow(clustered, aspect="auto", vmin=0, vmax=1, cmap="plasma")
axes[1].set_title("Unnormalized Clustered")
axes[1].set_xticks(np.arange(clustered.shape[1]))
axes[1].set_yticks([])

# Layout and display
plt.tight_layout()
# plt.show()

## =============================================================================
## How to determine if a gene is of interest
## Base it on the % change from the first measurement to the final measurement. Since
## all the data is already normalied, this should identify the genes that have the
## largest change. Talking the absolute value of the genes will provide us with the
## largest 230 changes, link that using the label to the % change table, which will
## be used to derive the true list of largest changes. Create the list with the gene
## label, and the % its expression is increased or decreased. 
perchange = []

for i in range(data.shape[0]):
    change = ((data[i,-1]-data[i,0])/data[i,0])*100    
    perchange.append(change)

perchange_df = pd.DataFrame({
    'Gene': names,
    'Percent Change': perchange
})

perchange_df['AbsoluteChange'] = perchange_df['Percent Change'].abs()
perchange_df = perchange_df.sort_values(by='AbsoluteChange', ascending=False)
perchange_df.to_csv('Lab2/percent_change_full.txt', index=False, sep=',')

top230_df = perchange_df.head(230)
top230_df = top230_df.reset_index(drop=True)
top230_df.to_csv('Lab2/top_230_genes.txt', index=False, sep=',')

