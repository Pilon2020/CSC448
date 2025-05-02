"""

CSC448 - Bioinformatics
Project 2. Wine making yeast Part 1
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

df = pd.read_csv('Lab2/data/diauxic_raw_ratios.txt',sep='\t',usecols=[1,2,3,4,5,6,7,8,9])
df['Name'] = df['Name'] + df['Extension'].fillna('').apply(lambda x: f" {x}" if x else '')
df = df.drop(columns=['Extension'])

names = df['Name'].tolist()
ratios = df.columns[1:].tolist()        # ['R1.Ratio', …, 'R7.Ratio']
data   = df.iloc[:, 1:].astype(float).to_numpy()


fig, axes = plt.subplots(1, 2, figsize=(16, 10))

axes[0].imshow(data, aspect="auto", vmin=0, vmax=1, cmap="plasma")
axes[0].set_title("Unnormalized Pre-Cluster")
axes[0].set_xticks(np.arange(data.shape[1]))
axes[0].set_yticks([])

n_samples = data.shape[0]
cluster_range = range(2, 15)

silhouette_scores = []
c_data = linkage(data, method="average")
for k in cluster_range:
    labels = fcluster(c_data, k, criterion='maxclust')
    score = silhouette_score(data, labels)
    silhouette_scores.append(score)

best_k = list(cluster_range)[np.argmax(silhouette_scores)]
print(f"Best number of clusters by silhouette score: {best_k}")
print(f"Silhouette score: {max(silhouette_scores):.10f}")

ordered = leaves_list(c_data)
clustered = data[ordered, :]
axes[1].imshow(clustered, aspect="auto", vmin=0, vmax=1, cmap="plasma")
axes[1].set_title("Unnormalized Clustered")
axes[1].set_xticks(np.arange(clustered.shape[1]))
axes[1].set_yticks([])

# Layout and display
plt.tight_layout()

ratio = np.log2(data[:, 1:] / data[:, :-1])

average_column = ratio.mean(axis=1).reshape(-1, 1)
average_with_avg = np.hstack((ratio, average_column))

ratio_df = pd.DataFrame({
    'Name': names,
    '1 - 2': average_with_avg[:, 0],
    '2 - 3': average_with_avg[:, 1],
    '3 - 4': average_with_avg[:, 2],
    '4 - 5': average_with_avg[:, 3],
    '5 - 6': average_with_avg[:, 4],
    '6 - 7': average_with_avg[:, 5],
    'Average': average_with_avg[:, 6],
})

ratio_df['Absolute'] = ratio_df['Average'].abs()
ratio_df = ratio_df.sort_values(by='Absolute', ascending=False)
ratio_df.to_csv('Lab2/fullgenesranked.txt', index=False, sep=',')

top230_names = ratio_df['Name'].head(230).tolist()
name_to_index = {name: idx for idx, name in enumerate(names)}
top230_data_rows = np.array([np.log2(data[name_to_index[name]]) for name in top230_names])

top230_df = pd.DataFrame(
    top230_data_rows,
    columns=['R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7','Absolute']
)
top230_df.insert(0, 'Name', top230_names)

top230_df.to_csv('Lab2/top_230_genes_log.txt', index=False, sep=',')

names = top230_df["Name"]
# Split off the first space
split_df = top230_df["Name"].str.split(" ", n=1, expand=True)
split_df.columns = ["name", "extension"]

# Drop rows where extension is missing, empty, or literally '""'
mask = (
    split_df["extension"].notna()                             # not NaN
    & (split_df["extension"].str.strip() != "")               # not blank
    & (split_df["extension"] != '""')                         # not the literal two-quote string
)
filtered_extensions = split_df.loc[mask, "extension"]

# Save to file
filtered_extensions.to_csv('Lab2/230names.txt', index=False, sep=',')

plt.savefig('Lab2/full_presort.png', bbox_inches='tight')
plt.show()