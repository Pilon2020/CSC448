import numpy as np
import pandas as pd

# 1. Load normalized Smith–Waterman similarity scores
norm_scores = np.loadtxt('nSW_scores_2.csv', delimiter=',')

# 2. Load protein IDs
ids = []
with open('protein_sequences.txt') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 1:
            ids.append(parts[0])
ids = np.array(ids)

# 3. Extract upper-triangle indices (i < j)
triu_idx = np.triu_indices_from(norm_scores, k=1)
scores = norm_scores[triu_idx]

# 4. Identify two most similar (highest) and two most dissimilar (lowest) pairs
sorted_idx = np.argsort(scores)
lowest_idx = sorted_idx[:2]
highest_idx = sorted_idx[-2:]

# 5. Compile results
results = []
for category, idx_list in [('Most Dissimilar', lowest_idx), ('Most Similar', highest_idx)]:
    for idx in idx_list:
        i, j = triu_idx[0][idx], triu_idx[1][idx]
        results.append({
            "Category": category,
            "Protein 1": ids[i],
            "Protein 2": ids[j],
            "Similarity": f"{scores[idx]:.18f}"
        })

# 6. Display results
df = pd.DataFrame(results)
print(df.to_string(index=False))
