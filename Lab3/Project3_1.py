import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from itertools import combinations

# 1) Load the PPI graph
G = nx.Graph()
with open('interacting_proteins.txt') as f:
    for line in f:
        u, v = line.strip().split()
        G.add_edge(u, v)

# 2) Load known drivers
with open('onco_genes.txt') as f:
    oncogenes = [g.strip() for g in f if g.strip()]
oncogenes = [g for g in oncogenes if g in G]
if len(oncogenes) < 2:
    raise RuntimeError("Need at least two oncogenes present in the graph")

# 3) Helper to compute avg shortest-path over all pairs
def avg_shortest_path(nodes):
    dists = []
    for a, b in combinations(nodes, 2):
        try:
            dists.append(nx.shortest_path_length(G, a, b))
        except nx.NetworkXNoPath:
            pass
    return np.mean(dists)

# 4) Observed score
obs_score = avg_shortest_path(oncogenes)
print(f"Observed avg shortest-path: {obs_score:.3f}")

# 5) Null distribution (using numpy.random)
all_nodes  = np.array(list(G.nodes()))
mask       = ~np.isin(all_nodes, oncogenes)
candidates = all_nodes[mask]

# fix seed for reproducibility
rng = np.random.default_rng(42)

null_scores = []
for _ in range(1000):
    # sample without replacement
    sample = rng.choice(candidates, size=len(oncogenes), replace=False)
    null_scores.append(avg_shortest_path(sample))
null_scores = np.array(null_scores)

# 6) Compute p-value
p_value = np.mean(null_scores <= obs_score)
print(f"Null mean: {null_scores.mean():.3f}")
print(f"P-value (fraction ≤ observed): {p_value:.4f}")

# 7) Plot the null distribution
plt.figure(figsize=(8,5))
plt.hist(null_scores, bins=30, alpha=0.75, edgecolor='k')
plt.axvline(obs_score, color='red', linewidth=2,
            label=f'Observed = {obs_score:.2f}')
plt.xlabel('Average shortest-path')
plt.ylabel('Count (out of 1000)')
plt.title('Null distribution of avg shortest-path\n(random gene-sets via NumPy)')
plt.legend()
plt.tight_layout()
plt.show()

