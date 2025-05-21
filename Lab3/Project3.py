import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from itertools import combinations

G = nx.Graph()
with open('interacting_proteins.txt') as f:
    for line in f:
        p1, p2 = line.strip().split()
        G.add_edge(p1, p2)  

print(f"Loaded graph with {G.number_of_nodes()} proteins and {G.number_of_edges()} interactions")

with open('onco_genes.txt') as f:
    oncogenes = {line.strip() for line in f if line.strip()}

present = set(G.nodes()) & oncogenes
missing = oncogenes - present
oncogenes = list(present)

def avg_shortest_path(nodes):
    dists = []
    for u, v in combinations(nodes, 2):
        try:
            d = nx.shortest_path_length(G, u, v)  
            dists.append(d)
        except nx.NetworkXNoPath:
            pass
    return np.mean(dists)

obs_score = avg_shortest_path(oncogenes)
print(f"Observed avg. shortest‐path among known drivers: {obs_score:.3f}")

n_iter = 1000
null_scores = []
all_nodes = list(G.nodes())
k = len(oncogenes)

for _ in range(n_iter):
    sample = random.sample(all_nodes, k)
    null_scores.append(avg_shortest_path(sample))

plt.hist(null_scores, bins=30, edgecolor='k')
plt.axvline(obs_score, color='r', linewidth=2, label=f'Observed = {obs_score:.2f}')
plt.xlabel('Avg. shortest‐path')
plt.ylabel('Frequency')
plt.legend()
plt.show()

p_val = sum(1 for s in null_scores if s <= obs_score) / n_iter
print(f"Empirical p‐value: {p_val:.3f}")
