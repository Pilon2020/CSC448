import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import hypergeom, norm
from tqdm import tqdm
import pandas as pd

# Parameters
total_genes = 6000
colored_green_total = 4260
sample_size = 50       # <-- this is the 'N' for hypergeom
observed_green = 38

def simulate_draws(n_simulations):
    return hypergeom.rvs(
        M=total_genes,
        n=colored_green_total,
        N=sample_size,
        size=n_simulations
    )


trial_sizes = np.linspace(10, 10000, 1000, dtype=int)


# run sims
means, variances = [], []
for X in tqdm(trial_sizes, desc="Stability sims"):
    res = simulate_draws(X)
    means.append(np.mean(res >= observed_green))
    variances.append(np.var(res >= observed_green))
means = np.array(means)
variances = np.array(variances)

# polynomial fit for visualization
degree = 5
log_x = np.log10(trial_sizes)
poly = np.poly1d(np.polyfit(log_x, means, degree))
x_cont = np.linspace(trial_sizes.min(), trial_sizes.max(), 200)
y_cont = poly(np.log10(x_cont))

plt.figure()
plt.plot(trial_sizes, means, 'o', label='Estimated p-value')
plt.plot(x_cont, y_cont, '--', label=f'Poly deg {degree}')
plt.axhline(y=y_cont[-1], color='r', label=f'Final y = {y_cont[-1]:.4f}')
plt.xscale('log')
plt.xlabel('Number of simulations')
plt.ylabel('Estimated p-value')
plt.title('Stabilization with polynomial fit on log(X)')
plt.legend()
plt.grid(True)


threshold = 1e-4
diffs = np.abs(np.diff(means))
stable_idxs = np.where(diffs < threshold)[0]
if stable_idxs.size > 0:
    chosen_X = trial_sizes[stable_idxs[0] + 1]
    print(f"Stability reached at X = {chosen_X} (Δp < {threshold})")
else:
    chosen_X = trial_sizes[-1]
    print(f"Threshold never reached; using max trials: {chosen_X}")

final_res = simulate_draws(chosen_X)

# z‐score for 95% CI
z = norm.ppf(0.975)

for obs in (35, 38, 40):

    p_hat = np.mean(final_res >= obs)           # https://en.wikipedia.org/wiki/Population_proportion
    se    = np.sqrt(p_hat*(1-p_hat)/chosen_X)   # https://external-content.duckduckgo.com/iu/?u=https%3A%2F%2Fgetcalc.com%2Fformula%2Fstatistics%2Fstandard-error-observed-sample-proportion.png&f=1&nofb=1&ipt=cabe5d4ed51b3b4fac06311b18c545b97c0cb31d32312b0f0637b5503f5c902a
    ci_lo = p_hat - z*se
    ci_hi = p_hat + z*se
    if obs == 38:
        p1 = p_hat
    else: pass    
    
    exact_p = hypergeom.sf(
        obs - 1,
        total_genes,
        colored_green_total,
        sample_size
    )

    sig = "YES" if exact_p < 0.05 else "NO"
    print(f"\n=== Observed = {obs} greens ===")
    print(f"Simulated P-Value = {p_hat:.4f} ± {se:.4f}")
    print(f"95% CI      = [{ci_lo:.4f}, {ci_hi:.4f}]")
    print(f"Exact P-Value= {exact_p:.6f}")
# plt.savefig('Lab2/Pvalue_regression.png', bbox_inches='tight')
# plt.show()


# =================================================================================
# Question 2.4
N = 50 # number of items in list
M = 6000  # number of items in background
pvalue = [0.003, 0.01] # P-value required range
step = 10

ranges = []
for n in range(50, M+1, step):
    for k in range(0,min(n,N)+1):
        p = hypergeom.sf(k-1, M, n, N)
        if 0.003 <= p <= 0.01:
            ranges.append((n, k, p))
            
df = pd.DataFrame(ranges, columns=['n_background', 'k_list', 'p_value'])
df = df.sort_values(by='p_value').reset_index(drop=True)
df.to_csv('Lab2/nk_probelm.txt', index=False, sep='\t')
print(f'\n{df.iloc[1]}')


# =================================================================================
# Question 2.5
p1 = 0.2703
p2 = 0.003009
pvals = {'Green': p1, 'Purple': p2}
alpha = 0.05
m = 10

#H0 = aplha_g = alpha_p = alpha_n ...
# Need to multiply pvalues by m #https://middleprofessor.github.io/benchbiostats/chapters/multiple-testing.html
print("=== Bonferroni Correction ===") # https://en.wikipedia.org/wiki/Bonferroni_correction
alpha_bonf = alpha / m
# Reject the null hypothesis if pval is less than alpha_bonf
print(alpha_bonf)
for name, p in pvals.items():
    newp = p * 10
    sig = newp < alpha_bonf
    print(f'{name:6s}: old p = {p:.4f}, new p = {newp:.4f}, still significant? {sig}')
    
    
# Try using the Benjamini-Hochberg method
