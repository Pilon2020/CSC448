"""

CSC448 - Bioinformatics
Project 2. Wine making yeast Part 1
@Author: Oliver Pilon
Date: April 22, 2025

"""

# Library Imports
import numpy as np      #NumPy Documentation Link: https://numpy.org/doc/2.2/numpy-user.pdf       
import scipy as sp      #SciPy Documentation Link: https://docs.scipy.org/doc/scipy/tutorial/index.html#user-guide
import sklearn as sk    #sklearn Documentation Link: https://scikit-learn.org/stable/user_guide.html
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from enum import IntEnum

df = pd.read_csv('Lab2/data/diauxic_raw_ratios.txt',sep='\t',usecols=[1,3,4,5,6,7,8,9])


names = df['Name'].tolist()
ratios = df.columns[1:].tolist()        # ['R1.Ratio', …, 'R7.Ratio']
data   = df.iloc[:, 1:].astype(float).to_numpy()


fig, ax = plt.subplots()
im = ax.imshow(data, aspect='auto')

# X‐axis: your NAMES across the bottom
ax.set_xticks(np.arange(data.shape[1]))

# Y‐axis: your RATIO titles on the left
ax.set_yticks(np.arange(data.shape[0]))

ax.set_title(f"Heat Map")
fig.tight_layout()
plt.show()
