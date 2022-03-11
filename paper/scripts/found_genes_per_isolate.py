#!/usr/bin/env python3
# SK

import sys
from matplotlib.pyplot import axis
import pandas as pd


data = pd.read_csv(sys.argv[1], sep='\t')
othercols = ['Assembly number', 'Sample', 'Assembly method', 'FastANI %', '# Genes found']
genecols = [c for c in data.columns if c not in othercols]

data.drop('Assembly number', axis=1).to_latex('found_genes_per_isolate.tex', index=False)


# tex summary
booldf = data[genecols] != '.'
data = pd.concat([data[['Assembly method']], booldf], axis=1)

gb = data.groupby(['Assembly method']).sum()
gb = gb.T
gb = gb[['shortreads', 'longreads', 'hybrid']]

gb.to_latex(sys.stdout)



