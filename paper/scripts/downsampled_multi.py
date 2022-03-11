#!/usr/bin/env python3
# SK

import sys
import numpy as np
import pandas as pd
import matplotlib
# from cycler import cycler
matplotlib.use('Agg')  # do not require X window
import matplotlib.pyplot as plt
import seaborn as sns


aggdata = sys.argv[1]



####
# aggregated results data

data = pd.read_csv(aggdata, sep='\t', index_col=0)
data = data.sort_index()


data['SAMPLE'] = [f.rsplit('_',1)[0] for f in data.index]
data['ISO'], data['PERCENT'], data['REPLICATE'] = zip(*[f.split('_') for f in data['SAMPLE']])
data['PERCENT'] = data['PERCENT'].apply(lambda x: x.strip('pct'))

othercols = ['NUM_FOUND', 'SAMPLE', 'ISO', 'PERCENT', 'REPLICATE']
genecols = [col for col in data.columns if col not in othercols]

print(data)


####
# get stats

# ground truth for Iso02507
gt_genes = ["APH(3'')-Ib", "APH(6)-Id", "CRP", "Vibrio cholerae varG", "almG", "catB9", "dfrA1", "floR", "sul2", "rsmA"]
skip_genes = ["Escherichia coli parE conferring resistance to fluoroquinolones"]
cols_corr = [g for g in genecols if g in gt_genes]
cols_false = [g for g in genecols if (g not in gt_genes and g not in skip_genes)]

vals = data.copy()



for row in vals.index:

    # vals.at[row, '# genes found'] = data.loc[data['PERCENT'] == row, ['NUM_FOUND']].iat[0,0]
    vals['# genes found'] = vals['NUM_FOUND']

    # true / false positives
    vals.at[row, '# genes correct'] = (vals.loc[[row], cols_corr] != '.').sum(axis=1).iloc[0]
    vals.at[row, '# genes false'] = (vals.loc[[row], cols_false] != '.').sum(axis=1).iloc[0]

    # # coverage of genes
    # # take each copy separately
    # pct_list = []
    # for field in data.loc[data['PERCENT'] == row, genecols].iloc[0]:
    #     if field == '.':
    #         continue
    #     pct_list.extend(field.split(';'))
    # print(pct_list)
    # vals.at[row, 'mean % coverage'] = np.mean([float(p) for p in pct_list]) / 100

print(vals)

pdata = vals[['PERCENT', '# genes correct', '# genes false']].copy()
pdata.columns = ['% of data subsampled', 'correct', 'false']
pdata['sample'] = vals.index

pdata = pdata.melt(id_vars=['sample', '% of data subsampled'], value_vars=['correct', 'false'], var_name='correctness', value_name='# genes')

# crop to 10%
maxpercent = 10
keep = [int(p)<=maxpercent for p in pdata['% of data subsampled']]
pdata = pdata[keep]

print(pdata)

fig, ax = plt.subplots(figsize=(6,4))

sns.barplot(data=pdata, x='% of data subsampled', y='# genes', hue='correctness')

# ax.set_title('AMR genes found in subsampled data (Iso02507)')

ax.set_ylim([0,10.5])
ax.set_xlabel('% of data / total bases [Mb] / coverage')
ax.set_ylabel('# detected AMR genes')
# for i, cb in enumerate(novels.index):
#     total = len(tpm[tpm['tissues_expressed'] == cb].loc[:,'is_novel'])
#     plt.text(i-0.07, 1.02, total, rotation='vertical')

total_bases = 864_385_799
asm_length = 4_144_502

xlab_perc = [lab.lstrip('0') if lab != '100' else '100' for lab in [t.get_text() for t in ax.get_xticklabels()]]
xlab_bases = [(float(p)/100*total_bases) for p in xlab_perc]
xlab_cover = [tb/asm_length for tb in xlab_bases]


xlabs = [f'{p}\n{b/1e6:.1f}\n{c:.0f}X' for p, b, c in zip(xlab_perc, xlab_bases, xlab_cover)]
ax.set_xticklabels(xlabs)
fig.tight_layout()

plt.savefig('genes_found_subsampled_multi.pdf', bbox_inches='tight')
plt.savefig('genes_found_subsampled_multi.png', bbox_inches='tight')


######
# second axis

if len(sys.argv) <= 2:
    exit(0)


compl = pd.read_csv(sys.argv[2], sep='\t')
compl['SAMPLE'] = [row.split('/samples/',1)[-1].rsplit('/longreads/')[0] for row in compl['File']]
compl.index = compl['SAMPLE']
compl['ISO'], compl['PERCENT'], compl['REPLICATE'] = zip(*[f.split('_') for f in compl['SAMPLE']])
compl['PERCENT'] = compl['PERCENT'].apply(lambda x: x.strip('pct'))
print(compl)

cmeans = compl.groupby(['PERCENT']).mean()
cmeans['percent_int'] = [int(p) for p in cmeans.index]
keep = cmeans['percent_int'] <= maxpercent
cmeans = cmeans[keep]
print(cmeans)


col = 'red'
ax2 = ax.twinx()
ax2.plot(cmeans['completeness'], color=col, marker='.')


ax2.set_ylabel('% assembly completeness', color=col) 
ax2.tick_params(axis='y', labelcolor=col)
ax2.set_ylim([0, 1.05])
ax2.set_yticklabels(['0', '20', '40', '60', '80', '100'])

ax.set_ylabel('# detected AMR genes', color='#1f77b4')
ax.tick_params(axis='y', labelcolor='#1f77b4')
ax2.set_xticklabels(xlabs)

ax.legend(loc='lower right')

plt.savefig('genes_found_subsampled_multi_twinx.pdf', bbox_inches='tight')
plt.savefig('genes_found_subsampled_multi_twinx.png', bbox_inches='tight')




