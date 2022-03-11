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
data['ISO'], data['TIME'] = zip(*[f.split('_') for f in data['SAMPLE']])
data['TIME'] = data['TIME'].apply(lambda x: int(x.strip('min')))

othercols = ['NUM_FOUND', 'SAMPLE', 'ISO', 'TIME']
genecols = [col for col in data.columns if col not in othercols]

print(data)


####


# ground truth for Iso02507
gt_genes = ["APH(3'')-Ib", "APH(6)-Id", "CRP", "Vibrio cholerae varG", "almG", "catB9", "dfrA1", "floR", "sul2", "rsmA"]
skip_genes = ["Escherichia coli parE conferring resistance to fluoroquinolones"]
cols_corr = [g for g in genecols if g in gt_genes]
print(cols_corr)
cols_false = [g for g in genecols if (g not in gt_genes and g not in skip_genes)]
print(cols_false)

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

pdata = vals[['TIME', '# genes correct', '# genes false']].copy()
pdata.columns = ['time after start of sequencing [min]', 'correct', 'false']
pdata['sample'] = vals.index

pdata.index = pdata['time after start of sequencing [min]']
pdata.drop('time after start of sequencing [min]', axis=1, inplace=True)

# pdata = pdata.melt(id_vars=['sample', 'time after start of sequencing [min]'], value_vars=['true positive', 'false positive'], var_name='correctness', value_name='# genes')

print(pdata)

fig, ax = plt.subplots(figsize=(6,4))

pdata.plot(ax=ax, marker=".")
# plt.plot(x=pdata['time after start of sequencing [min]'], scalex='log')
# sns.boxplot(data=pdata, x='time after start of sequencing [min]', y='# genes', hue='correctness')
# g = sns.relplot(data=pdata, x='time after start of sequencing [min]', y='# genes', hue="correctness") #, size="mass", sizes=(10, 200))

# ax.set_title('AMR genes by sequencing time (Iso02507)')

# ax.set_ylabel('# genes')
# for i, cb in enumerate(novels.index):
#     total = len(tpm[tpm['tissues_expressed'] == cb].loc[:,'is_novel'])
#     plt.text(i-0.07, 1.02, total, rotation='vertical')
# xlabs = [lab.lstrip('0') if lab != '100' else '100' for lab in [t.get_text() for t in ax.get_xticklabels()]]
# ax.set_xticklabels(xlabs)
ax.set_xscale('log')
# g.set(xscale="log")

ax.set_ylim([-0.5,10.5])

fig.tight_layout()
plt.savefig('genes_found_timeline_multi.pdf', bbox_inches='tight')
plt.savefig('genes_found_timeline_multi.png', bbox_inches='tight')


######
# second axis

if len(sys.argv) <= 2:
    exit(0)


compl = pd.read_csv(sys.argv[2], sep='\t')
compl['SAMPLE'] = [row.split('/samples/',1)[-1].rsplit('/longreads/')[0] for row in compl['File']]
compl.index = compl['SAMPLE']
compl['ISO'], compl['TIME'] = zip(*[f.split('_') for f in compl['SAMPLE']])
compl['TIME'] = compl['TIME'].apply(lambda x: int(x.strip('min')))
print(compl)



col = 'red'
ax2 = ax.twinx()
ax2.set_xscale('log')

# plot it to other ax to ensure same y alignment
ax.plot(compl['TIME'], compl['completeness']*10, color=col, marker='.', label='completeness')


ax2.set_ylabel('% assembly completeness', color=col) 
ax2.tick_params(axis='y', labelcolor=col)
ax2.set_ylim([-0.05, 1.05])
ax2.set_yticklabels(['PAD', '0', '20', '40', '60', '80', '100'])

ax.set_ylabel('# detected AMR genes', color='#1f77b4')
ax.tick_params(axis='y', labelcolor='#1f77b4')
# ax2.set_xticklabels(xlabs)

ax.legend()

plt.savefig('genes_found_timeline_multi_twinx.pdf', bbox_inches='tight')
plt.savefig('genes_found_timeline_multi_twinx.png', bbox_inches='tight')






