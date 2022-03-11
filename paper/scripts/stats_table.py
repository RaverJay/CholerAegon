#!/usr/bin/env python3
# SK

import sys
from glob import glob
import numpy as np
import pandas as pd
import matplotlib
# from cycler import cycler
matplotlib.use('Agg')  # do not require X window
import matplotlib.pyplot as plt
import seaborn as sns


readstats_lr = '/data/mahlzeitlocal/sebastian/cholera_schuldt/assembly_stats/reads_stats/longreads/*.txt'
readstats_sr = '/data/mahlzeitlocal/sebastian/cholera_schuldt/assembly_stats/reads_stats/shortreads/*.txt'
assemblystats = '/data/mahlzeitlocal/sebastian/cholera_schuldt/assembly_stats/general_stats.txt'
coveragestats_lr = '/data/mahlzeitlocal/sebastian/cholera_schuldt/coverage/coverage_longreads/*.mosdepth.summary.txt'
coveragestats_sr = '/data/mahlzeitlocal/sebastian/cholera_schuldt/coverage/coverage_shortreads/*.mosdepth.summary.txt'
fastanidata = '/data/mahlzeitlocal/sebastian/cholera_schuldt/vibrio_nf/results_rebasecall_noconta2_all/fastANI_aggregated_results.tsv'
amr_results = '/data/mahlzeitlocal/sebastian/cholera_schuldt/vibrio_nf/results_rebasecall_noconta2_all/AMR_combined_aggregated_results.csv'

####
# table


data = pd.DataFrame()


####
# readstats

lrstatsfiles = glob(readstats_lr)
srstatsfiles = glob(readstats_sr)


for file in lrstatsfiles:
    with open(file) as infh:
        for line in infh:
            if line.startswith('--- Stats'):
                iso = line.strip('- \n').split(' ')[-1].split('_')[0]
                assert iso.startswith('Iso'), line

            if line.startswith('Num reads:'):
                data.at[iso, 'lr_num_reads'] = int(line.strip().split()[-1].replace(',', ''))
            if line.startswith('Median'):
                data.at[iso, 'lr_median_len'] = float(line.strip().split()[-1].replace(',', ''))
            if line.startswith('Total'):
                data.at[iso, 'lr_total_bases'] = int(line.strip().split()[-1].replace(',', ''))


for file in srstatsfiles:
    with open(file) as infh:
        for line in infh:
            if line.startswith('--- Stats'):
                iso = line.strip('- \n').split(' ')[-1].split('_')[0]
                assert iso.startswith('Iso'), line

            if line.startswith('Num reads:'):
                data.at[iso, 'sr_num_reads'] = int(line.strip().split()[-1].replace(',', ''))
            if line.startswith('Median'):
                data.at[iso, 'sr_median_len'] = float(line.strip().split()[-1].replace(',', ''))
            if line.startswith('Total'):
                data.at[iso, 'sr_total_bases'] = int(line.strip().split()[-1].replace(',', ''))

data.sort_index(inplace=True)
data.index.name = 'isolate'


####
# asm stats

asmtypes = {'pilon_polished': 'hybrid', 'medaka_consensus': 'longreads', 'spades_assembly': 'shortreads'}

with open(assemblystats) as infh:

    iso = None
    for line in infh:
        if line.startswith('--- Stats'):
            iso, asmtype = line.strip('- \n').split(' ')[-1].rsplit('/')[-1].replace('_assembly.fasta', '').rsplit('_', 1)
            assert iso.startswith('Iso'), line
            # asmtype = [asmtypes[asm_p] for asm_p in asmtypes if asm_p == asm_prog][0]
            
        if line.startswith('Num reads:'):
            data.at[iso, f'asm_{asmtype}_numfrag'] = int(line.strip().split()[-1].replace(',', ''))
        if line.startswith('Total'):
            data.at[iso, f'asm_{asmtype}_len'] = int(line.strip().split()[-1].replace(',', ''))


####
# coverage data

lr_covfiles = glob(coveragestats_lr)
sr_covfiles = glob(coveragestats_sr)

for file in lr_covfiles:

    iso = file.rsplit('/', 1)[-1].split('_',1)[0]
    assert iso.startswith('Iso'), iso

    with open(file) as infh:
        for line in infh:
            if line.startswith('total'):
                _, length, bases, mean, mincov, maxcov = line.strip().split()
                data.at[iso, 'lr_meancoverage'] = float(mean)

for file in sr_covfiles:

    iso = file.rsplit('/', 1)[-1].split('_',1)[0]
    assert iso.startswith('Iso'), iso

    with open(file) as infh:
        for line in infh:
            if line.startswith('total'):
                _, length, bases, mean, mincov, maxcov = line.strip().split()
                data.at[iso, 'sr_meancoverage'] = float(mean)




####
# fastani results data

anidata = pd.read_csv(fastanidata, sep='\t', index_col=0, header=None)
anidata.columns = ['reference', 'ANI', 'frags_aligned', 'frags_total']
anidata['sample'] = anidata.index
anidata['asmtype'] = [i.rsplit('_',2)[-2] for i in anidata.index]
anidata.index = [i.rsplit('_',2)[-0] for i in anidata.index]

for asmt in ['longreads', 'hybrid', 'shortreads']:
    data[f'asm_{asmt}_ANI'] = anidata.loc[anidata['asmtype']==asmt, 'ANI']


####
# amr results

amrdata = pd.read_csv(amr_results, sep='\t', index_col=0, header=0)
amrdata['sample'] = amrdata.index
amrdata['asmtype'] = [i.rsplit('_',1)[-1] for i in amrdata.index]
amrdata.index = [i.split('_',1)[0] for i in amrdata.index]

for asmt in ['longreads', 'hybrid', 'shortreads']:
    data[f'asm_{asmt}_numgenesfound'] = amrdata.loc[amrdata['asmtype']==asmt, 'NUM_FOUND']


data['lr_median_len'] = data['lr_median_len'].astype(int)
data['lr_num_reads'] = data['lr_num_reads'].astype(int)
data['sr_num_reads'] = data['sr_num_reads'].astype(int)

data['asm_hybrid_len'] = data['asm_hybrid_len'].astype(int)
data['asm_longreads_len'] = data['asm_longreads_len'].astype(int)
data['asm_shortreads_len'] = data['asm_shortreads_len'].astype(int)

data['asm_hybrid_numfrag'] = data['asm_hybrid_numfrag'].astype(int)
data['asm_longreads_numfrag'] = data['asm_longreads_numfrag'].astype(int)
data['asm_shortreads_numfrag'] = data['asm_shortreads_numfrag'].astype(int)


print(data)
print(data.columns)


############


data.to_csv('all_stats.csv', sep='\t')

# write tex table
data.to_latex('all_stats_latex.txt')



####
# plots

sns.set_theme(style="ticks")

fig, axes = plt.subplots(nrows=1, ncols=5, figsize=(12,4))


data1 = data[['sr_num_reads', 'lr_num_reads', 'lr_median_len']]
ax = axes[0]
# data1.boxplot(ax=ax)
data1.melt()
sns.boxplot(data=data1, ax=ax)
ax.set_yscale('log')
ax.set_title('Read statistics')
# ax.set_xlabel('% of data in subsampling')
ax.set_ylabel('# reads / long read length')
ax.set_xticklabels(['#SR', '#LR', 'LenLR'])


data2 = data[['sr_total_bases', 'lr_total_bases']]
ax = axes[1]
# data2.boxplot(ax=ax)
data2.melt()
sns.boxplot(data=data2, ax=ax)
# ax.set_yscale('log')
ax.set_title('Total bases')
# ax.set_xlabel('% of data in subsampling')
ax.set_ylabel('# bases [Gb]')
ax.set_xticklabels(['SR', 'LR'])


data3 = data[['asm_shortreads_len', 'asm_longreads_len', 'asm_hybrid_len']]
ax = axes[2]
# data3.boxplot(ax=ax)
data3.melt()
sns.boxplot(data=data3, ax=ax)
# ax.set_yscale('log')
ax.set_title('Assembly lengths')
# ax.set_xlabel('% of data in subsampling')
ax.set_ylabel('length [Mb]')
ax.set_xticklabels(['SR', 'LR', 'hybrid'])


data4 = data[['sr_meancoverage', 'lr_meancoverage']]
ax = axes[3]
# data4.boxplot(ax=ax)
data4.melt()
sns.boxplot(data=data4, ax=ax)
# ax.set_yscale('log')
ax.set_title('Assembly coverage')
# ax.set_xlabel('% of data in subsampling')
ax.set_ylabel('mean coverage')
ax.set_xticklabels(['SR', 'LR'])


data5 = data[['asm_shortreads_ANI', 'asm_longreads_ANI', 'asm_hybrid_ANI']]
# filter other serotypes
keep = [(row not in ['Iso02538', 'Iso02539', 'Iso02583']) for row in data5.index]
data5 = data5[keep]
ax = axes[4]
# data5.boxplot(ax=ax)
data5.melt()
sns.boxplot(data=data5, ax=ax)
# ax.set_yscale('log')
ax.set_title('Assembly ANI %')
# ax.set_xlabel('% of data in subsampling')
ax.set_ylabel('% ANI')
ax.set_xticklabels(['SR', 'LR', 'hybrid'])


plt.tight_layout()
plt.savefig('stats_overall.pdf', bbox_inches='tight')
plt.savefig('stats_overall.png', bbox_inches='tight')

plt.close()



table = pd.DataFrame(index=['longreads', 'shortreads', 'hybrid'])

table['mean_num_reads'] = [data['lr_num_reads'].mean(), data['sr_num_reads'].mean(), '-']
table['median_len_reads'] = [data['lr_median_len'].mean(), data['sr_median_len'].mean(), '-']
table['mean_total_bases'] = [data['lr_total_bases'].mean(), data['sr_total_bases'].mean(), '-']
table['mean_coverage'] = [data['lr_meancoverage'].median(), data['sr_meancoverage'].median(), '-']
table['mean_asm_len'] = [data['asm_longreads_len'].mean(), data['asm_shortreads_len'].mean(), data['asm_hybrid_len'].mean()]
table['median_fragments'] = [data['asm_longreads_numfrag'].median(), data['asm_shortreads_numfrag'].median(), data['asm_hybrid_numfrag'].median()]
table['mean_ANI'] = [data['asm_longreads_ANI'].median(), data['asm_shortreads_ANI'].median(), data['asm_hybrid_ANI'].median()]

table = table.T
table = table[['shortreads', 'longreads', 'hybrid']]
print(table)

table.to_csv('table_stats_overall.csv')
table.to_latex('table_stats_overall.tex')

