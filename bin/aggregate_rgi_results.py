#!/usr/bin/env python3
# SK

import sys
import pandas as pd

infiles = sys.argv[1:]


table = None
samples = []

# add sample column and concatenate data
for file in infiles:

    samplename = file.rsplit('.', 1)[0]
    samples.append(samplename)

    data = pd.read_table(file, header=0)

    # handle empty results
    if data.shape[0] == 0:
        continue

    data.insert(0, 'Sample', samplename)

    if table is None:
        table = data
    else:
        table = pd.concat((table, data), ignore_index=True)

print(samples, file=sys.stderr)
# table.to_csv('table.csv', sep='\t')



# transform to abricate aggregated format
hits = {}
genes = []


for row in table.index.tolist():

    sample = table.at[row, 'Sample']
    gene = table.at[row, 'Best_Hit_ARO']
    pcoverage = table.at[row, 'Percentage Length of Reference Sequence']

    if gene not in genes:
        genes.append(gene)

    if sample not in hits:
        hits[sample] = {gene: [pcoverage]}
    else:
        if gene not in hits[sample]:
            hits[sample][gene] = [pcoverage]
        else:
            hits[sample][gene].append(pcoverage)

print(hits, file=sys.stderr)
print(genes, file=sys.stderr)



agg = pd.DataFrame(index=hits.keys())
agg.index.name = '#FILE'


# include samples with no hits
for sample in samples:

    # num genes found
    agg.at[sample, 'NUM_FOUND'] = len(hits[sample].keys()) if sample in hits else 0

    for gene in genes:

        if (sample in hits) and (gene in hits[sample]):
            agg.at[sample, gene] = ';'.join([str(g) for g in hits[sample][gene]])
        else:
            agg.at[sample, gene] = '.'

agg['NUM_FOUND'] = agg['NUM_FOUND'].astype('int16')
agg.to_csv(sys.stdout, sep='\t')








