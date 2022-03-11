#!/usr/bin/env python3
# SK

import sys
import argparse
import pandas as pd


def error(string, error_type=1):
    sys.stderr.write(f'ERROR: {string}\nAborting.\n')
    sys.exit(error_type)


def log(string, newline_before=False):
    if newline_before:
        sys.stderr.write('\n')
    sys.stderr.write(f'LOG: {string}\n')

###

log('Started combine_results.py ...')

default_coverage_threshold = 80.
default_identity_threshold = 80.

### parse args
parser = argparse.ArgumentParser(description='Combine results from RGI and abricate')
parser.add_argument("--rgi_results", help="rgi results files", required=True)
parser.add_argument("--abricate_results", help="abricate results files", required=True)
parser.add_argument("--output", help="output file prefix", required=True)
parser.add_argument("--coverage_threshold", help="minimum coverage% to keep", default=default_coverage_threshold)
parser.add_argument("--identity_threshold", help="minimum identity% to keep", default=default_identity_threshold)
args = parser.parse_args()



coverage_threshold = float(args.coverage_threshold)
identity_threshold = float(args.identity_threshold)



### build table
# with consistent column names
column_names = ['File', 'Gene', 'Coverage', 'Identity', 'Model_type', 'Prediction_tool', 'Contig', 'Start', 'Stop', 'Strand', 'Notes']


### parse RGI results

# ORF_ID	Contig	Start	Stop	Orientation	Cut_Off	Pass_Bitscore	Best_Hit_Bitscore	Best_Hit_ARO	Best_Identities
# ARO	Model_type	SNPs_in_Best_Hit_ARO	Other_SNPs	Drug Class	Resistance Mechanism	AMR Gene Family	Predicted_DNA
# 	Predicted_Protein	CARD_Protein_Sequence	Percentage Length of Reference Sequence	ID	Model_ID	Nudged	Note
df = pd.read_csv(args.rgi_results, header=0, sep='\t')

# get filename
basename = args.rgi_results.rsplit('/', 1)[-1]
assert basename.startswith('RGI_')
assert basename.endswith('.txt')
samplename = basename.split('RGI_',1)[-1].rsplit('.txt',1)[0]

# set sample name
df['FILE'] = samplename

# set fixed columns
df['Prediction_tool'] = 'RGI'
df['Notes'] = ''

# select columns from results and rename
df = df[['FILE', 'Best_Hit_ARO', 'Percentage Length of Reference Sequence', 'Best_Identities', 'Model_type', 'Prediction_tool', 'Contig', 'Start', 'Stop', 'Orientation', 'Notes']]
df.columns = column_names

data = df.copy()
print(data)


### parse and merge Abricate results

# #FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
df = pd.read_csv(args.abricate_results, header=0, sep='\t')

# check if there is Abricate data
if df.shape[0] > 0:

    # get filename
    df['#FILE'] = df['#FILE'].iloc[0].rsplit('.fasta',1)[0].rsplit('_assembly',1)[0]

    # set fixed columns
    df['Model_type'] = 'protein homolog model'
    df['Prediction_tool'] = 'Abricate_CARD'
    df['Notes'] = ''

    # select columns from results and rename
    df = df[['#FILE', 'GENE', '%COVERAGE', '%IDENTITY', 'Model_type', 'Prediction_tool', 'SEQUENCE', 'START', 'END', 'STRAND', 'Notes']]
    df.columns = column_names

    # CARE: Abricate changes spaces to underscores in gene names, revert that
    df['Gene'] = [g.replace('_',' ') for g in df['Gene']]

    print(df)


    ### merge
    # Strategy:
    # - Leave genes from gene variant models as is (only RGI uses these)
    # - For each gene from protein homology models that has Abricate and RGI results, use the Abricate results if they have equal/higher coverage
    #       (they usually have, as Abricate merges hits much more reliably)

    visited_genes = []

    for row in df.index:

        gene = df.at[row, 'Gene']

        # handled already?
        if gene in visited_genes:
            continue
        visited_genes.append(gene)

        # not homology model? (this is redundant atm, as only RGI uses other models)
        if df.at[row, 'Model_type'] != 'protein homolog model':
            continue


        # handle Abricate data
        abr_subset = df[df['Gene']==gene].copy()

        # is there RGI data? maybe replace
        if (data['Gene']==gene).any():

            rgi_subset = data[data['Gene']==gene].copy()
            max_rgi_cov = max(rgi_subset['Coverage'])

            max_abr_cov = max(abr_subset['Coverage'])
            
            if max_abr_cov >= max_rgi_cov:
                # replace
                abr_subset['Notes'] += 'Replaced RGI results with Abricate results;'
                data.drop(list(rgi_subset.index), inplace=True)
                data = pd.concat([data, abr_subset], ignore_index=True)

        else:
            # only Abricate found this, add it
            data = pd.concat([data, abr_subset], ignore_index=True)

        

data.sort_values(by='Gene', inplace=True)
print(data)



### filter
# Strategy: (default values)
# - reject hits under 80% reference coverage
# - reject hits under 80% reference identity


data['Coverage_passed'] = data['Coverage'] >= coverage_threshold
data['Identity_passed'] = data['Identity'] >= identity_threshold
data.loc[~data['Coverage_passed'], ['Notes']] += f'Failed coverage threshold {coverage_threshold};'
data.loc[~data['Identity_passed'], ['Notes']] += f'Failed identity threshold {identity_threshold};'
data.to_csv(f'{args.output}_prefilter.csv', sep='\t', index=False)


data_pass = data[ data['Coverage_passed'] & data['Identity_passed'] ]
data_pass.to_csv(f'{args.output}.csv', sep='\t', index=False)









