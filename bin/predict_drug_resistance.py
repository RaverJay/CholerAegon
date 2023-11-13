#!/usr/bin/env python3
# SK


import sys
import argparse
import pandas as pd
import obonet

def error(string, error_type=1):
    sys.stderr.write(f'ERROR: {string}\n')
    sys.exit(error_type)

def log(string, newline_before=False):
    if newline_before:
        sys.stderr.write('\n')
    sys.stderr.write(f'LOG: {string}\n')


### parse args

parser = argparse.ArgumentParser(description='Predict resistances from detected AMR genes')
# parser.add_argument("--sampleslist", help="path to big samples list csv")
parser.add_argument("--aro_ontology", required=True, help="path to aro.obo")
parser.add_argument("--gene_results", required=True, help="aggregated gene results tsv")
parser.add_argument("--output", help='output tsv file', default='resistance_prediction.tsv')

args = parser.parse_args()




### load CARD

# ontology obo
obo = obonet.read_obo(args.aro_ontology)
log(f'Loaded ARO ontology from file {args.aro_ontology}.')


# extract ARO accessions and names
accs_dict = {}
names_dict = {}

for node in obo.nodes(data=True):
    gene_acc = node[0]
    name = node[1]['name']

    accs_dict[name] = gene_acc
    names_dict[gene_acc] = name

def get_acc(name):
    return accs_dict[name]
def get_name(acc):
    return names_dict[acc]



### data

data = pd.read_csv(args.gene_results, sep='\t', index_col='#FILE')
      

def stringify(strlst):
    return ', '.join(strlst)

log(f'Loaded gene detection results.')
log(f'- List of files: {stringify(data.index.tolist())}')


cols_other = ['NUM_FOUND']
cols_amr = [col for col in data.columns if col not in cols_other]

data = data[cols_other + cols_amr].sort_index()



### predictions

###
# WHO suggested antibiotics

antibiotics_primary = ['doxycycline']
antibiotics_secondary = ['azithromycin', 'ciprofloxacin']
who_antibiotics = antibiotics_primary + antibiotics_secondary



if cols_amr == []:
    log('No AMR genes were found.')
else:

    for row in data.index:
        # print(row)
        

        ## get genes and check for multiple component genes
        resistance_genes_acc = []
        multi_genes_acc = {}
        
        
        for gene in cols_amr:
            if data.loc[row, gene] != '.':

                gene_acc = get_acc(gene)
                resistance_genes_acc.append(gene_acc)

                # find potential multi genes by
                # checking edges to parent nodes
                for child, parent, key in obo.out_edges(gene_acc, keys=True):

                    if key == 'part_of':
                        
                        # print(f'• {get_name(child)} ⟶ {key} ⟶ {get_name(parent)}')
                        multi_genes_acc[parent] = []


        # check multigene candidates
        for multi_gene_acc in multi_genes_acc:
            parts_acc = []
            has_all = True

            for child, parent, key in obo.in_edges(multi_gene_acc, keys=True):

                if key == 'part_of':

                    # print(f'• {get_name(child)} ⟶ {key} ⟶ {get_name(parent)}')
                    if child in resistance_genes_acc:
                        parts_acc.append(child)
                    else:
                        has_all = False
                        break
            if has_all:
                # has all parts of multi-gene
                # only add and consider further if complete
                multi_genes_acc[multi_gene_acc] = parts_acc
                if multi_gene_acc not in resistance_genes_acc:
                    resistance_genes_acc.append(multi_gene_acc)


        
        ## collect resistances

        resisted_drugs_acc = {}
        resisted_drug_classes_acc = {}

        for gene_acc in resistance_genes_acc:

            # find potential resistances by
            # checking edges to parent nodes
            for child, parent, key in obo.out_edges(gene_acc, keys=True):

                if key == 'confers_resistance_to_antibiotic':
                    
                    # print(f'• {get_name(child)} ⟶ {key} ⟶ {get_name(parent)}')
                    if parent not in resisted_drugs_acc:
                        resisted_drugs_acc[parent] = []
                    if child in multi_genes_acc:
                        resisted_drugs_acc[parent] += [get_name(child)+':'+('+'.join(sorted([get_name(acc) for acc in multi_genes_acc[child]])))]
                    else:
                        resisted_drugs_acc[parent] += [get_name(child)]

                elif key == 'confers_resistance_to_drug_class':

                    if parent not in resisted_drug_classes_acc:
                        resisted_drug_classes_acc[parent] = []
                    if child in multi_genes_acc:
                        resisted_drug_classes_acc[parent] += [get_name(child)+':'+('+'.join(sorted([get_name(acc) for acc in multi_genes_acc[child]])))]
                    else:
                        resisted_drug_classes_acc[parent] += [get_name(child)]



        ## check for multiple component drugs
        multi_drugs_acc = []
        resisted_multi_drugs_acc = {}

        # get all candidates
        for drug_acc in resisted_drugs_acc:
            for child, parent, key in obo.in_edges(drug_acc, keys=True):
                if key == 'has_part':
                    # print(f'• {get_name(child)} ⟶ {key} ⟶ {get_name(parent)}')
                    if child not in multi_drugs_acc:
                        multi_drugs_acc.append(child)

        # check them
        for multi_drug_acc in multi_drugs_acc:
            parts_acc = []
            for child, parent, key in obo.out_edges(multi_drug_acc, keys=True):
                if key == 'has_part':
                    parts_acc.append(parent)

            if parts_acc != []:
                resist_all = True
                for part_acc in parts_acc:
                    if part_acc not in resisted_drugs_acc:
                        resist_all = False
                        break

                if resist_all:
                    # resistant against all parts
                    if multi_drug_acc not in resisted_multi_drugs_acc:
                        resisted_multi_drugs_acc[multi_drug_acc] = []
                    resisted_multi_drugs_acc[multi_drug_acc] += ['+'.join([get_name(part)+('('+','.join(resisted_drugs_acc[part]))+')' for part in parts_acc])]

        # add
        for multi_drug_acc in resisted_multi_drugs_acc:
            if multi_drug_acc not in resisted_drugs_acc:
                resisted_drugs_acc[multi_drug_acc] = []
            resisted_drugs_acc[multi_drug_acc] += resisted_multi_drugs_acc[multi_drug_acc]

        print(resisted_drugs_acc)
        print(resisted_drug_classes_acc)


        ## check antibiotics

        who_antibs_resisted = {}


        for antib in who_antibiotics:

            # get ARO accession
            antib_acc = get_acc(antib)


            # drug class
            for child, parent, key in obo.out_edges(antib_acc, keys=True):
                if key == 'is_a':
                    if parent in resisted_drug_classes_acc:
                        # print(f'-- resistant to {get_name(parent)} ({parent})')
                        
                        # add antib and gene
                        if antib not in who_antibs_resisted:
                            who_antibs_resisted[antib_acc] = []
                        who_antibs_resisted[antib_acc] += [get_name(child)]
            
            
            # drug
            if antib_acc in resisted_drugs_acc:
                # print(f'-- resistant to {antib} ({antib_acc})')
                if antib not in who_antibs_resisted:
                    who_antibs_resisted[antib_acc] = []
                for gene in resisted_drugs_acc[antib_acc]:
                    if get_acc(gene) not in who_antibs_resisted[antib_acc]:
                        who_antibs_resisted[antib_acc] += [gene]

        
        print(who_antibs_resisted)


        # write to datatable
        drug_res_str = ';'.join([f'{get_name(drug)} ({",".join(resisted_drugs_acc[drug])})' for drug in resisted_drugs_acc])
        data.loc[row, 'DRUG_RESISTANCES'] = '.' if drug_res_str == '' else drug_res_str

        drug_class_res_str = ';'.join([f'{get_name(drugc)} ({",".join(resisted_drug_classes_acc[drugc])})' for drugc in resisted_drug_classes_acc])
        data.loc[row, 'DRUG_CLASS_RESISTANCES'] = '.' if drug_class_res_str == '' else drug_class_res_str

        who_drug_res_str = ';'.join([f'{get_name(drug)} ({",".join(who_antibs_resisted[drug])})' for drug in who_antibs_resisted])
        data.loc[row, 'WHO_SUGGESTED_DRUG_RESISTANCES'] = '.' if who_drug_res_str == '' else who_drug_res_str



### save
data.to_csv(args.output, sep='\t')



