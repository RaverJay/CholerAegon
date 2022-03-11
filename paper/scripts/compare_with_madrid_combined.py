#!/usr/bin/env python3

import sys
import pandas as pd
import argparse

from pandas.io.pytables import IndexCol

def error(string, error_type=1):
    sys.stderr.write(f'ERROR: {string}\n')
    sys.exit(error_type)

def log(string, newline_before=False):
    if newline_before:
        sys.stderr.write('\n')
    sys.stderr.write(f'LOG: {string}\n')


### parse args

parser = argparse.ArgumentParser(description='compare resistance predictions to lab results')
# parser.add_argument("--sampleslist", help="path to big samples list csv")
parser.add_argument("--predictions", required=True, help="path to .csv")
parser.add_argument("--madrid_results", required=True, help="path to madrid results with Iso numbers .xlsx")
args = parser.parse_args()


### madrid data

rawdata = pd.read_excel(args.madrid_results, index_col='Iso-Nummer')

data = rawdata.loc[:,('AMP.1', 'C.1', 'CIP.1', 'CN.1', 'TE.1', 'SxT.1', 'NA.1')]
data.columns = ('AMP', 'C', 'CIP', 'CN', 'TE', 'SxT', 'NA')
print(data)

# ARO accessions of the Madrid set:
madrid_antibs = ['ampicillin', 'chloramphenicol', 'ciprofloxacin', 'nalidixic acid', 'gentamicin A', 'tetracycline', 'trimethoprim-sulfamethoxazole']
madrid_aros = ['ARO:3000637', 'ARO:3000385', 'ARO:0000036', 'ARO:3000661', 'ARO:3004015', 'ARO:0000051', 'ARO:3004024']

map_to_aros = {'AMP':'ARO:3000637', 'C':'ARO:3000385', 'CIP': 'ARO:0000036', 'CN': 'ARO:3004015', 'TE': 'ARO:0000051', 'SxT': 'ARO:3004024', 'NA': 'ARO:3000661'}
map_to_names = {'AMP':'ampicillin', 'C':'chloramphenicol', 'CIP': 'ciprofloxacin', 'CN': 'gentamicin A', 'TE': 'tetracycline', 'SxT': 'trimethoprim-sulfamethoxazole', 'NA': 'nalidixic acid'}


### prediction data
      
pred = pd.read_csv(args.predictions, sep='\t')
print(pred)


### compare

rcolumns = []
for c in data.columns:
    rcolumns += [c+'_test', c+'_pred', c+'_comp']
results = pd.DataFrame(columns=rcolumns)


for row in pred.index:
    iso = pred.loc[row, '#FILE'].rsplit('_',1)[0]
    preds = pred.loc[row, 'DRUG_RESISTANCES']

    pred_names = [] if type(preds) is not str else [p.split(' (')[0] for p in preds.split(';')]

    for antib in data.columns:
        results.loc[iso, antib+'_test'] = data.loc[iso, antib]

        pred_bool = map_to_names[antib] in pred_names
        results.loc[iso, antib+'_pred'] = 'yes' if pred_bool else 'no'
        
        comp = 'unhandled'
        if type(data.loc[iso, antib]) is not str:
            comp = 'nodata'
        else:
            comp = data.loc[iso, antib]

        comp += '_' + ('yes' if pred_bool else 'no')
        results.loc[iso, antib+'_comp'] = comp
        
print('RESULTS')
print(results)

###
# resistance genes corresponding to the antibiotics

antibs = {}
genes = {}

for field in pred['DRUG_RESISTANCES']:

    for pair in field.split(';'):

        antib, bracket = pair.split(' (')
        bracket = bracket.rsplit(')',1)[0]

        if '+' not in bracket:
            # single part
            gs = bracket.split(',')
            for gene in gs:
                if antib not in genes:
                    genes[antib] = [gene]
                else:
                    if gene not in genes[antib]:
                        genes[antib] += [gene]

        else:
            # multiple parts
            antibs[antib] = bracket.split('+')
                
print(antibs)
print(genes)



#####
# report

comp = results.loc[:,[c for c in results.columns if 'comp' in c]]


htmlheader = '''<!DOCTYPE html><html><head>
<title>Comparison to Madrid results</title>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">

<style>
* {
    font-family:"Helvetica Neue",Helvetica,"Segoe UI",Arial,freesans,sans-serif
}

.content {
max-width: 1700px;
margin: auto;
}

table.tablestyle {
background-color: #FFFFFF;
width: 1000px;
text-align: center;
border-collapse: collapse;
}
table.tablestyle td, table.tablestyle th {
border: 2px solid #8B8B8B;
padding: 5px 5px;
}
table.tablestyle tbody td {
font-size: 20px;
color: #000000;
}
# table.tablestyle tr:nth-child(even) {
# background: #E6F5FF;
# }
table.tablestyle thead {
background: #E6F5FF;
}
table.tablestyle thead th {
font-size: 20px;
font-weight: bold;
color: #000000;
text-align: center;
}
table.tablestyle tfoot td {
font-size: 13px;
}
table.tablestyle tfoot .links {
text-align: right;
}
table.tablestyle tfoot .links a{
display: inline-block;
background: #FFFFFF;
color: #398AA4;
padding: 2px 8px;
border-radius: 5px;
}
</style>
</head>

<body>
<div class="content">'''


htmlfooter = '''
</div>
</body></html>
'''


color_markup = '#8a006d'
color_good_green = '#046907'
color_warn_orange = '#ac7800'
color_error_red = '#a50500'

pred_corr = ('R_yes', 'SD_yes', 'I_yes', 'S_no')
pred_false = ('R_no', 'SD_no', 'I_no')
pred_falseish = ('S_yes',)



# summary
rownames = comp.index.tolist()

for col in comp.columns:
    corr, false = 0, 0

    for str_in in comp[col]:
        if str_in in pred_corr:
            corr += 1
        elif (str_in in pred_false) or (str_in in pred_falseish):
            false += 1
    
    total = corr + false

    comp.loc['Prediction correctness', col] = f'<b><font color="{color_good_green}">{corr/total*100:.0f}%</font> ' + \
        f'<font color="{color_error_red}">{false/total*100:.0f}%</font></b>'

    comp.loc['Correct', col] = f'<b><font color="{color_good_green}">{corr}</font>'
    comp.loc['False', col] = f'<b><font color="{color_error_red}">{false}</font>'



comp = comp.reindex(['Prediction correctness', 'Correct', 'False'] + rownames)


# markup
def formatter(str_in):
    color = None
    if str_in in pred_corr:
        color = color_good_green
    elif str_in in pred_false:
        color = color_error_red
    elif str_in in pred_falseish:
        color = color_markup

    return f'<b><font color="{color}">{str_in}</font></b>'

print('COMP')
print(comp)


# tex markup
def texcolor(str_in):
    color = None
    if str_in in pred_corr:
        return '\color{green}{' + str_in + '}'
    elif str_in in pred_false:
        return '\color{red}{' + str_in + '}'
    elif str_in in pred_falseish:
        return '\color{purple}{' + str_in + '}'

    return str_in



# csv data
csvcomp = comp[3:]
csvcomp.to_csv('comparison_madrid_result.csv')
csvcomp.to_latex('comparison_madrid_result.tex')




csvpred = pred.copy()
csvpred.index = [row.rsplit('_',1)[0] for row in csvpred['#FILE']]

csvdata = pd.concat((csvcomp, csvpred.loc[:,pred.columns[1:]]), axis=1)
csvdata.to_csv('comparison_madrid_result_full.csv')



comp.columns = [c.rsplit('_comp',1)[0] + '\n' + map_to_names[c.rsplit('_comp',1)[0]] for c in comp.columns]
comp[1:] = comp[1:].applymap(formatter)


### genes row
for col in comp.columns:
    antib = map_to_names[col.split('\n')[0]]
    genestr = ''
    if antib in antibs:
        genestr = ', '.join([f'{", ".join(genes[a])} ({a})' for a in antibs[antib]])
    else:
        if antib in genes:
            genestr = ', '.join(genes[antib])
        else:
            genestr = '-'
            
    comp.loc['Genes conferring resistances', col] = genestr

comp = comp.reindex(['Genes conferring resistances', 'Prediction correctness', 'Correct', 'False'] + rownames)


# write out
with open('comparison_madrid.html', 'w') as outfh:
    outfh.write(htmlheader)
    outfh.write('<h1 class="header" id="main-header">Comparison to Madrid results</h1>\n')

    # # general
    # outfh.write('<h2 class="header" id="params-header">Overall summary</h2>\n')
    # outfh.write('Prediction')

    # results table
    # outfh.write('<h2 class="header" id="table-header">Abricate + CARD</h2>\n')
    outfh.write(comp.to_html(classes=['tablestyle'], escape=False, bold_rows=False, \
            na_rep=f'<font color="{color_error_red}">n/a</font>'))

    outfh.write(htmlfooter)
log(f'Wrote report to comparison_madrid.html')




