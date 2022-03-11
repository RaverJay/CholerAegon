#!/usr/bin/env python3
'''
Generate a summary report for multiple samples
./summary_report.py
'''
# SK

import sys
import time
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

class SummaryReport():

    report_time = time.localtime()
    report_name = f'CholerAegon_summary_report_{time.strftime("%Y-%m-%d--%H-%M-%S", report_time)}'
    output_filename = report_name + '.html'
    # CholerAegon_params = {}
    # tool_versions = {}

    tabledata = None
    tabledataraw = None
    col_formatters = {}
    col_descriptions = []
    
    sample_QC_status = None
    sample_QC_info = {}
    control_string_patterns = ['control', 'negative']

    # colors
    color_markup = '#8a006d'
    color_good_green = '#046907'
    color_warn_orange = '#ac7800'
    color_error_red = '#a50500'


    def init(self, report_name):
        if report_name is not None:
            self.report_name = report_name
            self.output_filename = report_name + '.html'
        log(f'Created report object: {self.report_name}')


    def add_column(self, column_name, pandas_series):
        self.tabledata[column_name] = pandas_series

    def add_column_raw(self, column_name, pandas_series):
        self.tabledataraw[column_name] = pandas_series


    def add_col_formatter(self, colname, colformatter):
        assert colname not in self.col_formatters, f'Duplicate column formatter: {colname}'
        self.col_formatters[colname] = colformatter


    def add_param(self, param_name, param_value):
        assert param_name not in self.porecov_params, f'Duplicate parameter: {param_name}'
        self.porecov_params[param_name] = param_value
        log(f'Added porecov param: {param_name}: {param_value}')


    def add_QC_info(self, info_name, info_value):
        assert info_name not in self.sample_QC_info, f'Duplicate QC info: {info_name}'
        self.sample_QC_info[info_name] = info_value
        log(f'Added QC info: {info_name}: {info_value}')


    def add_time_param(self):
        self.add_param('Report created', f'{time.strftime("%Y-%m-%d %H:%M:%S %Z", self.report_time)}')


    def add_version_param(self, porecov_version):
        pc_param = '<a href="https://github.com/RaverJay/CholerAegon"><b>poreCov</b></a> version'
        warning_msg = 'Warning: Not an official release version of CholerAegon. Use parameter \'-r\' to specify a release version.'
        revision, commitID, scriptID = porecov_version.split(':')
        if revision != 'null':
            self.add_param(pc_param, revision)
        else:
            if commitID != 'null':
                self.add_param(pc_param, commitID + ' (git commitID) - ' + warning_msg)
            else:
                self.add_param(pc_param, scriptID + ' (local scriptID) - ' + warning_msg)


    # def parse_version_config(self, version_config_file):
    #     version_dict = {}
    #     try:
    #         with open(version_config_file) as infh:
    #             for line in infh:
    #                 lt = line.strip().replace(' ','')
    #                 if lt.startswith('withLabel:'):
    #                     tname, tinfo = lt.split(':',1)[1].split('{',1)
    #                     tcontainer = tinfo.split("'")[1] if "'" in tinfo else tinfo.split('"')[1]
    #                     tversion = tcontainer.split(':')[1].split('-')[0].lstrip('v')

    #                     assert tname not in version_dict, f'Duplicate tool in version config: {tname}'
    #                     version_dict[tname] = tversion
    #     except:
    #         print(version_dict)
    #         error(f'Failed to parse version config file: {version_config_file}')
    #     self.tool_versions = version_dict
    #     log('Parsed version config file.')



    def check_and_init_table_with_samples(self, samples):
        if samples != 'samples_list.csv':
            log('No sample list input.')
        else:
            log('Using samples input.')
            s_list = [s.strip() for s in open(samples).readlines()]
            log(f'Samples: {s_list}')

            s_table = pd.DataFrame(index=s_list)
            self.force_index_dtype_string(s_table)
            self.check_and_init_tabledata(s_table.index)


    def check_and_init_tabledata(self, t_index):
        '''If tabledata is None, initialize it now. Then check if all new index values are in the existing table.
        Thus samples input or adding the results with the most samples (kraken2) first is required.'''
        if self.tabledata is None:
            self.tabledata = pd.DataFrame(index=sorted(t_index))
            self.tabledata.columns.name = 'Sample'
            self.force_index_dtype_string(self.tabledata)

            # TODO
            # self.add_col_description(f'Missing values (\'<font color="{self.color_error_red}">.</font>\') denote cases where abricate did not find the gene.')
            # self.tabledataraw = self.tabledata.copy()
        else:
            for item in t_index:
                assert item in self.tabledata.index, f'Index not found in existing table: {item}. Available: {self.tabledata.index}'


    def add_col_description(self, desc):
        self.col_descriptions.append(f'{desc}')


    def reindex_tabledata(self):
        log('Setting a number range as table index ...')
        self.tabledata.index = pd.Series(range(1, self.tabledata.shape[0]+1))
        self.tabledata.columns.name = 'Assembly number'

    
    def reorder_tablecolumns(self):
        log('Reordering table columns ...')
        prioritized = ['Sample', 'Assembly method', 'FastANI %', '# Genes found']
        frontcols = [col for col in prioritized if col in self.tabledata.columns]
        othercols = [col for col in self.tabledata.columns if col not in frontcols]
        self.tabledata = self.tabledata[frontcols + othercols]


    ### html writing functions


    def write_column_descriptions(self, filehandle):
        for desc in self.col_descriptions:
            filehandle.write(f'{desc}<br>\n')


    def write_html_table(self, filehandle):
        filehandle.write(self.tabledata.to_html(classes=['tablestyle'], escape=False, bold_rows=False, \
            na_rep=f'<font color="{self.color_error_red}">n/a</font>', formatters=self.col_formatters, float_format=lambda f: f'{f:.2f}'))


    def write_html_report(self):
        '''Write the html report to a file'''

        htmlheader = '''<!DOCTYPE html><html><head>
        <title>''' + self.report_name + '''</title>
        <meta charset="utf-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">

        <style>
        * {
            font-family:"Open Sans",freesans,sans-serif
        }

        .content {
        max-width: 1700px;
        margin: auto;
        }

        table.tablestyle {
        background-color: #FFFFFF;
        text-align: center;
        border-collapse: collapse;
        }
        table.tablestyle td, table.tablestyle th {
        border: 2px solid #8B8B8B;
        padding: 5px 5px;
        }
        table.tablestyle tbody td {
        font-size: 16px;
        color: #000000;
        }
        table.tablestyle tr:nth-child(even) {
        background: #E6F5FF;
        }
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

        with open(self.output_filename, 'w') as outfh:
            outfh.write(htmlheader)
            outfh.write('<h1 class="header" id="main-header">CholerAegon Summary Report</h1>\n')

            # # general
            # outfh.write('<h2 class="header" id="params-header">Run information</h2>\n')
            # for info, value in self.sample_QC_info.items():
            #     outfh.write(value + '<br>\n')
            # outfh.write('<br>\n')
            # for param, value in self.porecov_params.items():
            #     outfh.write(param + ': ' + value + '<br>\n')

            # results table
            outfh.write('<h2 class="header" id="table-header">Sample results</h2>\n')
            self.write_html_table(outfh)
            self.write_column_descriptions(outfh)

            outfh.write(htmlfooter)
        log(f'Wrote report to {self.output_filename}.')


    ### functions to add columns

    def force_index_dtype_string(self, dataframe):
        dataframe.index = dataframe.index.astype('string')


    def add_amr_results(self, amr_results):
        log('Adding RGI results ...')

        res_data = pd.read_csv(amr_results, index_col='#FILE', sep='\t', dtype={'#FILE': str})
        self.force_index_dtype_string(res_data)

        # parsing
        data = res_data.copy().drop('NUM_FOUND', axis=1)

        assemblies = []
        assembly_types = []


        for row in data.index:

            assembly = row

            # try to infer assembly type
            sample, asmtype = row.rsplit('_', 1)
            if asmtype in ['longreads', 'shortreads', 'hybrid']:

                # pipeline-assembled
                if asmtype not in assembly_types:
                    assembly_types.append(asmtype)

                data.loc[row, 'SAMPLE'] = sample
                data.loc[row, 'ASSEMBLY_TYPE'] = asmtype

            else:
                # pre-assembled fasta input
                data.loc[row, 'SAMPLE'] = assembly
                data.loc[row, 'ASSEMBLY_TYPE'] = 'pre-assembled'
                
            if assembly not in assemblies:
                assemblies.append(assembly)


        def stringify(strlst):
            return ', '.join(strlst)

        log(f'Found assembly types: {stringify(assembly_types)}')
        log(f'List of assemblies: {stringify(assemblies)}')

        # sort
        data['NUM_FOUND'] = res_data['NUM_FOUND']
        data.sort_values(by=['SAMPLE', 'ASSEMBLY_TYPE'], inplace=True)
        

        # select columns
        cols_other = ['SAMPLE', 'ASSEMBLY_TYPE', 'NUM_FOUND']
        cols_amr = [col for col in data.columns if col not in cols_other]

        data = data[cols_other + cols_amr]

        # set tabledata
        self.check_and_init_tabledata(data.index)

        renaming = {'SAMPLE': 'Sample', 'ASSEMBLY_TYPE': 'Assembly method', 'NUM_FOUND': '# Genes found'}
        for col in data.columns:
            self.add_column(col if col not in renaming else renaming[col], data[col])

        self.add_col_description('Resistance genes were determined with <a href="https://card.mcmaster.ca/analyze/rgi">RGI (CARD resistance gene identifier)</a> ' + \
            'and <a href="https://github.com/tseemann/abricate">Abricate</a> (with CARD database). When both identified the same gene, the higher %coverage hits were chosen.')



    # unused
    def add_abricate_results(self, abricate_results):
        log('Adding abricate results ...')

        res_data = pd.read_csv(abricate_results, index_col='#FILE', sep='\t', dtype={'#FILE': str})
        self.force_index_dtype_string(res_data)

        # parsing
        data = res_data.copy().drop('NUM_FOUND', axis=1)

        databases = []
        isolates = []
        assembly_types = []

        for ab_res in data.index:
            fname, ext = ab_res.split('.')
            if ext != 'csv':
                log(f'Encountered abnormal abricate result file extension: {ab_res}')
            
            abr, database, isolate, asmtype = fname.split('_', 3)
            if abr != 'abricate':
                log(f'Encountered non-abricate result filename: {abr}')

            if database not in databases:
                databases.append(database)
            if isolate not in isolates:
                isolates.append(isolate)
            if asmtype not in assembly_types:
                assembly_types.append(asmtype)

            data.loc[ab_res, 'DATABASE'] = database
            data.loc[ab_res, 'ISOLATE'] = isolate
            data.loc[ab_res, 'ASSEMBLY_TYPE'] = asmtype

        def stringify(strlst):
            return ', '.join(strlst)

        log(f'Found databases: {stringify(databases)}')
        log(f'Found assembly types: {stringify(assembly_types)}')
        log(f'List of isolates: {stringify(isolates)}')

        # sort
        data['NUM_FOUND'] = res_data['NUM_FOUND']
        data.sort_values(by=['ISOLATE', 'ASSEMBLY_TYPE', 'DATABASE'], inplace=True)

        cols_other = ['ISOLATE', 'ASSEMBLY_TYPE', 'DATABASE', 'NUM_FOUND']
        cols_amr = [col for col in data.columns if col not in cols_other]

        data = data[cols_other + cols_amr]
        
        # TODO REFACTOR THIS - consistent assembly names in the pipeline
        self.check_and_init_tabledata(data.index)

        renaming = {'ISOLATE': 'Isolate', 'ASSEMBLY_TYPE': 'Assembly method', 'DATABASE': 'Database', 'NUM_FOUND': '# Genes found'}

        for col in data.columns:
            self.add_column(col if col not in renaming else renaming[col], data[col])

        self.add_col_description('Resistance genes were determined with abricate.')


    def add_fastani_results(self, fastani_results):
        log('Adding FastANI results ...')

        res_data = pd.read_csv(fastani_results, index_col=0, sep='\t', header=None)
        self.force_index_dtype_string(res_data)

        res_data.columns = ['reference', 'ani_value', 'frag_mapped', 'frag_total']

        # rename index
        res_data.index = [name.rsplit('.',1)[0].rsplit('_assembly',1)[0] for name in res_data.index]

        self.add_column('FastANI %', res_data['ani_value'])

        self.add_col_description('ANI% was determined with <a href="https://github.com/ParBLiSS/FastANI">FastANI</a> against "Vibrio cholerae O1 biovar El Tor str. N16961" (<a href="https://www.ncbi.nlm.nih.gov/genome/?term=NC_002505.1">NC_002505.1</a>).')



    def check_if_control(self, sample_name):
        for pattern in self.control_string_patterns:
            if pattern in sample_name:
                return True
        return False


    def add_QC_status_info(self):
        if self.sample_QC_status is None:
            error('sample_QC_status was not set before calling add_QC_status_info().')

        n_realsamples = 0
        n_controls = 0
        n_passrealsamples = 0
        n_passcontrols = 0
        for sample, status in self.sample_QC_status.items():
            if self.check_if_control(sample):
                n_controls += 1
                if status == 'pass':
                    n_passcontrols += 1
            else:
                n_realsamples += 1
                if status == 'pass':
                    n_passrealsamples += 1

        # add status info
        if n_passrealsamples > 0:
            self.add_QC_info('Passed samples', f'<font color="{self.color_good_green}"><b>{n_passrealsamples} / {n_realsamples} of samples passed QC criteria.</b></font>')
        if n_passrealsamples < n_realsamples:
            self.add_QC_info('Failed samples', f'<font color="{self.color_error_red}"><b>{n_realsamples-n_passrealsamples} / {n_realsamples} of samples failed QC criteria.</b></font>')
        if n_controls > 0:
            if n_passcontrols < n_controls:
                self.add_QC_info('Negative controls', f'<font color="{self.color_good_green}"><b>{n_controls-n_passcontrols} / {n_controls} of control samples correctly did not produce an assembly that passed QC criteria.</b></font>')
            if n_passcontrols > 0:
                self.add_QC_info('Bad controls', f'<font color="{self.color_error_red}"><b>{n_passcontrols} / {n_controls} of control samples wrongly produced an assembly that passed QC criteria.</b></font>')
        else:
            self.add_QC_info('Negative control', f'<font color="{self.color_warn_orange}"><b>No negative control samples were found by automatic detection.</b></font>')
        patterns = "'" + "', '".join(self.control_string_patterns) + "'"
        self.add_QC_info('Note', f'Note: samples are considered negative controls if their name contains certain keywords ({patterns}) - please check if these assignments were correct.')

        # mark control samples
        def mark_controls(sample_name):
            if self.check_if_control(sample_name):
                return sample_name + f'<br>(<font color="{self.color_spike_markup}">considered control</font>)'
            else:
                return sample_name

        self.tabledata.index = [mark_controls(sn) for sn in self.tabledata.index]


    def write_table_output(self):
        if self.tabledataraw is None:
            self.tabledataraw = self.tabledata.copy()
        # self.tabledataraw.to_excel(self.report_name + '_datatable.xlsx' , sheet_name='poreCov', index_label='sample')
        self.tabledataraw.to_csv(self.report_name + '_datatable.tsv', index_label='Assembly number', sep='\t')



###

if __name__ == '__main__':

    log('Started summary_report.py ...')

    parser = argparse.ArgumentParser(description='Generate a summary report for multiple samples run with CholerAegon')
    # parser.add_argument("-v", "--version_config", help="version config", required=True)
    # parser.add_argument("--CholerAegon_version", help="CholerAegon version", required=True)
    parser.add_argument("--amr_results", help="rgi results table", required=True)
    parser.add_argument("--abricate_results", help="abricate results table")
    parser.add_argument("--fastani_results", help="fastani results table")
    args = parser.parse_args()


    ### build report
    report = SummaryReport()


    # # check for samples input 
    # if args.samples:
    #     report.check_and_init_table_with_samples(args.samples)


    # results table
    if args.amr_results:
        report.add_amr_results(args.amr_results)
    if args.fastani_results:
        report.add_fastani_results(args.fastani_results)

    # finalize table
    report.reindex_tabledata()
    report.reorder_tablecolumns()

    # metadata

    # # total QC status
    # report.add_QC_status_info()


    report.write_table_output()
    report.write_html_report()

