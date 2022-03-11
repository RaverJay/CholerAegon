#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import gzip
import numpy as np
import pandas as pd

def error(string, error_type=1):
    sys.stderr.write('ERROR: ' + string + '\n')
    exit(error_type)

def log(string):
    sys.stderr.write('LOG: ' + string + '\n')

#####

data = pd.DataFrame(index=['File', 'longest', 'second_longest', 'completeness'])

for file in sys.argv[1:]:

    read_lengths = []
    gzipped = False

    if file.endswith('.gz'):
        gzipped = True

    if file.rsplit('.gz',1)[0].endswith('.sam'):

        with gzip.open(file, 'rt') if gzipped else open(file) as fh:

            for line in fh:
                if line.startswith('@'):
                    continue
                lt = line.strip().split('\t')
                read_name = lt[0]
                flags = lt[1]
                # check if primary alignment (0 or 16 (rev)) or unmapped (4)
                if flags not in ['0', '4', '16']:
                    continue
                
                read_lengths.append(int(len(lt[9])))

    elif file.rsplit('.gz',1)[0].endswith(('.fasta', '.fa')):

        with gzip.open(file, 'rt') if gzipped else open(file) as fh:

            first = True
            seq = ''
            for line in fh:

                if line.startswith('>'):
                    # new seq
                    if first:
                        first = False
                    else:
                        read_lengths.append(len(seq))
                        seq = ''
                    
                else:
                    # continue seq
                    seq += line.strip()
                    
            # end, last seq
            read_lengths.append(len(seq))

    elif file.rsplit('.gz',1)[0].endswith(('.fastq', '.fq')):
        
        with gzip.open(file, 'rt') if gzipped else open(file) as fh:

            for line in fh:

                l4 = [line.strip(), fh.readline().strip(), fh.readline().strip(), fh.readline().strip()]
                read_name = l4[0].split()[0][1:]
                # if read_name in read_names:
                #     error(f'Duplicate read names in {file}')
                
                read_lengths.append(int(len(l4[1])))

    if read_lengths == []:
        error(f'ERROR: file contained no reads: {file}')

    read_lengths = sorted(read_lengths, reverse=True)

    reflen = 4_144_502
    maxlen = read_lengths[0]
    max2len = read_lengths[1] if len(read_lengths) > 1 else 0
    completeness = (maxlen + max2len) / reflen

    data[file] = [file, maxlen, max2len, completeness]



data.T.to_csv(sys.stdout, sep='\t', index=False)
