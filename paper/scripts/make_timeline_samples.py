#!/usr/bin/env python3
# SK

import os
import sys
import time
import numpy as np
import pandas as pd


def error(string, error_type=1):
    sys.stderr.write(f'ERROR: {string}\nAborting.\n')
    sys.exit(error_type)


def log(string, newline_before=False):
    if newline_before:
        sys.stderr.write('\n')
    sys.stderr.write(f'LOG: {string}\n')


fastq_in = sys.argv[1]

sample = os.path.basename(fastq_in).split('_')[0]
log(sample)

# get readids and timestamps
readids = []
timestamps = {}

log('Parse reads ...')

with open(fastq_in) as infh:

    ln = 0
    for line in infh:

        # only use header lines
        if ln==0:
            assert line.startswith('@')
            lt = line.strip().split()
            readid = lt[0][1:]
            timestr = lt[5]
            assert timestr.startswith('start_time=')
            timestr = timestr.split('=')[1]
            # e.g. 2020-05-06T16:18:29Z
            timestamp = time.mktime(time.strptime(timestr, '%Y-%m-%dT%H:%M:%SZ'))
            
            readids.append(readid)
            timestamps[readid] = timestamp
        
        ln = (ln + 1) % 4

log(f'Total reads: {len(readids)}')
log('Sort readids ...')

readids.sort(key=lambda x: timestamps[x])

for i in range(5):
    log(readids[i], timestamps[readids[i]])

global_starttime = timestamps[readids[0]]
latest_read = timestamps[readids[-1]] - global_starttime

log(f'Latest read: {latest_read}')
log(f'Max minutes: {latest_read//60}')

outdir = f'samples/{sample}'
os.makedirs(outdir, exist_ok=True)


limits_minutes = list(range(1, 20)) + list(range(20, 50, 5)) + list(range(50, 200, 10)) + list(range(200, 1000, 50)) + list(range(1000, int(latest_read)//60, 100))
log(limits_minutes)

for maxtime in limits_minutes:

    outfile = f'{outdir}/{sample}_{maxtime:04d}min.ids'
    with open(outfile, mode='w') as outfh:

        for readid in readids:
            if (timestamps[readid] - global_starttime) <= maxtime*60:
                outfh.write(readid + '\n')
            else:
                break
    log(f'Wrote {outfile}')

log('Done.')
