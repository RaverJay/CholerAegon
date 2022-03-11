#!/usr/bin/env python3
# SK

import os
import sys

with open('run_readstats.sh', 'w') as outfh:

    outfh.write('#!/bin/bash\n\nset -euo pipefail\n\n')

    with open(sys.argv[1]) as infh:

        for line in infh:
            iso, lr, sr1, sr2 = line.strip().split(',')

            sys.stderr.write(f'{iso} longreads\n')
            outfh.write(f'echo {iso} longreads\n')
            outfh.write(f'cat {lr} > {iso}_lr.fq\n')
            outfh.write(f'read_stats.py {iso}_lr.fq >longreads/{iso}_lr_stats.txt\n')
            outfh.write(f'rm {iso}_lr.fq\n\n')

            sys.stderr.write(f'{iso} shortreads\n')
            outfh.write(f'echo {iso} shortreads\n')
            outfh.write(f'cat {sr1} {sr2} > {iso}_sr.fq.gz\n')
            outfh.write(f'read_stats.py {iso}_sr.fq.gz >shortreads/{iso}_sr_stats.txt\n')
            outfh.write(f'rm {iso}_sr.fq.gz\n\n\n')

sys.stderr.write('Done.\n')
