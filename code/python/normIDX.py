#!/usr/bin/env python
# USAGE: normIDX.py input_file output_file 100
# Normalizes read abundances to gene and read length

import sys
import math

with open(sys.argv[1], 'r') as idxStats:

        outFile = str(sys.argv[1]).rstrip('cstvx') + 'norm.tsv'
        outFile = open(outFile, 'w')
        #outFile = open(sys.argv[2], 'w')
        #outFile.write('target_name\ttarget_length\tnormalized_abundance\n')
        outFile.write('target_name\tnormalized_abundance\n')

        readLen = float(sys.argv[2])
        #readLen = float(sys.argv[3])

        for line in idxStats:
                if line[0] == '*': continue
                else:
                        line = line.split()
                        seqName = str(line[0])
                        seqLen = float(line[1])
                        readAbund = float(line[2])
                        normAbund = (readAbund * readLen) / seqLen

                        #seqLen = str(int(seqLen))
                        normAbund = str(int(math.ceil(normAbund)))

                        #entry = seqName + '\t' + seqLen + '\t' + normAbund + '\n'
                        entry = seqName + '\t' + normAbund + '\n'
                        outFile.write(entry)

outFile.close()
