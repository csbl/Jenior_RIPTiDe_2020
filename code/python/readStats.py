#!/usr/bin/python
'''USAGE: fastqStats.py fastq
Reports stats for fastq files
'''
import sys
import numpy

q_score = {'!': 0, '\"': 1, '#': 2, '$': 3, '%': 4, '&': 5, '\'': 6, '(': 7, ')': 8, '*': 9, 
'+': 10, ',': 11, '-': 12, '.': 13, '/': 14, '0': 15, '1': 16, '2': 17, '3': 18, 
'4': 19, '5': 20, '6': 21, '7': 22, '8': 23, '9': 24, ':': 25, ';': 26, '<': 27, 
'=': 28, '>': 29, '?': 30, '@': 31, 'A': 32, 'B': 33, 'C': 34, 'D': 35, 'E': 36, 
'F': 37, 'G': 38, 'H': 39, 'I': 40, 'J': 41, 'K': 42}


def read_fastq(fastq):

        lenList = []
        qualList = 0

        currLine = 1
        for line in fastq:
                currLine += 1
                if currLine == 5:
                        quals = list(line.strip())
                        currLen = len(quals)
                        lenList.append(currLen)

                        qual_ave = 0
                        for score in quals:
                                qual_ave += q_score[score]
                        qual_ave = float(qual_ave) / float(currLen)
                        qualList.append(round(qual_ave, 3))

                        currLine = 1

        return lenList,  qualList


def calc_stats(lengths):

        lengths.sort()

        total_reads = len(lengths) # Total number of sequences
        total_Mb = sum(lengths)/1000000.00 # Total number of residues (Mb)
        shortest_read = lengths[0]
        longest_read = lengths[-1]

        if total_reads >= 4:
                median_len = lengths[int(round(total_reads/2))]
                mode_len = max(set(lengths), key=lengths.count)
                mode_count = lengths.count(mode_len)
                
        else:
                mode_len = 'Too few sequences to calculate'
                mode_count = 'Too few sequences to calculate'
                median_len = 'Too few sequences to calculate'

        return([total_reads, total_Mb, shortest_read, longest_read, median_len, mode_len, mode_count])


read_lens, read_quals = read_fastq(open(sys.argv[1], 'r'))

stat_lst = calc_stats(read_lens)
ave_qual = sum(read_quals) / len(read_quals)

output_str = """
# Input file name:\t{filename}
# Total reads:\t\t{reads}
# Average Q-score:\t\t{qual}
# Total bases (Mb):\t{mb} 
# Shortest read length:\t{short}
# Longest read length:\t{long}
# Median read length:\t{med}
# Mode read length:\t{mode}
# Mode frequency:\t{mode_freq}
""".format(filename = str(sys.argv[1]).split('/')[-1],
        reads = stat_lst[0],
        qual = ave_qual,
        mb = "%.2f" % stat_lst[1],
        short = stat_lst[2],
        long = stat_lst[3],
        med = stat_lst[4],
        mode = stat_lst[5],
        mode_freq = stat_lst[6])

print output_str

