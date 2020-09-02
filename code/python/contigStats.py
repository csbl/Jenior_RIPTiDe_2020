#!/usr/bin/python2.7
'''USAGE: contigStats.py seqFile
This script calculates various statistics about the provided fasta file.
'''
import sys
import math

# This function reads in fasta file, appends the length of each sequence to a list, and counts all Gs & Cs.
        # It returns a sorted list of sequence lengths with the G+C % as the last element.
def read_lengths(fasta):

        len_lst = []
        gc_cnt = 0
        n_cnt = 0
        all_cnt = 0
        seq = ''
        firstLine = 'placeholder'
        while firstLine[0] != '>' or firstLine == '\n':
                firstLine = fasta.readline()

        for line in fasta:
                if line[0] == '>': 
                        gc_cnt += seq.count('G')
                        gc_cnt += seq.count('C')
                        n_cnt += seq.count('N')
                        all_cnt += len(seq)
                        len_lst.append(len(seq))
                        seq = ''
                else:
                        seq += line.strip().upper()

        gc_cnt += seq.count('G')
        gc_cnt += seq.count('C')
        n_cnt += seq.count('N')
        all_cnt += len(seq)
        len_lst.append(len(seq))

        gc_prcnt = round((float(gc_cnt) / float(all_cnt)) * 100.0, 2)
        len_lst.sort()

        return(len_lst, gc_prcnt, n_cnt)


# Function t0 calculate satdard deviation for a list of numbers
def standDev(values):
        x_mean = sum(values) / len(values)
        sd_list = []
        for x in values:
                y = (x - x_mean) ** 2
                sd_list.append(y)
        y_mean = sum(sd_list) / len(sd_list)   
        sd = math.sqrt(y_mean)
        return(sd)


# This function calculates and returns all the printed statistics.
def calc_stats(lengths):

        shortest = lengths[0]
        longest = lengths[-1]
        total_contigs = len(lengths) # Total number of sequences
        len_sum = sum(lengths) # Total number of residues
        total_Mb = len_sum/1000000.00 # Total number of residues expressed in Megabases
        mid_pos = int(round(total_contigs/2))

        median_len = lengths[mid_pos] # Median sequence length
        q1 = lengths[0:mid_pos][int(len(lengths[0:mid_pos])/2)]
        q3 = lengths[mid_pos:-1][int(len(lengths[mid_pos:-1])/2)]
        iqr = q3 - q1

        mean_len = sum(lengths) / len(lengths)
        sd = standDev(lengths)

        # Pearsons Coefficient of Skewness
        skew = ((mean_len - median_len) * 3) / sd
        mean_len = round(mean_len, 2)
        sd = round(sd, 2)
        skew = round(skew, 3)
 
        current_bases = 0
        n50 = 0
        n90 = 0
        seqs_1000 = 0
        seqs_5000 = 0
        seqs_10000 = 0
        percent50_bases = int(round(len_sum*0.5))
        percent90_bases = int(round(len_sum*0.1))

        for x in lengths:

                current_bases += x

                if x > 1000:
                        seqs_1000 += 1
                if x > 5000:
                        seqs_5000 += 1
                if x > 10000:
                        seqs_10000 += 1

                if current_bases >= percent50_bases and n50 == 0:
                        n50 = x
                if current_bases >= percent90_bases and n90 == 0:
                        n90 = x

        l50 = lengths.count(n50)
        short_contigs = total_contigs - seqs_1000

        return(total_contigs, total_Mb, n50, l50, n90, median_len, q1, q3, iqr, seqs_1000, seqs_5000, seqs_10000, shortest, longest, mean_len, sd, skew, short_contigs)

#----------------------------------------------------------------------------------------#

with open(sys.argv[1], 'r') as contigs:
        contig_lengths, gc_content, n_count = read_lengths(open(sys.argv[1], 'r'))

if len(contig_lengths) < 5:
        print('\nToo few contigs in ' + str(sys.argv[1]) + ' to calculate useful statistics. Exiting.\n')
        sys.exit()

stat_list = calc_stats(contig_lengths)
stat_list = list(stat_list)
stat_list.append(gc_content)
stat_list.append(n_count)

output_string = """
Assembly:       {filename}
Contigs:        {total_contigs}
Bases(Mb):      {total_mb}
Ns:             {ns}
G+C(%):         {gc}
-----------------------------------
N50:            {n50}
L50:            {l50}
N90:            {n90}
-----------------------------------
Median:         {median_len}
Q1:             {q1}
Q3:             {q3}
IQR:            {iqr}
Mean:           {mean}
Std:            {sd}
Skewness:       {skewness}
-----------------------------------
Shortest:       {short}
Longest:        {long}
Contigs<1kb:    {short_contigs}
Contigs>1kb:    {seqs_1k}
Contigs>5kb:    {seqs_5k}
Contigs>10kb:   {seqs_10k}

###################################""".format(filename = str(sys.argv[1]).split('/')[-1],  
        total_contigs = stat_list[0], 
        total_mb = "%.2f" % stat_list[1], 
        n50 = stat_list[2],
        l50 = stat_list[3],
        n90 = stat_list[4],
        median_len = stat_list[5], 
        q1 = stat_list[6], 
        q3 = stat_list[7],
        iqr = stat_list[8],
        mean = stat_list[14],
        sd = stat_list[15],
        skewness = stat_list[16],
        short = stat_list[12], 
        long = stat_list[13],
        short_contigs = stat_list[17], 
        seqs_1k = stat_list[9], 
        seqs_5k = stat_list[10],
        seqs_10k = stat_list[11],
        gc = stat_list[18],
        ns = stat_list[19])

print(output_string)
