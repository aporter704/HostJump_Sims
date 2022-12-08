#!/usr/bin/env python3
import argparse, re

parser = argparse.ArgumentParser()
parser.add_argument('-s', '--s', help = 'alignment file in fasta format', required = True)

args = parser.parse_args()
aln_file_name = args.s

date_file_name = re.sub('[.]fasta', '.date', aln_file_name)
aln = open(aln_file_name, 'r').readlines()

out = open(date_file_name, 'w')
tip_names_dates = []
for i in aln:
    if '>' in i:
        tip_name = re.sub('\'', '_', re.sub('>|\n', '', i))
        date = re.sub('\'', '', re.sub('.+_|_|\n', '', i))
        tip_names_dates.append(tip_name+'\t'+date)

out.write(str(len(tip_names_dates))+'\n')
for i in tip_names_dates:
    out.write(i+'\n')
out.close()
