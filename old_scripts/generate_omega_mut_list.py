#! /usr/bin/env python

import sys, os, glob

mut_dir = sys.argv[1]
sj_dir = sys.argv[2]
tag = sys.argv[3]

mut_files = glob.glob(mut_dir + "/*/*.genomon_mutation.result.filt.txt")
sj_files = glob.glob(sj_dir + "/*/*.SJ.out.tab")

sample2mut_file = {}
for mut_file in sorted(mut_files):
    sample = os.path.basename(os.path.dirname(mut_file))
    sample = sample[0:15]
    if not sample.endswith(tag): continue
    sample2mut_file[sample] = mut_file

sample2sj_file = {}
for sj_file in sorted(sj_files):
    sample = os.path.basename(os.path.dirname(sj_file))
    sample = sample[0:15]
    if not sample.endswith(tag): continue
    sample2sj_file[sample] = sj_file

for sample in sorted(sample2mut_file):
    if sample not in sample2sj_file: continue
    print sample + '\t' + sample2mut_file[sample] + '\t' + sample2sj_file[sample]



