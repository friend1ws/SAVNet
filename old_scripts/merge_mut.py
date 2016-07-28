#! /usr/bin/env python

import sys

mut_file_list = sys.argv[1]
mut2sample_file = sys.argv[2]

hout = open(mut2sample_file, 'w')

mut2sample = {}
sample_num = "1"
with open(mut_file_list, 'r') as hin1:
    for line1 in hin1:
        F1 = line1.rstrip('\n').split('\t')
        sample = F1[0]
        mut_file = F1[1]
        with open(mut_file, 'r') as hin2:
            for line2 in hin2:
                F2 = line2.rstrip('\n').split('\t')
                if F2[0].startswith('#'): continue
                if F2[0] == "Chr": continue

                key = '\t'.join(F2[0:5])
        
                if key not in mut2sample: 
                    mut2sample[key] = []
               
                mut2sample[key].append(sample_num)

        sample_num = str(int(sample_num) + 1)

for mut in sorted(mut2sample):
    print >> hout, mut + '\t' + ','.join(mut2sample[mut]) 


hout.close()
        
