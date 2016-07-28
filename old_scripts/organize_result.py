#! /usr/bin/env python

import sys

input_file = sys.argv[1]
sample_list_file = sys.argv[2]
mut_info_file = sys.argv[3]
SJ_info_file = sys.argv[4]

id2sample = {}
temp_id = "1"
with open(sample_list_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        id2sample[temp_id] = F[0]
        temp_id = str(int(temp_id) + 1)

mut_id2mut_info = {}
with open(mut_info_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        mut_id2mut_info[F[0] + '\t' + F[1]] = F[2] + '\t' + F[3]

SJ_id2SJ_info = {}
with open(SJ_info_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        SJ_id2SJ_info[F[0] + '\t' + F[1]] = F[2] + '\t' + F[3] + '\t' + F[4]


with open(input_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')

        gene = F[0]
        mutation_states = F[1].split(';')
        link_vector = F[3].split(';')
        BIC_min = F[4]
        BIC0 = F[5]
        conf_min_vector = F[6].split(',')

        mut_id2sample_id = {}
        for mut_state in mutation_states:
            mut_id, sample_ids = mut_state.split(':')
            mut_id2sample_id[mut_id] = sample_ids


        active_link_vector = [link_vector[i] for i in range(len(link_vector)) if conf_min_vector[i] == "1"]

        for active_link in active_link_vector:
            mut_id, SJ_id = active_link.split(',')
            
            # get sample names
            sample_names = []
            for sample_id in mut_id2sample_id[mut_id].split(','):            
                sample_names.append(id2sample[sample_id])
    
            # get mutation info
            mut_info = mut_id2mut_info[gene + '\t' + mut_id]
        
            # get SJ info
            SJ_info = SJ_id2SJ_info[gene + '\t' + SJ_id]

            print gene + '\t' + ';'.join(sample_names) + '\t' + mut_info + '\t' + SJ_info
 
