#! /usr/bin/env python

import sys

input_file = sys.argv[1]
mut2sample_file = sys.argv[2]
output_count_file = sys.argv[3]
output_mut_file = sys.argv[4]
output_SJ_file = sys.argv[5]

def get_mut_sample_info(mut_info, mut2sample):
    mut_infos = mut_info.split(',')
    tchr, tstart, tend, tref, talt = mut_infos[0], mut_infos[1], mut_infos[1], mut_infos[3], mut_infos[4]
    # deletion
    if len(tref) > 1:
        tref = tref[1:]
        talt = "-"
        tstart = str(int(tstart) + 1)
        tend = str(int(tstart) + len(tref) - 1)
    if len(talt) > 1:
        tref = "-"
        talt = talt[1:]

    mut = '\t'.join([tchr, tstart, tend, tref, talt])
    return(mut2sample[mut])


hout1 = open(output_count_file, 'w')
hout2 = open(output_mut_file, 'w')
hout3 = open(output_SJ_file, 'w')
    
mut2sample = {}
with open(mut2sample_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        mut = '\t'.join(F[0:5])
        mut2sample[mut] = F[5]


header2ind = {}
with open(input_file, 'r') as hin:

    header = hin.readline().rstrip('\n').split('\t')
    for i in range(len(header)):
        header2ind[header[i]] = i

    # print >> hout0, "Gene" + '\t' + '\t'.join(header)
    temp_gene = ""
    temp_mut2info = {} 
    temp_mut2sample = []
    temp_mut2id = {}
    temp_id2mut = {}
    temp_SJ2info = {} 
    temp_SJ_count = []
    temp_SJ2id = {} 
    temp_id2SJ = {}
    temp_mut_SJ_link = []
    temp_mut_id = "1"
    temp_SJ_id = "1"
    for line in hin:
        F = line.rstrip('\n').split('\t')

        if F[header2ind["Gene_Symbol"]] != temp_gene:
            if temp_gene != "":

                # flush the result
                print >> hout1, temp_gene + '\t' + ';'.join(temp_mut2sample) + '\t' + ';'.join(temp_SJ_count) + '\t' + ';'.join(temp_mut_SJ_link)
            
                for id in sorted(temp_id2mut):
                    print >> hout2, temp_gene + '\t' + id + '\t' + temp_id2mut[id] + '\t' + ';'.join(temp_mut2info[temp_id2mut[id]])

                for id in sorted(temp_id2SJ):
                    print >> hout3, temp_gene + '\t' + id + '\t' + temp_id2SJ[id] + '\t' + temp_SJ2info[temp_id2SJ[id]]

            temp_gene = F[header2ind["Gene_Symbol"]] 
            temp_mut2info = {} 
            temp_mut2sample = []
            temp_mut2id = {}
            temp_id2mut = {}
            temp_SJ2info = {} 
            temp_SJ_count = []
            temp_SJ2id = {} 
            temp_id2SJ = {}
            temp_mut_SJ_link = []
            temp_mut_id = "1"
            temp_SJ_id = "1"

        mut = F[header2ind["Mutation_Info"]]
        sample = get_mut_sample_info(F[header2ind["Mutation_Info"]], mut2sample)
        mut_info = F[header2ind["Splicing_Motif_Pos"]] + ',' + F[header2ind["Splicing_Mutation_Type"]] + ',' + F[header2ind["Is_Cannonical"]]

        if mut not in temp_mut2info:
            temp_mut2info[mut] = []
            temp_mut2id[mut] = temp_mut_id
            temp_id2mut[temp_mut_id] = mut
            temp_mut2sample.append(temp_mut_id + ':' + sample)
            temp_mut_id = str(int(temp_mut_id) + 1)

        if mut_info not in temp_mut2info[mut]: temp_mut2info[mut].append(mut_info) 


        SJ = F[header2ind["SJ_1"]] + ':' + F[header2ind["SJ_2"]] + '-' + F[header2ind["SJ_3"]]
        if SJ not in temp_SJ2info:
            temp_SJ2info[SJ] = '\t'.join([F[header2ind[x]] for x in ["Splicing_Class", "Is_Inframe", "Gene_1", "Exon_Num_1", "Is_Boundary_1", "Gene_2", "Exon_Num_2", "Is_Boundary_2"]])
            temp_SJ2id[SJ] = temp_SJ_id
            temp_id2SJ[temp_SJ_id] = SJ
            temp_SJ_count.append(F[header2ind["SJ_4"]])
            temp_SJ_id = str(int(temp_SJ_id) + 1)

        if temp_mut2id[mut] + ',' + temp_SJ2id[SJ] not in temp_mut_SJ_link:
            temp_mut_SJ_link.append(temp_mut2id[mut] + ',' + temp_SJ2id[SJ])


print >> hout1, temp_gene + '\t' + ';'.join(temp_mut2sample) + '\t' + ';'.join(temp_SJ_count) + '\t' + ';'.join(temp_mut_SJ_link)

for id in sorted(temp_id2mut):
    print >> hout2, temp_gene + '\t' + id + '\t' + temp_id2mut[id] + '\t' + ';'.join(temp_mut2info[temp_id2mut[id]]) 

for id in sorted(temp_id2SJ):
    print >> hout3, temp_gene + '\t' + id + '\t' + temp_id2SJ[id] + '\t' + temp_SJ2info[temp_id2SJ[id]]

hout1.close()
hout2.close()
hout3.close()

