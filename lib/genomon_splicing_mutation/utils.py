#! /usr/bin/env python

import sys, glob, subprocess, re, pysam

def merge_SJ(SJ_list_file, output_file, control_file, junc_num_thres):

    if control_file is not None:
        control_db = pysam.TabixFile(control_file)

    file_list = [] 
    with open(SJ_list_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            file_list.append(F[2])


    junc2list = {}
    for file in file_list:
        print >> sys.stderr, "proceccing: " + file
        with open(file, 'r') as hin:
            for line in hin:

                F = line.rstrip('\n').split('\t')
                if F[5] != "0": continue
                if int(F[6]) < junc_num_thres: continue
                key = F[0] + '\t' + F[1] + '\t' + F[2]
                if key not in junc2list: junc2list[key] = 1


    junc2count = {}
    junc2control = {}
    ind = 0
    for file in file_list:
        print >> sys.stderr, "processing: " + file
        with open(file, 'r') as hin:
            for line in hin:
                F = line.rstrip('\n').split('\t')
                key = F[0] + '\t' + F[1] + '\t' + F[2]

                if key not in junc2list: continue

                if key in junc2control:
                    if junc2control[key] == 1: continue
                else:

                    # remove control files
                    tabixErrorFlag = 0
                    if control_file is not None:
                        try:
                            records = control_db.fetch(F[0], int(F[1]) - 5, int(F[1]) + 5)
                        except Exception as inst:
                            print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
                            tabixErrorMsg = str(inst.args)
                            tabixErrorFlag = 1

                    control_flag = 0;
                    if tabixErrorFlag == 0:
                        for record_line in records:
                            record = record_line.split('\t')
                            if F[0] == record[0] and F[1] == record[1] and F[2] == record[2]:
                                control_flag = 1

                    if control_flag == 1: 
                        junc2control[key] = 1
                        continue
                    else: junc2control[key] = 0

                if key not in junc2count: junc2count[key] = ["0"] * len(file_list)
                junc2count[key][ind] = F[6]

        ind = ind + 1


    hout = open(output_file, 'w') 
    for junc in sorted(junc2count):
        print >> hout, junc + '\t' + ','.join(junc2count[junc])

    hout.close()



def merge_mut(sample_list_file, output_file):


    mut2sample = {}
    sample_num = "1"
    with open(sample_list_file, 'r') as hin1:
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

    hout = open(output_file, 'w')
    for mut in sorted(mut2sample):
        print >> hout, mut + '\t' + ','.join(mut2sample[mut]) 

    hout.close()



def add_gene_symbol(input_file, output_file):

    header2ind = {}
    header = ""
    hout0 = open(output_file + ".tmp0", 'w')
    hout1 = open(output_file + ".tmp1", 'w')
    with open(input_file, 'r') as hin:

        header = hin.readline().rstrip('\n').split('\t')
        for i in range(len(header)):
            header2ind[header[i]] = i

        print >> hout0, "Gene_Symbol" + '\t' + '\t'.join(header)

        for line in hin:
            F = line.rstrip('\n').split('\t')
            genes = F[header2ind["Gene_1"]].split(';') + F[header2ind["Gene_2"]].split(';')
            genes = list(set(genes))

            if "---" in genes: genes.remove("---")
            if len(genes) > 0: 
                genes_nm = filter(lambda x: x.find("(NM_") > 0, genes)
                if len(genes_nm) > 0: genes = genes_nm

            gene = genes[0]
            gene = re.sub(r"\(N[MR]_\d+\)", "", gene)

            print >> hout1, gene + '\t' + '\t'.join(F)

    hout0.close()
    hout1.close()

    hout2 = open(output_file + ".tmp2", 'w')
    subprocess.call(["sort", "-k1", output_file + ".tmp1"], stdout = hout2)
    hout2.close()

    hout3 = open(output_file, 'w')
    subprocess.call(["cat", output_file + ".tmp0", output_file + ".tmp2"], stdout = hout3)
    hout3.close()

    subprocess.call(["rm", "-rf", output_file + ".tmp0"])
    subprocess.call(["rm", "-rf", output_file + ".tmp1"])
    subprocess.call(["rm", "-rf", output_file + ".tmp2"])



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


def organize_mut_SJ_count(input_file, mut2sample_file, output_count_file, output_mut_file, output_SJ_file):

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


    # last flush
    print >> hout1, temp_gene + '\t' + ';'.join(temp_mut2sample) + '\t' + ';'.join(temp_SJ_count) + '\t' + ';'.join(temp_mut_SJ_link)

    for id in sorted(temp_id2mut):
        print >> hout2, temp_gene + '\t' + id + '\t' + temp_id2mut[id] + '\t' + ';'.join(temp_mut2info[temp_id2mut[id]]) 

    for id in sorted(temp_id2SJ):
        print >> hout3, temp_gene + '\t' + id + '\t' + temp_id2SJ[id] + '\t' + temp_SJ2info[temp_id2SJ[id]]

    hout1.close()
    hout2.close()
    hout3.close()

