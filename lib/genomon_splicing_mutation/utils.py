#! /usr/bin/env python

import sys, glob, subprocess, re, math, copy, random, pysam

def merge_SJ2(SJ_file_list, output_file, control_file, junc_num_thres):

    # list up junctions to pick up
    junc2list = {}
    for SJ_file in SJ_file_list:
        with open(SJ_file, 'r') as hin:
            for line in hin:
                F = line.rstrip('\n').split('\t')
                if F[5] != "0": continue
                if int(F[6]) < junc_num_thres: continue
                key = F[0] + '\t' + F[1] + '\t' + F[2]
                if key not in junc2list: junc2list[key] = 1


    temp_id = 0
    hout = open(output_file + ".tmp.unsorted.txt", 'w')
    for SJ_file in SJ_file_list:
        with open(SJ_file, 'r') as hin:
            for line in hin:
                F = line.rstrip('\n').split('\t')
                if F[0] + '\t' + F[1] + '\t' + F[2] in junc2list:
                    print >> hout, F[0] + '\t' + F[1] + '\t' + F[2] + '\t' + str(temp_id) + '\t' + F[6]
      
            temp_id = temp_id + 1 

    hout.close()


    hout = open(output_file + '.tmp.sorted.txt', 'w')
    subprocess.call(["sort", "-k1,1", "-k2,2n", "-k3,3n", output_file + ".tmp.unsorted.txt"], stdout = hout)
    hout.close()


    if control_file is not None:
        control_db = pysam.TabixFile(control_file)

    temp_chr = ""
    temp_start = ""
    temp_end = ""
    temp_count = ["0"] * temp_id
    hout = open(output_file, 'w')
    with open(output_file + '.tmp.sorted.txt', 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')

            if F[1] != temp_start or F[2] != temp_end or F[0] != temp_chr:

                # if not the first line 
                if temp_chr != "":

                    # skip if the junction is included in the control file
                    tabixErrorFlag = 0
                    if control_file is not None:
                        try:
                            records = control_db.fetch(temp_chr, int(temp_start) - 5, int(temp_start) + 5)
                        except Exception as inst:
                            # print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
                            # tabixErrorMsg = str(inst.args)
                            tabixErrorFlag = 1

                    control_flag = 0;
                    if tabixErrorFlag == 0:
                        for record_line in records:
                            record = record_line.split('\t')
                            if temp_chr == record[0] and temp_start == record[1] and temp_end == record[2]:
                                control_flag = 1

                    if control_flag == 0:
                        print >> hout, temp_chr + '\t' + temp_start + '\t' + temp_end + '\t' + ','.join(temp_count)

                temp_chr = F[0]
                temp_start = F[1]
                temp_end = F[2]
                temp_count = ["0"] * temp_id

            if F[1] == "10003993" and F[2] == "10009695":
                print '\t'.join(F)
                print ','.join(temp_count)

            temp_count[int(F[3])] = F[4]

    # last check 

    # skip if the junction is included in the control file
    tabixErrorFlag = 0
    if control_file is not None:
        try:
            records = control_db.fetch(temp_chr, int(temp_start) - 5, int(temp_start) + 5)
        except Exception as inst:
            # print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
            # tabixErrorMsg = str(inst.args)
            tabixErrorFlag = 1
            
    control_flag = 0;
    if tabixErrorFlag == 0:
        for record_line in records:
            record = record_line.split('\t')
            if temp_chr == record[0] and temp_start == record[1] and temp_end == record[2]:
                control_flag = 1
                
    if control_flag == 0:
        print >> hout, temp_chr + '\t' + temp_start + '\t' + temp_end + '\t' + ','.join(temp_count)

    hout.close()
 
    # remove intermediate files
    subprocess.call(["rm", "-rf", output_file + ".tmp.unsorted.txt"])
    subprocess.call(["rm", "-rf", output_file + ".tmp.sorted.txt"])

    if control_file is not None:
        control_db.close()


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



def merge_intron_retention(IR_file_list, output_file, control_file, ratio_thres, num_thres):

    # list up junctions to pick up
    intron_retention2list = {}
    header2ind = {}
    target_header = ["Chr", "Boundary_Pos", "Gene_Symbol", "Motif_Type", "Strand",
                     "Junction_List", "Gene_ID_List", "Exon_Num_List"]

    for IR_file in IR_file_list:
        with open(IR_file, 'r') as hin:
            header = hin.readline().rstrip('\n').split('\t')
            for (i, cname) in enumerate(header):
                header2ind[cname] = i

            for line in hin:
                F = line.rstrip('\n').split('\t')
                if int(F[header2ind["Intron_Retention_Read_Count"]]) < num_thres: continue
                ratio = 0
                if F[header2ind["Edge_Read_Count"]] != "0":
                    ratio = float(F[header2ind["Intron_Retention_Read_Count"]]) / float(F[header2ind["Edge_Read_Count"]])
                if ratio < ratio_thres: continue

                key = '\t'.join([F[header2ind[x]] for x in target_header])

                if key not in intron_retention2list: intron_retention2list[key] = 1


    temp_id = 0
    hout = open(output_file + ".tmp.unsorted.txt", 'w')
    for IR_file in IR_file_list:
        with open(IR_file, 'r') as hin:
            for line in hin:
                F = line.rstrip('\n').split('\t')
                key = '\t'.join([F[header2ind[x]] for x in target_header])
                if key in intron_retention2list:
                    print >> hout, key + '\t' + str(temp_id) + '\t' + F[header2ind["Intron_Retention_Read_Count"]]

        temp_id = temp_id + 1

    hout.close()


    hout = open(output_file + '.tmp.sorted.txt', 'w')
    subprocess.call(["sort", "-k1,1", "-k2,2n", output_file + ".tmp.unsorted.txt"], stdout = hout)
    hout.close()

    if control_file is not None:
        control_db = pysam.TabixFile(control_file)

    temp_chr = ""
    temp_pos = ""
    temp_key = ""
    temp_count = ["0"] * temp_id
    hout = open(output_file, 'w')
    print >> hout, '\t'.join(target_header) + '\t' + "Read_Count_Vector"

    with open(output_file + '.tmp.sorted.txt', 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')

            if F[1] != temp_pos or F[0] != temp_chr:

                # if not the first line 
                if temp_chr != "":

                    # skip if the junction is included in the control file
                    tabixErrorFlag = 0
                    if control_file is not None:
                        try:
                            records = control_db.fetch(temp_chr, int(temp_pos) - 5, int(temp_pos) + 5)
                        except Exception as inst:
                            # print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
                            # tabixErrorMsg = str(inst.args)
                            tabixErrorFlag = 1

                    control_flag = 0;
                    if tabixErrorFlag == 0:
                        for record_line in records:
                            record = record_line.split('\t')
                            if temp_chr == record[0] and temp_pos == record[1]:
                                control_flag = 1

                    if control_flag == 0:
                        print >> hout, temp_key + '\t' + ','.join(temp_count)

                temp_chr = F[0]
                temp_pos = F[1]
                temp_key = '\t'.join([F[header2ind[x]] for x in target_header])
                temp_count = ["0"] * temp_id

            temp_count[int(F[8])] = F[9]


    # last check 
    # skip if the junction is included in the control file
    tabixErrorFlag = 0
    if control_file is not None:
        try:
            records = control_db.fetch(temp_chr, int(temp_pos) - 5, int(temp_pos) + 5)
        except Exception as inst:
            # print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
            # tabixErrorMsg = str(inst.args)
            tabixErrorFlag = 1

    control_flag = 0;
    if tabixErrorFlag == 0:
        for record_line in records:
            record = record_line.split('\t')
            if temp_chr == record[0] and temp_pos == record[1]:
                control_flag = 1

    if control_flag == 0:
        print >> hout, '\t'.join(F[:8]) + '\t' + ','.join(temp_count)

    hout.close()
 
    # remove intermediate files
    subprocess.call(["rm", "-rf", output_file + ".tmp.unsorted.txt"])
    subprocess.call(["rm", "-rf", output_file + ".tmp.sorted.txt"])

    if control_file is not None:
        control_db.close()





def merge_mut(mutation_file_list, output_file):


    mut2sample = {}
    sample_num = "1"
    for mut_file in mutation_file_list:
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


def merge_SJ_IR_files(SJ_input_file, IR_input_file, output_file):

    header2ind = {}
    header = ""
    hout = open(output_file + ".unsorted", 'w')

    with open(SJ_input_file, 'r') as hin:

        header = hin.readline().rstrip('\n').split('\t')
        for i in range(len(header)):
            header2ind[header[i]] = i

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

            splicing_key = F[header2ind["SJ_1"]] + ':' + F[header2ind["SJ_2"]] + '-' + F[header2ind["SJ_3"]]

            print >> hout, gene + '\t' + splicing_key + '\t' + F[header2ind["Splicing_Class"]] + '\t' + F[header2ind["Is_Inframe"]] + '\t' + \
                            F[header2ind["SJ_4"]] + '\t' + F[header2ind["Mutation_Key"]] + '\t' + F[header2ind["Motif_Pos"]] + '\t' + \
                            F[header2ind["Mutation_Type"]] + '\t' + F[header2ind["Is_Canonical"]]


    with open(IR_input_file, 'r') as hin:

        header = hin.readline().rstrip('\n').split('\t')
        for i in range(len(header)):
            header2ind[header[i]] = i

        for line in hin:
            F = line.rstrip('\n').split('\t')
            splicing_key = F[header2ind["Chr"]] + ':' + F[header2ind["Boundary_Pos"]] + '-' + F[header2ind["Boundary_Pos"]]
            splicing_class = "intron-retention" if F[header2ind["Intron_Retention_Type"]] == "direct-impact" else "opposite-side-intron-retention"
            print >> hout, F[header2ind["Gene_Symbol"]] + '\t' + splicing_key + '\t' + splicing_class + '\t' + '---' + '\t' + \
                           F[header2ind["Read_Count_Vector"]] + '\t' + F[header2ind["Mutation_Key"]] + '\t' + \
                           F[header2ind["Motif_Pos"]] + '\t' + F[header2ind["Mutation_Type"]] + '\t' + F[header2ind["Is_Canonical"]]

    hout.close()

    hout = open(output_file, 'w')
    print >> hout, '\t'.join(["Gene_Symbol", "Splicing_Key", "Splicing_Class", "Is_Inframe", "Read_Counts",
                              "Mutation_Key", "Motif_Pos", "Mutation_Type", "Is_Canonical"])
    hout.close()

    hout = open(output_file, 'a')
    subprocess.call(["sort", "-k1", output_file + ".unsorted"], stdout = hout)
    hout.close()

    subprocess.call(["rm", "-rf", output_file + ".unsorted"])


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
    tchr, tstart, tend, tref, talt = mut_infos[0], mut_infos[1], mut_infos[1], mut_infos[2], mut_infos[3]
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


def organize_mut_splicing_count(input_file, mut2sample_file, output_count_file, output_mut_file, output_splicing_file):

    hout1 = open(output_count_file, 'w')
    hout2 = open(output_mut_file, 'w')
    hout3 = open(output_splicing_file, 'w')
        
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
        temp_splicing2info = {} 
        temp_splicing_count = []
        temp_splicing2id = {} 
        temp_id2splicing = {}
        temp_mut_splicing_link = []
        temp_mut_id = "1"
        temp_splicing_id = "1"
        for line in hin:

            F = line.rstrip('\n').split('\t')

            if F[header2ind["Gene_Symbol"]] != temp_gene:
                if temp_gene != "":

                    # flush the result
                    print >> hout1, temp_gene + '\t' + ';'.join(temp_mut2sample) + '\t' + ';'.join(temp_splicing_count) + '\t' + ';'.join(temp_mut_splicing_link)
                
                    for id in sorted(temp_id2mut):
                        print >> hout2, temp_gene + '\t' + id + '\t' + temp_id2mut[id] + '\t' + ';'.join(temp_mut2info[temp_id2mut[id]])

                    for id in sorted(temp_id2splicing):
                        print >> hout3, temp_gene + '\t' + id + '\t' + temp_id2splicing[id] + '\t' + temp_splicing2info[temp_id2splicing[id]]

                temp_gene = F[header2ind["Gene_Symbol"]] 
                temp_mut2info = {} 
                temp_mut2sample = []
                temp_mut2id = {}
                temp_id2mut = {}
                temp_splicing2info = {} 
                temp_splicing_count = []
                temp_splicing2id = {} 
                temp_id2splicing = {}
                temp_mut_splicing_link = []
                temp_mut_id = "1"
                temp_splicing_id = "1"


            mut = F[header2ind["Mutation_Key"]]
            sample = get_mut_sample_info(F[header2ind["Mutation_Key"]], mut2sample)
            mut_info = F[header2ind["Motif_Pos"]] + ',' + F[header2ind["Mutation_Type"]] + ',' + F[header2ind["Is_Canonical"]]

            if mut not in temp_mut2info:
                temp_mut2info[mut] = []
                temp_mut2id[mut] = temp_mut_id
                temp_id2mut[temp_mut_id] = mut
                temp_mut2sample.append(temp_mut_id + ':' + sample)
                temp_mut_id = str(int(temp_mut_id) + 1)

            if mut_info not in temp_mut2info[mut]: temp_mut2info[mut].append(mut_info) 


            splicing_key = F[header2ind["Splicing_Key"]]
            if splicing_key not in temp_splicing2info:
                temp_splicing2info[splicing_key] = '\t'.join([F[header2ind[x]] for x in ["Splicing_Class", "Is_Inframe"]])
                temp_splicing2id[splicing_key] = temp_splicing_id
                temp_id2splicing[temp_splicing_id] = splicing_key 
                temp_splicing_count.append(F[header2ind["Read_Counts"]])
                temp_splicing_id = str(int(temp_splicing_id) + 1)

            if temp_mut2id[mut] + ',' + temp_splicing2id[splicing_key] not in temp_mut_splicing_link:
                temp_mut_splicing_link.append(temp_mut2id[mut] + ',' + temp_splicing2id[splicing_key])


    # last flush
    print >> hout1, temp_gene + '\t' + ';'.join(temp_mut2sample) + '\t' + ';'.join(temp_splicing_count) + '\t' + ';'.join(temp_mut_splicing_link)

    for id in sorted(temp_id2mut):
        print >> hout2, temp_gene + '\t' + id + '\t' + temp_id2mut[id] + '\t' + ';'.join(temp_mut2info[temp_id2mut[id]]) 

    for id in sorted(temp_id2splicing):
        print >> hout3, temp_gene + '\t' + id + '\t' + temp_id2splicing[id] + '\t' + temp_splicing2info[temp_id2splicing[id]]

    hout1.close()
    hout2.close()
    hout3.close()


def simple_link_effect_check(mutation_state, splicing_count, link, weight_vector, pseudo_count = 0.1):

    mutation_states = mutation_state.split(';')
    splicing_counts = splicing_count.split(';')
    link_vector = link.split(';')

    sample_num = len(splicing_counts[0].split(','))

    """
    mut_vector = [0] * sample_num
    for j in range(len(mutation_states)):
        mut_id, sample_id_str = mutation_states[j].split(':')
        for sample_id in sample_id_str.split(','):
            mut_vector[int(sample_id) - 1] = int(mut_id)
    """

    # pass_links 
    effect_size_vector = [0] * len(link_vector)

    # simple check for each link
    for i in range(len(link_vector)):
        mut_id, sp_id = link_vector[i].split(',')
        splicing_cont_vector = splicing_counts[int(sp_id) - 1].split(',') 

        # extract samples with the mutation of the link in consideration
        mut_vector = [0] * sample_num
        for j in range(len(mutation_states)):
            tmut_id, sample_id_str = mutation_states[j].split(':')
            if tmut_id == mut_id:  
                for sample_id in sample_id_str.split(','):
                    mut_vector[int(sample_id) - 1] = 1

        weight_sum_null = sum([float(weight_vector[j]) for j in range(sample_num) if mut_vector[j] == 0])
        weight_sum_target = sum([float(weight_vector[j]) for j in range(sample_num) if mut_vector[j] == 1])

        count_sum_null = sum([int(splicing_cont_vector[j]) for j in range(sample_num) if mut_vector[j] == 0])
        count_sum_target = sum([int(splicing_cont_vector[j]) for j in range(sample_num) if mut_vector[j] == 1])

        mean_null = float(count_sum_null) / weight_sum_null if weight_sum_null > 0.0 and count_sum_null > 0 else 0.0
        mean_target = float(count_sum_target) / weight_sum_target if weight_sum_target > 0 and count_sum_target > 0 else 0.0
    
        effect_size_vector[i] = (float(mean_target) + pseudo_count) / (float(mean_null) + pseudo_count)

    return effect_size_vector 


def cluster_link(link):

    link_vector = link.split(';')
    mut_id2sp_id = {}
    for i in range(len(link_vector)):
        mut_id, sp_id = link_vector[i].split(',')
        if mut_id not in mut_id2sp_id: mut_id2sp_id[mut_id] = []
        mut_id2sp_id[mut_id].append(sp_id)

    link_str = []
    mut_ids = mut_id2sp_id.keys()
    active_ids = copy.deepcopy(mut_ids) # deep copy
    for i in range(len(mut_ids)):

        if mut_ids[i] not in active_ids: continue

        clustered_ids = [mut_ids[i]]

        no_more_cluster = 0
        while no_more_cluster == 0:
            
            no_more_cluster = 1
            for j in range(len(active_ids)):
                if active_ids[j] in clustered_ids: continue
                is_cluster = 0
                for k in range(len(clustered_ids)):
                    if len(set(mut_id2sp_id[clustered_ids[k]]) & set(mut_id2sp_id[active_ids[j]])) > 0:
                        clustered_ids.append(active_ids[j])
                        no_more_cluster = 0
                        break
 
        link_strs = []
        for j in range(len(clustered_ids)):
            link_strs.append(';'.join([clustered_ids[j] + ',' + x for x in mut_id2sp_id[clustered_ids[j]]]))
            active_ids.remove(clustered_ids[j])
            
        link_str.append(';'.join(link_strs))

    return link_str


def convert_pruned_file(input_file, output_file, weight_vector, margin):

    hout = open(output_file, 'w')
    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            gene = F[0]
            mutation_state = F[1]
            splicing_count = F[2]
            link = F[3]

            effect_size_vector = simple_link_effect_check(mutation_state, splicing_count, link, weight_vector)
            
            # print F[0]
            # print '\t'.join([str(x) for x in effect_size_vector])
 
            link_vector = link.split(';')

            link2effect_size = {}
            for i in range(len(link_vector)):
                link2effect_size[link_vector[i]] = effect_size_vector[i]

            pass_link = [link_vector[i] for i in range(len(link_vector)) if effect_size_vector[i] >= margin]

            if len(pass_link) > 0:
                clustered_sets = cluster_link(';'.join(pass_link))
                for sub_cluster in sorted(clustered_sets):
                    sub_cluster_link = sub_cluster.split(';')
                    # maximum number is 10
                    if len(sub_cluster_link) >= 10:

                        sub_cluster_effect_size_vector = [link2effect_size[x] for x in sub_cluster_link]
                        temp_margin = sorted(sub_cluster_effect_size_vector, reverse=True)[9]
                        sub_cluster_link = [sub_cluster_link[i] for i in range(len(sub_cluster_link)) if sub_cluster_effect_size_vector[i] >= temp_margin]

                    print >> hout, gene + '\t' + mutation_state + '\t' + splicing_count + '\t' + ';'.join(sub_cluster_link)

    hout.close()



def get_BIC(mutation_state, splicing_count, configuration, link):

    mutation_states = mutation_state.split(';')
    splicing_counts = splicing_count.split(';')
    configuration_vector = configuration.split(',')
    link_vector = link.split(';')

    # get the number of possible mutations
    possible_mut_num = 1
    possible_mutation = [0]
    for i in range(len(mutation_states)):
        mut_id, sample_id = mutation_states[i].split(':')
        if mut_id not in possible_mutation: possible_mutation.append(int(mut_id))

    possible_mut_num = len(possible_mutation)


    sample_num = len(splicing_counts[0].split(','))

    loglikelihood = 0
    params = []
    param_num = 0
    for i in range(len(splicing_counts)):

        splicing_cont_vector = splicing_counts[i].split(',')
        active_mut_list = []
        # get mutations associated with the current splicing
        for j in range(len(link_vector)):
            if int(configuration_vector[j]) == 0: continue
            mut_id, sp_id = link_vector[j].split(',')
            if int(sp_id) == i + 1:
                active_mut_list.append(int(mut_id))

        param_num = param_num + len(active_mut_list) + 1
        current_mut_vector = [0] * sample_num
    
        for j in range(len(mutation_states)):
            mut_id, sample_id_str = mutation_states[j].split(':')
            if int(mut_id) in active_mut_list:
                for sample_id in sample_id_str.split(','):
                    current_mut_vector[int(sample_id) - 1] = int(mut_id)

        current_params = ["0.000"] * possible_mut_num
        for k in range(possible_mut_num):
            sample_sum = len([j for j in range(sample_num) if current_mut_vector[j] == k])
            count_sum = sum([int(splicing_cont_vector[j]) for j in range(sample_num) if current_mut_vector[j] == k])
            if sample_sum > 0: current_params[k] = str(round(float(count_sum) / sample_sum, 3))
            if count_sum > 0:
                loglikelihood = loglikelihood - count_sum * (1 + math.log(sample_sum)) + count_sum * math.log(count_sum)

        params.append(','.join(current_params))

    BIC = round(-2 * loglikelihood + 2 * math.log(sample_num * len(splicing_counts)) * param_num, 4)

    return([BIC, ';'.join(params)])
   
     

def get_log_marginal_likelihood(mutation_state, splicing_count, configuration_vector, link, weight_vector, alpha0, beta0, alpha1, beta1):

    mutation_states = mutation_state.split(';')
    splicing_counts = splicing_count.split(';')
    link_vector = link.split(';')

    sample_num = len(splicing_counts[0].split(','))

    log_marginal_likelihood = 0
    for i in range(len(splicing_counts)):

        splicing_cont_vector = splicing_counts[i].split(',')
        active_mut_list = []

        # get mutations associated with the current splicing
        for e in range(len(link_vector)):
            if int(configuration_vector[e]) == 0: continue
            mut_id, sp_id = link_vector[e].split(',')
            if int(sp_id) == i + 1:
                active_mut_list.append(int(mut_id))

        # param_num = param_num + len(active_mut_list) + 1
        active_mut_vector = [0] * sample_num

        # set mutation status
        for n in range(len(mutation_states)):
            mut_id, sample_id_str = mutation_states[n].split(':')
            if int(mut_id) in active_mut_list:
                for sample_id in sample_id_str.split(','):
                    active_mut_vector[int(sample_id) - 1] = 1

        # for inactive mutations
        sample_sum = len([n for n in range(sample_num) if active_mut_vector[n] == 0])
        count_sum = sum([int(splicing_cont_vector[n]) for n in range(sample_num) if active_mut_vector[n] == 0])
        weight_sum = sum([weight_vector[n] for n in range(sample_num) if active_mut_vector[n] == 0])
        if sample_sum > 0:
            log_marginal_likelihood = log_marginal_likelihood + math.lgamma(count_sum + alpha0) - math.lgamma(alpha0) + \
                                      alpha0 * math.log(beta0) - (count_sum + alpha0) * math.log(weight_sum + beta0)


        for n in range(sample_num):
            if active_mut_vector[n] == 0: continue
            log_marginal_likelihood = log_marginal_likelihood + math.lgamma(int(splicing_cont_vector[n]) + alpha1) - math.lgamma(alpha1) + \
                                      alpha1 * math.log(beta1) - (int(splicing_cont_vector[n]) + alpha1) * math.log(weight_vector[n] + beta1)

    return(log_marginal_likelihood)



def generate_configurations(dim):

    conf = [[0], [1]]

    for k in range(dim - 1):

        new_conf = []
        for elm in conf:
            new_conf.append(elm + [0])
            new_conf.append(elm + [1])
            conf = new_conf

    return conf


def check_significance(input_file, output_file, weight_vector, log_BF_thres, alpha0, beta0, alpha1, beta1):

    def soft_max(x):
        x_max = max(x)
        return(x_max + math.log(sum([math.exp(y - x_max) for y in x]))) 

    hout = open(output_file, 'w')
    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            gene = F[0]
            mutation_state = F[1]
            splicing_count = F[2]
            link = F[3]

            conf_dim = len(link.split(';'))

            # log_MLs_nonnull = []
            # log_ML_null = float("-inf")
            # log_ML_max = float("-inf")
            # conf_max = [0] * conf_dim
            params_min = "---"
            conf_num = 0

            mut2log_ML_null = {}
            mut2log_ML_nonnull = {}
            mut2log_ML_nonnull_max = {}
            mut2conf_max = {}
            mut2log_BF = {}

            link_vector = link.split(';')
            # get mutations associated with the current splicing
            for e in range(len(link_vector)):
                mut_id, sp_id = link_vector[e].split(',')
                mut2log_ML_null[mut_id] = [] 
                mut2log_ML_nonnull[mut_id] = []
                mut2log_ML_nonnull_max[mut_id] = float("-inf") 
                mut2conf_max[mut_id] = [0] * conf_dim
                mut2log_BF[mut_id] = float("-inf")
 
            for conf in sorted(generate_configurations(conf_dim)):

                # get active mutation in the configuration in consideration
                active_mut_list = []
                # get mutations associated with the current splicing
                for j in range(len(link_vector)):
                    if int(conf[j]) == 0: continue
                    mut_id, sp_id = link_vector[j].split(',')
                    active_mut_list.append(mut_id)


                conf_num = conf_num + 1
                log_ML = get_log_marginal_likelihood(mutation_state, splicing_count, conf, link,
                                                     weight_vector, alpha0, beta0, alpha1, beta1)

                for mut_id in mut2log_ML_null:
                    if mut_id in active_mut_list:
                        mut2log_ML_nonnull[mut_id].append(log_ML)
                        if log_ML > mut2log_ML_nonnull_max[mut_id]:
                            mut2log_ML_nonnull_max[mut_id] = log_ML
                            mut2conf_max[mut_id] = conf
                    else:
                        mut2log_ML_null[mut_id].append(log_ML)

                """
                if conf == [0] * conf_dim:
                    log_ML_null = log_ML
                else:
                    log_MLs_nonnull.append(log_ML)

                if log_ML > log_ML_max:
                    log_ML_max = log_ML
                    conf_max = conf
                """

                """
                if gene == "DDX17":
                    print conf
                    print log_ML
                """

            # log_BF_sum = soft_max(log_MLs_nonnull) - math.log(float(conf_num - 1)) - log_ML_null

            for mut_id in mut2log_ML_null:

                mut2log_BF[mut_id] = soft_max(mut2log_ML_nonnull[mut_id]) - soft_max(mut2log_ML_null[mut_id]) - \
                                     math.log(len(mut2log_ML_nonnull[mut_id])) + math.log(len(mut2log_ML_null[mut_id]))

            # print gene + '\t' + str(conf_num) + '\t' + str(round((log_ML_max - log_ML_null), 4)) + '\t' + str(round(log_BF_sum, 4))

            for mut_id in mut2log_ML_null:

            # if log_BF_sum > log_BF_thres:
                print >> hout, '\t'.join(F) + '\t' + mut_id + '\t' + str(round(mut2log_BF[mut_id], 4)) + '\t' + ','.join([str(x) for x in mut2conf_max[mut_id]])
                # print >> hout, '\t'.join(F) + '\t' + str(log_BF_sum) + '\t' + ','.join([str(x) for x in conf_max]) 

            """
            # BIC based approach, deprecated 
            BIC0 = float("-inf") 
            BIC_min = float("inf")
            conf_min = [0] * conf_dim
            params_min = "---"
            for conf in sorted(generate_configurations(conf_dim)):

                BIC, params = get_BIC(mutation_state, splicing_count, ','.join([str(x) for x in conf]), link)
                if conf == [0] * conf_dim:
                    BIC0 = BIC
                    
                if BIC < BIC_min:
                    BIC_min = BIC
                    conf_min = conf
                    params_min = params

            if BIC0 - BIC_min > 10.0:
                print >> hout, '\t'.join(F) + '\t' + str(BIC_min) + '\t' + str(BIC0) + '\t' + ','.join([str(x) for x in conf_min]) + '\t' + params_min
            """

    hout.close()



def add_annotation(input_file, output_file, sample_name_list, mut_info_file, sp_info_file):

    id2sample = {}
    temp_id = "1"
    for sample_name in sample_name_list:
        id2sample[temp_id] = sample_name 
        temp_id = str(int(temp_id) + 1)

    """
    with open(sample_list_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            id2sample[temp_id] = F[0]
            temp_id = str(int(temp_id) + 1)
    """

    mut_id2mut_info = {}
    with open(mut_info_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')

            FF = F[3].split(';')
            info1, info2, info3 = [], [], []
            for i in range(len(FF)):
                FFF = FF[i].split(',')
                info1.append(FFF[0] + ',' + FFF[1])
                info2.append(FFF[2])
                info3.append(FFF[3])

            mut_id2mut_info[F[0] + '\t' + F[1]] = F[2] + '\t' + ';'.join(info1) + '\t' + ';'.join(info2) + '\t' + ';'.join(info3)

    sp_id2sp_info = {}
    with open(sp_info_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            sp_id2sp_info[F[0] + '\t' + F[1]] = F[2] + '\t' + F[3] + '\t' + F[4]

    hout = open(output_file, 'w')
    print >> hout, "Gene_Symbol" + '\t' + "Sample_Name" + '\t' + "Mutation_Key" + '\t' + "Motif_Pos" + '\t' + \
                   "Mutation_Type" + '\t' + "Is_Canonical" + '\t' + "Splicing_Key" +'\t' + "Splicing_Class" + '\t' + \
                   "Is_Inframe" + '\t' + "Score"

    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')

            gene = F[0]
            mutation_states = F[1].split(';')
            splicing_count_vector = F[2].split(';')
            link_vector = F[3].split(';')
            mut_id = F[4]
            log_BF = F[5]
            # log_ML0 = F[5]
            conf_max_vector = F[6].split(',')

            mut_id2sample_id = {}
            for mut_state in mutation_states:
                tmut_id, sample_ids = mut_state.split(':')
                mut_id2sample_id[tmut_id] = sample_ids


            active_link_vector = [link_vector[i] for i in range(len(link_vector)) if conf_max_vector[i] == "1"]

            for active_link in active_link_vector:
                tmut_id, tsp_id = active_link.split(',')
                if tmut_id != mut_id: continue
                current_splicing_counts = splicing_count_vector[int(tsp_id) - 1].split(',')                


                # get sample names
                sample_names = []
                for sample_id in mut_id2sample_id[mut_id].split(','):            
                    if int(current_splicing_counts[int(sample_id) - 1]) > 0:
                        sample_names.append(id2sample[sample_id])

        
                # get mutation info
                mut_info = mut_id2mut_info[gene + '\t' + mut_id]
            
                # get SJ info
                sp_info = sp_id2sp_info[gene + '\t' + tsp_id]

                print >> hout, gene + '\t' + ';'.join(sample_names) + '\t' + mut_info + '\t' + sp_info + '\t' + \
                               str(round(float(log_BF), 4))

    hout.close()


def permute_mut_SJ_pairs(input_file, output_file):

    sample_num = ""
    permute = {}
    hout = open(output_file, 'w')
    with open(input_file, 'r') as hin:
        for line in hin:
    
            F = line.rstrip('\n').split('\t')
            gene = F[0]
            mutation_state = F[1]
            splicing_count = F[2]
            link = F[3]

            splicing_count_vector = F[2].split(';')
            if sample_num == "":
                sample_num = len(splicing_count_vector[0].split(','))
                shuffle_order = range(sample_num)
                no_overlap = 0
                while no_overlap == 0:
                    random.shuffle(shuffle_order)
                    overlap_num = sum([shuffle_order[i] == i for i in range(sample_num)])
                    if overlap_num == 0: no_overlap = 1
 
                for i in range(sample_num):
                    permute[str(i + 1)] = str(shuffle_order[i] + 1)

            permute_mutation_states = []
            for mut_sample in mutation_state.split(';'):
                mut_id, sample_ids = mut_sample.split(':')
                permute_mutation_states.append(mut_id + ':' + ','.join([permute[x] for x in sample_ids.split(',')]))

            print >> hout, gene + '\t' + ';'.join(permute_mutation_states) + '\t' + splicing_count + '\t' + link
    
    hout.close()


def calculate_q_value(input_file, permutation_file_prefix, output_file, permutation_num):

    logBF_values_null = []
    for i in range(permutation_num):
        header2ind = {}
        mut2logBF = {}
        with open(permutation_file_prefix + str(i) + ".txt", 'r') as hin:
            header = hin.readline().rstrip('\n').split('\t')
            for (i, cname) in enumerate(header):
                header2ind[cname] = i
            for line in hin:
                F = line.rstrip('\n').split('\t')
                mut2logBF[F[header2ind["Mutation_Key"]]] = float(F[header2ind["Score"]])

        logBF_values_null = logBF_values_null + mut2logBF.values()

    # count total mutation num
    mut_keys = []
    with open(input_file, 'r') as hin:
        header = hin.readline().rstrip('\n').split('\t')
        for (i, cnname) in enumerate(header):
            header2ind[cname] = i
        for line in hin:
            F = line.rstrip('\n').split('\t')
            mut_keys.append(F[header2ind["Mutation_Key"]])

    mut_keys = list(set(mut_keys))

    # add q-values for each key
    hout = open(output_file, 'w')
    header2ind = {}
    with open(input_file, 'r') as hin:
        header = hin.readline().rstrip('\n').split('\t')
        for (i, cnname) in enumerate(header):
            header2ind[cname] = i
        print >> hout, '\t'.join(header) + '\t' + "Q_Value"
        for line in hin:
            F = line.rstrip('\n').split('\t')
            qvalue = float(len([x for x in logBF_values_null if x >= float(F[header2ind["Score"]])])) / (len(mut_keys) * permutation_num)

            print >> hout, '\t'.join(F) + '\t' + str(round(qvalue, 4))

    hout.close()


