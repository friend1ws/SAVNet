#! /usr/bin/env python

import sys, glob, subprocess, re, math, copy, random, pysam

def merge_SJ2(sample_list_file, output_file, control_file, junc_num_thres):

    # list up junctions to pick up
    junc2list = {}
    with open(sample_list_file, 'r') as hin1:
        for line1 in hin1:
            sample_name, mut_file, SJ_file = line1.rstrip('\n').split('\t')
            with open(SJ_file, 'r') as hin2:
                for line2 in hin2:
                    F = line2.rstrip('\n').split('\t')
                    if F[5] != "0": continue
                    if int(F[6]) < junc_num_thres: continue
                    key = F[0] + '\t' + F[1] + '\t' + F[2]
                    if key not in junc2list: junc2list[key] = 1


    temp_id = 0
    hout = open(output_file + ".tmp.unsorted.txt", 'w')
    with open(sample_list_file, 'r') as hin1:
        for line1 in hin1:
            sample_name, mut_file, SJ_file = line1.rstrip('\n').split('\t')
            with open(SJ_file, 'r') as hin2:
                for line2 in hin2:
                    F = line2.rstrip('\n').split('\t')
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


def simple_link_effect_check(mutation_state, splicing_count, link, pseudo_count = 0.1):

    mutation_states = mutation_state.split(';')
    splicing_counts = splicing_count.split(';')
    link_vector = link.split(';')

    sample_num = len(splicing_counts[0].split(','))

    mut_vector = [0] * sample_num
    for j in range(len(mutation_states)):
        mut_id, sample_id_str = mutation_states[j].split(':')
        for sample_id in sample_id_str.split(','):
            mut_vector[int(sample_id) - 1] = int(mut_id)

    # pass_links 
    effect_size_vector = [0] * len(link_vector)

    # simple check for each link
    for i in range(len(link_vector)):
        mut_id, sp_id = link_vector[i].split(',')
        splicing_cont_vector = splicing_counts[int(sp_id) - 1].split(',') 
    
        sample_sum_null = len([j for j in range(sample_num) if mut_vector[j] == 0])
        sample_sum_target = len([j for j in range(sample_num) if mut_vector[j] == int(mut_id)])
    
        count_sum_null = sum([int(splicing_cont_vector[j]) for j in range(sample_num) if mut_vector[j] == 0])
        count_sum_target = sum([int(splicing_cont_vector[j]) for j in range(sample_num) if mut_vector[j] == int(mut_id)])
 
        mean_null = float(count_sum_null) / sample_sum_null if sample_sum_null > 0 and count_sum_null > 0 else 0.0
        mean_target = float(count_sum_target) / sample_sum_target if sample_sum_target > 0 and count_sum_target > 0 else 0.0

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


def convert_pruned_file(input_file, output_file, margin):

    hout = open(output_file, 'w')
    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            gene = F[0]
            mutation_state = F[1]
            splicing_count = F[2]
            link = F[3]

            effect_size_vector = simple_link_effect_check(mutation_state, splicing_count, link)
            
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
   
     

def generate_configurations(dim):

    conf = [[0], [1]]

    for k in range(dim - 1):

        new_conf = []
        for elm in conf:
            new_conf.append(elm + [0])
            new_conf.append(elm + [1])
            conf = new_conf

    return conf


def check_significance(input_file, output_file):

    hout = open(output_file, 'w')
    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            gene = F[0]
            mutation_state = F[1]
            splicing_count = F[2]
            link = F[3]

            conf_dim = len(link.split(';'))

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

    hout.close()



def summarize_result(input_file, output_file, sample_list_file, mut_info_file, SJ_info_file):

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

            FF = F[3].split(';')
            info1, info2, info3 = [], [], []
            for i in range(len(FF)):
                FFF = FF[i].split(',')
                info1.append(FFF[0] + ',' + FFF[1])
                info2.append(FFF[2])
                info3.append(FFF[3])

            mut_id2mut_info[F[0] + '\t' + F[1]] = F[2] + '\t' + ';'.join(info1) + '\t' + ';'.join(info2) + '\t' + ';'.join(info3)

    SJ_id2SJ_info = {}
    with open(SJ_info_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            SJ_id2SJ_info[F[0] + '\t' + F[1]] = F[2] + '\t' + F[3] + '\t' + F[4]

    hout = open(output_file, 'w')
    print >> hout, "Gene_Symbol" + '\t' + "Sample_Name" + '\t' + "Mutation_Pos" + '\t' + "Splicing_Motif_Pos" + '\t' + \
                   "Mutation_Type1" + '\t' + "Mutation_Type2" + '\t' + "Splicing_Pos" +'\t' + "Splicing_Type" + '\t' + \
                   "Is_Inframe" + '\t' + "Score"

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

                print >> hout, gene + '\t' + ';'.join(sample_names) + '\t' + mut_info + '\t' + SJ_info + '\t' + \
                               str(round(float(BIC0) - float(BIC_min), 4))

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


