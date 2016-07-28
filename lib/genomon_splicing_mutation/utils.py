#! /usr/bin/env python

import sys, glob, subprocess, re, math, pysam

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
            mut_id2mut_info[F[0] + '\t' + F[1]] = F[2] + '\t' + F[3]

    SJ_id2SJ_info = {}
    with open(SJ_info_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            SJ_id2SJ_info[F[0] + '\t' + F[1]] = F[2] + '\t' + F[3] + '\t' + F[4]

    hout = open(output_file, 'w')
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

                print >> hout, gene + '\t' + ';'.join(sample_names) + '\t' + mut_info + '\t' + SJ_info

    hout.close()

