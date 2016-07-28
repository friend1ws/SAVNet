#! /usr/bin/env python

import sys, math

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


input_file = sys.argv[1]

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
            print '\t'.join(F) + '\t' + str(BIC_min) + '\t' + str(BIC0) + '\t' + ','.join([str(x) for x in conf_min]) + '\t' + params_min


