#! /usr/bin/env python

from __future__ import print_function

import copy, math, random

from . import utils
from .sav import Sav

class Network(object):

    def __init__(self, gene, mutation_status, splicing_counts, link2info, sample_list, weight_vector):
        self.gene = gene
        self.mutation_status = mutation_status
        self.splicing_counts = splicing_counts
        self.link_vector = list(link2info)
        self.link2info = link2info
        self.sample_list = sample_list
        self.sample_num = len(sample_list)
        self.weight_vector = weight_vector

        self.pruned_link_vector = []

        self.link_vector2effect_size = {}
        # self.link_vector2median_count = {}

        self.link_vector2active_quartile_count = {}
        self.link_vector2inactive_quartile_count = {}

        self.clustered_link_vector = []

        self.mut2log_BF = {}
        self.mut2significant_links = {}


    def prune_link_vector(self, margin, active_zero_filter_prob, inactive_nonzero_filter_prob):

        self.__link_effect_size_scan()
        # self.__link_median_count_scan()
        self.__link_active_quartile_count_scan(active_zero_filter_prob)
        self.__link_inactive_quartile_count_scan(inactive_nonzero_filter_prob)
        # self.pruned_link_vector = [x for x in self.link_vector if self.link_vector2effect_size[x] >= margin and self.link_vector2median_count[x] == 0]
        self.pruned_link_vector = [x for x in self.link_vector if self.link_vector2effect_size[x] >= margin and \
                                                                  self.link_vector2active_quartile_count[x] > 0 and \
                                                                  self.link_vector2inactive_quartile_count[x] == 0]


    def cluster_link_vector(self):

        for sub_cluster_link in self.__simple_cluster_link():

            if len(sub_cluster_link) > 15:

                sub_cluster_effect_size_vector = [self.link_vector2effect_size[x] for x in sub_cluster_link]
                temp_margin = sorted(sub_cluster_effect_size_vector, reverse=True)[15]
                sub_cluster_link = [sub_cluster_link[i] for i in range(len(sub_cluster_link)) if sub_cluster_effect_size_vector[i] > temp_margin]

            # may be removed at the above filtering step
            if len(sub_cluster_link) > 0:
                self.clustered_link_vector.append(sub_cluster_link)

    
    def get_averaged_bayes_factors(self, alpha0, beta0, alpha1, beta1):

        for sub_cluster_link in self.clustered_link_vector:

            conf_dim = len(sub_cluster_link)

            mut2log_ML_null = {}
            mut2log_ML_nonnull = {}
            mut2log_ML_nonnull_max = {}
            mut2conf_max = {}
            mut2log_BF = {}

            # get mutations associated with the current splicing
            for mut_id, _ in sub_cluster_link:
                mut2log_ML_null[mut_id] = [] 
                mut2log_ML_nonnull[mut_id] = []
                mut2log_ML_nonnull_max[mut_id] = float("-inf") 
                mut2conf_max[mut_id] = [0] * conf_dim
                mut2log_BF[mut_id] = float("-inf")


            for conf in sorted(utils.generate_configurations(conf_dim)):

                # get active mutation in the configuration in consideration
                active_mut_list = []
                # get mutations associated with the current splicing
                for j in range(len(sub_cluster_link)):
                    if int(conf[j]) == 0: continue
                    mut_id, sp_id = sub_cluster_link[j]
                    active_mut_list.append(mut_id)

                log_ML = self.__get_log_marginal_likelihood(sub_cluster_link, conf, alpha0, beta0, alpha1, beta1)

                for mut_id in mut2log_ML_null:
                    if mut_id in active_mut_list:
                        mut2log_ML_nonnull[mut_id].append(log_ML)
                        if log_ML > mut2log_ML_nonnull_max[mut_id]:
                            mut2log_ML_nonnull_max[mut_id] = log_ML
                            mut2conf_max[mut_id] = conf
                    else:
                        mut2log_ML_null[mut_id].append(log_ML)


            for mut_id in mut2log_ML_null:

                self.mut2log_BF[mut_id] = utils.soft_max(mut2log_ML_nonnull[mut_id]) - utils.soft_max(mut2log_ML_null[mut_id]) - \
                                            math.log(len(mut2log_ML_nonnull[mut_id])) + math.log(len(mut2log_ML_null[mut_id]))

                max_conf = mut2conf_max[mut_id]
                significant_links = [sub_cluster_link[i] for i in range(len(max_conf)) if max_conf[i] == 1]
                self.mut2significant_links[mut_id] = significant_links


    def export_to_savs(self, log_BF_thres, sample_num_thres = 1):
   
        sav_list = []
        all_active_sample_list = []

        for mut_id in self.mut2log_BF:
            if self.mut2log_BF[mut_id] < log_BF_thres: continue

            active_sample_list = [self.sample_list[i] for i in self.mutation_status[mut_id]]
            all_active_sample_list = all_active_sample_list + active_sample_list

            active_splicing_counts_vector = []
            active_link_info_vector = [] 
            for active_link in self.mut2significant_links[mut_id]:
                cur_mut_id, cur_sp_id = active_link
                if cur_mut_id != mut_id: continue

                cur_splicing_counts = self.splicing_counts[cur_sp_id]
                active_splicing_counts_vector.append([cur_splicing_counts[i] for i in self.mutation_status[mut_id]])
                active_link_info_vector.append(self.link2info[active_link])

            tsav = Sav(self.gene, active_sample_list, active_link_info_vector, active_splicing_counts_vector, self.mut2log_BF[mut_id])
            sav_list.append(tsav)

        all_active_sample_list = list(set(all_active_sample_list))
        if len(all_active_sample_list) >= sample_num_thres:
            return sav_list
        else:
            return []
               
        
    def execute_permutation(self, seed):

        for mutation in self.mutation_status:
            shuffled_samples = [seed[x] for x in self.mutation_status[mutation]]
            self.mutation_status[mutation] = shuffled_samples

        
    def __link_effect_size_scan(self, pseudo_count = 0.1):

        # simple check for each link
        for mut_id, sp_id in self.link_vector:
            cur_splicing_count_vector = self.splicing_counts[sp_id] 

            # extract samples with the mutation of the link in consideration
            mut_vector = [0] * self.sample_num
            for sample_id in self.mutation_status[mut_id]:
                mut_vector[sample_id] = 1

            weight_sum_null = sum([float(self.weight_vector[j]) for j in range(self.sample_num) if mut_vector[j] == 0])
            weight_sum_target = sum([float(self.weight_vector[j]) for j in range(self.sample_num) if mut_vector[j] == 1])

            count_sum_null = sum([int(cur_splicing_count_vector[j]) for j in range(self.sample_num) if mut_vector[j] == 0])
            count_sum_target = sum([int(cur_splicing_count_vector[j]) for j in range(self.sample_num) if mut_vector[j] == 1])

            mean_null = float(count_sum_null) / weight_sum_null if weight_sum_null > 0.0 and count_sum_null > 0 else 0.0
            mean_target = float(count_sum_target) / weight_sum_target if weight_sum_target > 0 and count_sum_target > 0 else 0.0
    
            self.link_vector2effect_size[(mut_id, sp_id)] = (float(mean_target) + pseudo_count) / (float(mean_null) + pseudo_count)



    def __link_active_quartile_count_scan(self, active_zero_filter_prob = 0.5):

        # simple check for each link
        for mut_id, sp_id in self.link_vector:
            # cur_mut_id, cur_sp_id = self.link_vector[i]
            cur_splicing_count_vector = self.splicing_counts[sp_id]

            # extract samples with the mutation of the link in consideration
            mut_vector = [0] * self.sample_num
            for sample_id in self.mutation_status[mut_id]:
                mut_vector[sample_id] = 1

            # self.link_vector2median_count[(mut_id, sp_id)] = utils.median([int(cur_splicing_count_vector[j]) for j in range(self.sample_num) if mut_vector[j] == 0])
            self.link_vector2active_quartile_count[(mut_id, sp_id)] = utils.quartile([int(cur_splicing_count_vector[j]) for j in range(self.sample_num) if mut_vector[j] == 1], active_zero_filter_prob)


    """
    # def __link_median_count_scan(self):
    def __link_inactive_quartile_count_scan(self, inactive_nonzero_filter_prob = 0.5):

        # simple check for each link
        for mut_id, sp_id in self.link_vector:
            # cur_mut_id, cur_sp_id = self.link_vector[i]
            cur_splicing_count_vector = self.splicing_counts[sp_id]

            # extract samples with the mutation of the link in consideration
            mut_vector = [0] * self.sample_num
            for sample_id in self.mutation_status[mut_id]:
                mut_vector[sample_id] = 1

            # self.link_vector2median_count[(mut_id, sp_id)] = utils.median([int(cur_splicing_count_vector[j]) for j in range(self.sample_num) if mut_vector[j] == 0])
            self.link_vector2inactive_quartile_count[(mut_id, sp_id)] = utils.quartile([int(cur_splicing_count_vector[j]) for j in range(self.sample_num) if mut_vector[j] == 0], inactive_nonzero_filter_prob)
    """

    def __link_inactive_quartile_count_scan(self, inactive_nonzero_filter_prob = 0.5):

        sp_id2mut_vector = {}
        for mut_id, sp_id, in self.link_vector:

            if sp_id not in sp_id2mut_vector: sp_id2mut_vector[sp_id] = [0] * self.sample_num
            for sample_id in self.mutation_status[mut_id]:
                sp_id2mut_vector[sp_id][sample_id] = 1

            
        for mut_id, sp_id, in self.link_vector:
            cur_splicing_count_vector = self.splicing_counts[sp_id]
            self.link_vector2inactive_quartile_count[(mut_id, sp_id)] = utils.quartile([int(cur_splicing_count_vector[j]) for j in range(self.sample_num) if sp_id2mut_vector[sp_id][j] == 0], inactive_nonzero_filter_prob)



    # simply cluster links by their topologies
    def __simple_cluster_link(self):

        simple_clustered_link = []
        mut_id2sp_id = {}
        for mut_id, sp_id in self.pruned_link_vector:
            if mut_id not in mut_id2sp_id: mut_id2sp_id[mut_id] = []
            mut_id2sp_id[mut_id].append(sp_id)

        clustered_link = []
        mut_ids = list(mut_id2sp_id)

        # first, append all the mut_ids to the remaining_ids
        remaining_ids = copy.deepcopy(mut_ids) # deep copy

        for i in range(len(mut_ids)):

            # start from the i-th mut_id if this is still remaining
            if mut_ids[i] not in remaining_ids: continue
            clustered_ids = [mut_ids[i]]

            no_more_cluster = 0
            while no_more_cluster == 0:

                no_more_cluster = 1
                for check_mut_id in remaining_ids:
                    if check_mut_id in clustered_ids: continue
                    for cur_clustered_mut_id in clustered_ids:
                        # check whether there is overlapping splicing aberrations
                        if len(set(mut_id2sp_id[cur_clustered_mut_id]) & set(mut_id2sp_id[check_mut_id])) > 0:
                            clustered_ids.append(check_mut_id)
                            no_more_cluster = 0
                            break


            tmp_clink_vector = []
            for c_mut_id in sorted(clustered_ids):
                tmp_clink_vector = tmp_clink_vector + [(c_mut_id, x) for x in mut_id2sp_id[c_mut_id]]
                remaining_ids.remove(c_mut_id)

            # simple_clustered_link.append(tmp_clink_vector)
            yield tmp_clink_vector

        # return simple_clustered_link


    def __get_log_marginal_likelihood(self, sub_cluster_link, configuration_vector, alpha0, beta0, alpha1, beta1):

        active_sp_list = []
        # get mutations associated with the current splicing
        for mut_id, sp_id in sub_cluster_link:
            active_sp_list.append(sp_id)


        log_marginal_likelihood = 0
        for sp_id in range(len(self.splicing_counts)):

            if sp_id not in active_sp_list: continue

            cur_splicing_count_vector = self.splicing_counts[sp_id]
            active_mut_list = []

            # get mutations associated with the current splicing
            for link_id in range(len(sub_cluster_link)):
                if configuration_vector[link_id] == 0: continue
                mut_id, link_sp_id = sub_cluster_link[link_id]
                if link_sp_id == sp_id: active_mut_list.append(mut_id)

            # param_num = param_num + len(active_mut_list) + 1
            active_mut_vector = [0] * self.sample_num

            # set mutation status
            for mut_id in self.mutation_status: 
                if mut_id in active_mut_list:
                    sample_id_vector = self.mutation_status[mut_id]
                    for sample_id in sample_id_vector:
                        active_mut_vector[sample_id] = 1

            # for inactive mutations
            sample_sum = len([n for n in range(self.sample_num) if active_mut_vector[n] == 0])
            count_sum = sum([int(cur_splicing_count_vector[n]) for n in range(self.sample_num) if active_mut_vector[n] == 0])
            weight_sum = sum([self.weight_vector[n] for n in range(self.sample_num) if active_mut_vector[n] == 0])

            partial_log_marginal_likelihood = 0
            if sample_sum > 0:
                partial_log_marginal_likelihood = partial_log_marginal_likelihood + math.lgamma(count_sum + alpha0) - math.lgamma(alpha0) + \
                                                    alpha0 * math.log(beta0) - (count_sum + alpha0) * math.log(weight_sum + beta0)

            for n in range(self.sample_num):
                if active_mut_vector[n] == 0: continue
                partial_log_marginal_likelihood = partial_log_marginal_likelihood + math.lgamma(int(cur_splicing_count_vector[n]) + alpha1) - math.lgamma(alpha1) + \
                                                    alpha1 * math.log(beta1) - (int(cur_splicing_count_vector[n]) + alpha1) * math.log(self.weight_vector[n] + beta1)

            log_marginal_likelihood = log_marginal_likelihood + partial_log_marginal_likelihood

        return(log_marginal_likelihood)


if __name__ == "__main__":

    from collections import namedtuple

    mutation_status = {0: [2]}
    splicing_counts = [[0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                       [0,0,19,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                       [0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]]
    # link_vector = [(0, 0), (0, 1), (0, 2)]
    # sample_num = 26
    sample_list = ["A427", "A549", "ABC-1", "H1299", "H1437", "H1648", "H1650", "H1703", "H1819", 
                   "H1975", "H2126", "H2228", "H2347", "H322", "II-18", "LC2_ad", "PC-14", "PC-3", 
                   "PC-7", "PC-9", "RERF-LC-Ad1", "RERF-LC-Ad2", "RERF-LC-KJ", "RERF-LC-MS", "RERF-LC-OK", "VMRC-LCD"]

    weight_vector = [5.2935, 2.5843, 4.6718, 6.7032, 6.3449, 4.8819, 3.2902, 11.6925, 9.7603, 4.4712, 5.8023, 7.3304, 6.6022,
                     6.8118, 8.0286, 5.5164, 6.7205, 6.3547, 6.5036, 4.2754, 6.989, 5.63, 7.514, 6.6156, 4.1002, 6.0268]

    Link_info_mut = namedtuple("Link_info_mut", ("Mutation_Key", "Motif_Pos", "Mutation_Type", "Is_Canonical", "Splicing_Key", "Splicing_Class", "Is_Inframe"))
    link2info = {(0, 0): Link_info_mut("12,475272,T,A", "12:475270-475276,-", "splicing acceptor disruption", "canonical", "12:461491-493196", "exon-skip", "in-frame"),
                 (0, 1): Link_info_mut("12,475272,T,A", "12:475270-475276,-", "splicing acceptor disruption", "canonical", "12:472264-493196", "exon-skip", "in-frame"),
                 (0, 2): Link_info_mut("12,475272,T,A", "12:475270-475276,-", "splicing acceptor disruption", "canonical", "12:475260-493196", "alternative-3'-splice-site", "in-frame")}


    network = Network("KDM5A", mutation_status, splicing_counts, link2info, sample_list, weight_vector)

    network.prune_link_vector(3.0, 0.5, 0.5)


    print(network.link_vector2effect_size)
    # print network.link_vector2median_count
    print(network.link_vector2inactive_quartile_count)
    print(network.pruned_link_vector)


    network.cluster_link_vector()
    print(network.clustered_link_vector)


    network.get_averaged_bayes_factors(1.0, 1.0, 1.0, 0.01)
    print(network.mut2log_BF)
    print(network.mut2significant_links)

    for sav in network.export_to_savs(3.0):
        print(sav.print_records(with_fdr = False))

 
