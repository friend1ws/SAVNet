#! /usr/bin/env python

import copy


class Network(object):

    def __init__(self, gene, mutation_status, splicing_counts, link_vector, sample_num, weight_vector):
        self.gene = gene
        self.mutation_status = mutation_status
        self.splicing_counts = splicing_counts
        self.link_vector = link_vector
        self.sample_num = sample_num
        self.weight_vector = weight_vector
        self.pruned_link_vector = []

        self.link_vector2effect_size = {}
        self.link_vector2median_count = {}
        self.clustered_link_vector = []



    def prune_link_vector(self, margin):

        self.__link_effect_size_scan()
        self.__link_median_count_scan()
        self.pruned_link_vector = [x for x in self.link_vector if self.link_vector2effect_size[x] >= margin and self.link_vector2median_count[x] == 0]
        

    def cluster_link_vector(self):

        # self.__simple_cluster_link()
        # simple_clusterd_link_vector = self.clustered_link_vector
        # for sub_cluster_link in simple_clusterd_link_vector:
        for sub_cluster_link in self.__simple_cluster_link():

            if len(sub_cluster_link) > 15:

                sub_cluster_effect_size_vector = [self.link_vector2effect_size[x] for x in sub_cluster_link]
                temp_margin = sorted(sub_cluster_effect_size_vector, reverse=True)[15]
                sub_cluster_link = [sub_cluster_link[i] for i in range(len(sub_cluster_link)) if sub_cluster_effect_size_vector[i] > temp_margin]

            # may be removed at the above filtering step
            if len(sub_cluster_link) > 0:
                self.clustered_link_vector.append(sub_cluster_link)


    def __link_effect_size_scan(self, pseudo_count = 0.1):

        # simple check for each link
        for mut_id, sp_id in self.link_vector:
            cur_splicing_count_vector = self.splicing_counts[sp_id] 

            # extract samples with the mutation of the link in consideration
            mut_vector = [0] * self.sample_num
            for sample_id in self.mutation_status[mut_id]:
                mut_vector[sample_id] = 1

            weight_sum_null = sum([float(weight_vector[j]) for j in range(sample_num) if mut_vector[j] == 0])
            weight_sum_target = sum([float(weight_vector[j]) for j in range(sample_num) if mut_vector[j] == 1])

            count_sum_null = sum([int(cur_splicing_count_vector[j]) for j in range(sample_num) if mut_vector[j] == 0])
            count_sum_target = sum([int(cur_splicing_count_vector[j]) for j in range(sample_num) if mut_vector[j] == 1])

            mean_null = float(count_sum_null) / weight_sum_null if weight_sum_null > 0.0 and count_sum_null > 0 else 0.0
            mean_target = float(count_sum_target) / weight_sum_target if weight_sum_target > 0 and count_sum_target > 0 else 0.0
    
            self.link_vector2effect_size[(mut_id, sp_id)] = (float(mean_target) + pseudo_count) / (float(mean_null) + pseudo_count)



    def __link_median_count_scan(self):

        def median(numbers):
            if len(numbers) == 0:
                print >> sys.stderr, "Vector of zero length was put to median function. Return 0"
                return 0
            return (sorted(numbers)[int(round((len(numbers) - 1) / 2.0))] + sorted(numbers)[int(round((len(numbers) - 1) // 2.0))]) / 2.0

        # simple check for each link
        for mut_id, sp_id in link_vector:
            # cur_mut_id, cur_sp_id = self.link_vector[i]
            cur_splicing_count_vector = self.splicing_counts[sp_id]

            # extract samples with the mutation of the link in consideration
            mut_vector = [0] * self.sample_num
            mut_vector = [0] * self.sample_num
            for sample_id in self.mutation_status[mut_id]:
                mut_vector[sample_id] = 1

            self.link_vector2median_count[(mut_id, sp_id)] = median([int(cur_splicing_count_vector[j]) for j in range(sample_num) if mut_vector[j] == 0])
    

    # simply cluster links by their topologies
    def __simple_cluster_link(self):

        simple_clustered_link = []
        mut_id2sp_id = {}
        for mut_id, sp_id in self.pruned_link_vector:
            if mut_id not in mut_id2sp_id: mut_id2sp_id[mut_id] = []
            mut_id2sp_id[mut_id].append(sp_id)

        clustered_link = []
        mut_ids = mut_id2sp_id.keys()

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
                    for cur_clustered_mut_id in clusteted_ids:
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

if __name__ == "__main__":

    mutation_status = {0: [2]}
    splicing_counts = [[0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                       [0,0,19,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                       [0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]]
    link_vector = [(0, 0), (0, 1), (0, 2)]
    sample_num = 26
    weight_vector = [1] * 26

    network = Network("KDM5A", mutation_status, splicing_counts, link_vector, sample_num, weight_vector)
    network.prune_link_vector(3.0)


    print network.link_vector2effect_size
    print network.link_vector2median_count
    print network.pruned_link_vector


    network.cluster_link_vector()
    print network.clustered_link_vector



