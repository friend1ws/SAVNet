#! /usr/bin/env python

class Sav(object):

    def __init__(self, gene, sample_list, link_info_vector, supporting_read_num_vector, score):
        self.gene = gene
        self.sample_list = sample_list
        self.link_info_vector = link_info_vector
        self.link_num = len(link_info_vector)
        self.supporting_read_num_vector = supporting_read_num_vector
        self.score = score
        self.fdr = None

    def set_fdr(self, fdr):
        self.fdr = fdr


    def print_records(self):

        records = []
        for i in range(self.link_num):

            link_info = self.link_info_vector[i]
            supporting_read_num = self.supporting_read_num_vector[i]

            print_link_info = '\t'.join([link_info.Mutation_Key, link_info.Motif_Pos, link_info.Mutation_Type,
                                         link_info.Is_Canonical, link_info.Splicing_Key, link_info.Splicing_Class, link_info.Is_Inframe])

            records.append(self.gene + '\t' + ';'.join(self.sample_list) + '\t' + print_link_info + '\t' + \
                     ';'.join([str(x) for x in supporting_read_num]) + '\t' + str(round(float(self.score), 4)))

        return '\n'.join(records) 

