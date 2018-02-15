#! /usr/bin/env python

class Sav(object):

    def __init__(self, gene, sample_list, link_info, supporting_read_num_vector, score):
        self.gene = gene
        self.sample_list = sample_list
        self.link_info = link_info
        self.supporting_read_num_vector = supporting_read_num_vector
        self.score = score
        self.fdr = None

    def set_fdr(self, fdr):
        self.fdr = fdr


    def print_records(self):
        # return '\t'.join([self.gene, ';'.join(self.sample_list), self.link_info])
        return self.gene + '\t' + ';'.join(self.sample_list)

