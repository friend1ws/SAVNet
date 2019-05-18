#! /usr/bin/env python

from __future__ import print_function

class Sav(object):

    print_header_mut = '\t'.join(["Gene_Symbol", "Sample_Name", "Mutation_Key", "Motif_Pos", "Mutation_Type", "Is_Canonical",
                                  "Splicing_Key", "Splicing_Class", "Is_Inframe", "Supporting_Read_Num", "Score", "Q_Value"])

    print_header_sv = '\t'.join(["Gene_Symbol", "Sample_Name", "SV_Key", "SV_Type",
                                  "Splicing_Key", "Splicing_Class", "Is_Inframe", "Supporting_Read_Num", "Score", "Q_Value"])
 
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

    def print_records(self, sv_mode = False, with_fdr = True):

        records = []
        for i in range(self.link_num):

            link_info = self.link_info_vector[i]
            supporting_read_num = self.supporting_read_num_vector[i]

            if sv_mode == False:
                print_link_info = '\t'.join([link_info.Mutation_Key, link_info.Motif_Pos, link_info.Mutation_Type,
                                             link_info.Is_Canonical, link_info.Splicing_Key, link_info.Splicing_Class, link_info.Is_Inframe])
            else:
                print_link_info = '\t'.join([link_info.SV_Key, link_info.SV_Type, link_info.Splicing_Key, link_info.Splicing_Class, link_info.Is_Inframe])

            temp_record = self.gene + '\t' + ';'.join(self.sample_list) + '\t' + print_link_info + '\t' + \
                            ';'.join([str(x) for x in supporting_read_num]) + '\t' + str(round(self.score, 4))

            if with_fdr:
                temp_record = temp_record + '\t' + str(round(self.fdr, 4))
            else:
                temp_record = temp_record + '\t' + "---"

            records.append(temp_record)

        return records

