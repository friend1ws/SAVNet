#! /usr/bin/env python

from __future__ import print_function

import sys, gzip, subprocess, re, time 
import pysam

def merge_SJ2(SJ_file_list, output_file, control_file, junc_num_thres, is_keep_annotated):

    # get junctions registered in control
    control_db = {}
    if control_file is not None:
        with gzip.open(control_file, 'rt') as hin:
            for line in hin:
                F = line.rstrip('\n').split('\t')
                key = F[0] + '\t' + F[1] + '\t' + F[2]
                control_db[key] = 1


    # list up junctions to pick up
    junc2list = {}
    for SJ_file in SJ_file_list:
        with open(SJ_file, 'r') as hin:
            for line in hin:
                F = line.rstrip('\n').split('\t')
                if is_keep_annotated == False and F[5] != "0": continue
                if int(F[6]) < junc_num_thres: continue

                key = F[0] + '\t' + F[1] + '\t' + F[2]
                if key in control_db: continue
                if key not in junc2list: junc2list[key] = 1
                

    temp_id = 0
    hout = open(output_file + ".tmp.unsorted.txt", 'w')
    for SJ_file in SJ_file_list:
        with open(SJ_file, 'r') as hin:
            for line in hin:
                F = line.rstrip('\n').split('\t')
                if F[0] + '\t' + F[1] + '\t' + F[2] in junc2list:
                    print(F[0] + '\t' + F[1] + '\t' + F[2] + '\t' + str(temp_id) + '\t' + F[6], file = hout)
      
            temp_id = temp_id + 1 

    hout.close()


    hout = open(output_file + '.tmp.sorted.txt', 'w')
    subprocess.call(["sort", "-k1,1", "-k2,2n", "-k3,3n", output_file + ".tmp.unsorted.txt"], stdout = hout)
    hout.close()


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
                    print(temp_chr + '\t' + temp_start + '\t' + temp_end + '\t' + ','.join(temp_count), file = hout)

                temp_chr = F[0]
                temp_start = F[1]
                temp_end = F[2]
                temp_count = ["0"] * temp_id


            temp_count[int(F[3])] = F[4]

    # last check 

    if temp_chr != "":
        print(temp_chr + '\t' + temp_start + '\t' + temp_end + '\t' + ','.join(temp_count), file = hout)

    hout.close()
 
    # remove intermediate files
    subprocess.call(["rm", "-rf", output_file + ".tmp.unsorted.txt"])
    subprocess.call(["rm", "-rf", output_file + ".tmp.sorted.txt"])

    # if control_file is not None:
    # control_db.close()



def merge_intron_retention(IR_file_list, output_file, control_file, ratio_thres, num_thres):

    # get control intron retention info
    control_db = {}
    if control_file is not None:
        with gzip.open(control_file, 'rt') as hin:
            for line in hin:
                F = line.rstrip('\n').split('\t')
                control_db[F[0] + '\t' + F[1]] = 1


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

                # check the existence in control
                if F[0] + '\t' + F[1] in control_db: continue

                key = '\t'.join([F[header2ind[x]] for x in target_header])

                if key not in intron_retention2list: intron_retention2list[key] = 1


    temp_id = 0
    hout = open(output_file + ".tmp.unsorted.txt", 'w')
    for IR_file in IR_file_list:
        with open(IR_file, 'r') as hin:
            header = hin.readline().rstrip('\n').split('\t')

            for line in hin:
                F = line.rstrip('\n').split('\t')
                key = '\t'.join([F[header2ind[x]] for x in target_header])
                if key in intron_retention2list:
                    print(key + '\t' + str(temp_id) + '\t' + F[header2ind["Intron_Retention_Read_Count"]], file = hout)

        temp_id = temp_id + 1

    hout.close()


    hout = open(output_file + '.tmp.sorted.txt', 'w')
    subprocess.call(["sort", "-k1,1", "-k2,2n", output_file + ".tmp.unsorted.txt"], stdout = hout)
    hout.close()


    temp_chr = ""
    temp_pos = ""
    temp_key = ""
    temp_count = ["0"] * temp_id
    hout = open(output_file, 'w')
    print('\t'.join(target_header) + '\t' + "Read_Count_Vector", file = hout)

    with open(output_file + '.tmp.sorted.txt', 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')

            if F[1] != temp_pos or F[0] != temp_chr:

                # if not the first line 
                if temp_chr != "":
                    print(temp_key + '\t' + ','.join(temp_count), file = hout)

                temp_chr = F[0]
                temp_pos = F[1]
                temp_key = '\t'.join([F[header2ind[x]] for x in target_header])
                temp_count = ["0"] * temp_id

            temp_count[int(F[8])] = F[9]

    
    # last check 
    if temp_key != "":
        print(temp_key + '\t' + ','.join(temp_count) , file = hout)

    hout.close()
 
    # remove intermediate files
    subprocess.call(["rm", "-rf", output_file + ".tmp.unsorted.txt"])
    subprocess.call(["rm", "-rf", output_file + ".tmp.sorted.txt"])

    # if control_file is not None:
    # control_db.close()


def merge_chimera(chimera_file_list, output_file, control_file, num_thres, overhang_thres):

    # list up junctions to pick up
    chimera2list = {}
    header2ind = {}
    target_header = ["Chr_1", "Pos_1", "Dir_1", "Chr_2", "Pos_2", "Dir_2", "Inserted_Seq"]

    for chimera_file in chimera_file_list:
        with open(chimera_file, 'r') as hin:
            header = hin.readline().rstrip('\n').split('\t')
            for (i, cname) in enumerate(header):
                header2ind[cname] = i

            for line in hin:
                F = line.rstrip('\n').split('\t')
                if int(F[header2ind["Read_Pair_Num"]]) < num_thres: continue
                if int(F[header2ind["Max_Over_Hang_1"]]) < overhang_thres: continue
                if int(F[header2ind["Max_Over_Hang_2"]]) < overhang_thres: continue

                key = '\t'.join([F[header2ind[x]] for x in target_header])
                if key not in chimera2list: chimera2list[key] = 1


    temp_id = 0
    hout = open(output_file + ".tmp.unsorted.txt", 'w')
    for chimera_file in chimera_file_list:
        with open(chimera_file, 'r') as hin:
            for line in hin:
                F = line.rstrip('\n').split('\t')
                key = '\t'.join([F[header2ind[x]] for x in target_header])
                if key in chimera2list:
                    print(key + '\t' + str(temp_id) + '\t' + F[header2ind["Read_Pair_Num"]], file = hout)

        temp_id = temp_id + 1

    hout.close()


    hout = open(output_file + '.tmp.sorted.txt', 'w')
    subprocess.call(["sort", "-k1,1", "-k2,2n", "-k4,4", "-k5,5n", output_file + ".tmp.unsorted.txt"], stdout = hout)
    hout.close()

    if control_file is not None:
        control_db = pysam.TabixFile(control_file)
 
    temp_chr = ""
    temp_pos = ""
    temp_key = ""
    temp_count = ["0"] * temp_id
    hout = open(output_file, 'w')
    print('\t'.join(target_header) + '\t' + "Read_Count_Vector", file = hout)

    with open(output_file + '.tmp.sorted.txt', 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            key = '\t'.join([F[header2ind[x]] for x in target_header])
            if key != temp_key: 
 
                # if not the first line 
                if temp_key != "":
 
                    # skip if the junction is included in the control file
                    control_flag = 0
                    if control_file is not None:
                        tabixErrorFlag = 0
                        try:
                            records = control_db.fetch(temp_chr, int(temp_pos) - 5, int(temp_pos) + 5)
                        except Exception as inst:
                            # print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
                            # tabixErrorMsg = str(inst.args)
                            tabixErrorFlag = 1
 
                        if tabixErrorFlag == 0:
                            for record_line in records:
                                record = record_line.split('\t')
                                record_key = '\t'.join([record[header2ind[x]] for x in target_header])                            
                                if temp_key == record_key:
                                    control_flag = 1
 
                    if control_flag == 0:
                        print(temp_key + '\t' + ','.join(temp_count), file = hout)

                temp_chr = F[0]
                temp_pos = F[1] 
                temp_key = key
                temp_count = ["0"] * temp_id
 
            temp_count[int(F[7])] = F[8]

    # last check 
    # skip if the junction is included in the control file
    control_flag = 0
    if control_file is not None:
        tabixErrorFlag = 0
        try:
            records = control_db.fetch(temp_chr, int(temp_pos) - 5, int(temp_pos) + 5)
        except Exception as inst:
            # print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
            # tabixErrorMsg = str(inst.args)
            tabixErrorFlag = 1
 
        if tabixErrorFlag == 0:
            for record_line in records:
                record = record_line.split('\t')
                record_key = '\t'.join([record[header2ind[x]] for x in target_header])
                if temp_key == record_key:
                    control_flag = 1
 
    if control_flag == 0:
        print(temp_key + '\t' + ','.join(temp_count), file = hout)

    hout.close()
 
    # remove intermediate files
    subprocess.call(["rm", "-rf", output_file + ".tmp.unsorted.txt"])
    subprocess.call(["rm", "-rf", output_file + ".tmp.sorted.txt"])
 
    if control_file is not None:
        control_db.close()


def merge_mut(mutation_file_list, output_file):


    mut2sample = {}
    sample_ind = 0
    for mut_file in mutation_file_list:
        sample_ind = sample_ind + 1
        with open(mut_file, 'r') as hin2:
            for line2 in hin2:
                F2 = line2.rstrip('\n').split('\t')
                if F2[0].startswith('#'): continue
                if F2[0] == "Chr": continue

                key = '\t'.join(F2[0:5])
       
                if key not in mut2sample: 
                    mut2sample[key] = []
              
                mut2sample[key].append(str(sample_ind))

    sample_num = sample_ind

    hout = open(output_file, 'w')
    for mut in sorted(mut2sample):
        if len(mut2sample[mut]) == sample_num: continue
        print(mut + '\t' + ','.join(mut2sample[mut]), file = hout)

    hout.close()


def merge_mut2(mutation_file_list, output_file, reference):


    mut2sample = {}
    sample_ind = 0
    for mut_file in mutation_file_list:
        sample_ind = sample_ind + 1
        is_vcf = True if mut_file.endswith(".vcf") or mut_file.endswith(".vcf.gz") else False
        hin2 = gzip.open(mut_file, 'rt') if mut_file.endswith(".gz") else open(mut_file, 'r')

        for line2 in hin2:
            F2 = line2.rstrip('\n').split('\t')
            if F2[0].startswith('#'): continue
            if F2[0] == "Chr": continue

            if is_vcf == False:
                pos, ref, alt = F2[1], F2[3], F2[4]
        
                # insertion
                if F2[3] == "-":
                    # get the sequence for the reference base
                    seq = ""    
                    for item in pysam.faidx(reference, F2[0] + ":" + str(F2[1]) + "-" + str(F2[1])):
                        seq = seq + item.rstrip('\n')
                    seq = seq.replace('>', '')
                    seq = seq.replace(F2[0] + ":" + str(F2[1]) + "-" + str(F2[1]), '')
                    ref, alt = seq, seq + F2[4]

                # deletion
                if F2[4] == "-":
                    # get the sequence for the reference base
                    seq = ""    
                    for item in pysam.faidx(reference, F2[0] + ":" + str(int(F2[1]) - 1) + "-" + str(int(F2[1]) - 1)):
                        seq = seq + item.rstrip('\n')
                    seq = seq.replace('>', '')
                    seq = seq.replace(F2[0] + ":" + str(int(F2[1]) - 1) + "-" + str(int(F2[1]) - 1), '')
                    pos, ref, alt = str(int(F2[1]) - 1), seq + F2[3], seq

                QUAL = 60
                INFO = "SOMATIC"

                key = '\t'.join([F2[0], pos, '.', ref, alt, str(QUAL), "PASS", INFO])

            else:

                key = '\t'.join(F2[0:8])

            if key not in mut2sample:
                mut2sample[key] = []

            mut2sample[key].append(str(sample_ind))

        hin2.close()

    sample_num = sample_ind

    hout = open(output_file, 'w')
    for mut in sorted(mut2sample):
        if len(mut2sample[mut]) == sample_num: continue
        print(mut + '\t' + ','.join(mut2sample[mut]), file = hout)

    hout.close()


def merge_sv(sv_file_list, output_file):

    sv2sample = {}
    sample_num = "1"
    for sv_file in sv_file_list:
        with open(sv_file, 'r') as hin2:
            for line2 in hin2:
                F2 = line2.rstrip('\n').split('\t')
                if F2[0].startswith('#'): continue
                if F2[0] == "Chr_1": continue
                if F2[0] == "": continue

                key = '\t'.join(F2[0:7])

                if key not in sv2sample:
                    sv2sample[key] = []

                sv2sample[key].append(sample_num)

        sample_num = str(int(sample_num) + 1)

    hout = open(output_file, 'w')
    for sv in sorted(sv2sample):
        print(sv + '\t' + ','.join(sv2sample[sv]), file = hout)

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
            genes = sorted(list(set(genes)))

            if "---" in genes: genes.remove("---")
            if len(genes) > 0:
                genes_nm = list(filter(lambda x: x.find("(NM_") > 0, genes))
                if len(genes_nm) > 0: genes = genes_nm

            if len(genes) > 0:
                genes_single = list(filter(lambda x: x.find("-") == -1, genes))
                if len(genes_single) > 0: genes = genes_single
 
            gene = genes[0]
            gene = re.sub(r"\(N[MR]_\d+\)", "", gene)

            splicing_key = F[header2ind["SJ_1"]] + ':' + F[header2ind["SJ_2"]] + '-' + F[header2ind["SJ_3"]]

            print(gene + '\t' + splicing_key + '\t' + F[header2ind["Splicing_Class"]] + '\t' + F[header2ind["Is_Inframe"]] + '\t' + \
                  F[header2ind["SJ_4"]] + '\t' + F[header2ind["Mutation_Key"]] + '\t' + F[header2ind["Motif_Pos"]] + '\t' + \
                  F[header2ind["Mutation_Type"]] + '\t' + F[header2ind["Is_Canonical"]], file = hout)


    with open(IR_input_file, 'r') as hin:

        header = hin.readline().rstrip('\n').split('\t')
        for i in range(len(header)):
            header2ind[header[i]] = i

        for line in hin:
            F = line.rstrip('\n').split('\t')
            splicing_key = F[header2ind["Chr"]] + ':' + F[header2ind["Boundary_Pos"]] + '-' + F[header2ind["Boundary_Pos"]]
            splicing_class = "Intron retention" if F[header2ind["Intron_Retention_Type"]] == "Direct impact" else "Opposite side intron retention"
            print(F[header2ind["Gene_Symbol"]] + '\t' + splicing_key + '\t' + splicing_class + '\t' + '---' + '\t' + \
                  F[header2ind["Read_Count_Vector"]] + '\t' + F[header2ind["Mutation_Key"]] + '\t' + \
                  F[header2ind["Motif_Pos"]] + '\t' + F[header2ind["Mutation_Type"]] + '\t' + F[header2ind["Is_Canonical"]], file = hout)

    hout.close()

    hout = open(output_file, 'w')
    print('\t'.join(["Gene_Symbol", "Splicing_Key", "Splicing_Class", "Is_Inframe", "Read_Counts",
                     "Mutation_Key", "Motif_Pos", "Mutation_Type", "Is_Canonical"]), file = hout)
    hout.close()

    hout = open(output_file, 'a')
    subprocess.call(["sort", "-k1", output_file + ".unsorted"], stdout = hout)
    hout.close()

    subprocess.call(["rm", "-rf", output_file + ".unsorted"])


def get_sv_type(sv_key):

    sv_info = sv_key.split(',')
    sv_type = "dummy"
    if sv_info[0] != sv_info[3]:
        sv_type = "translocation"
    elif sv_info[2] == '+' and sv_info[5] == '-':
        sv_type = "deletion"
    elif sv_info[2] == '-' and sv_info[5] == '+':
        sv_type = "tandem_duplication"
    else:
        sv_type = "inversion"
    return(sv_type)


def merge_SJ_IR_chimera_files_sv(SJ_input_file, IR_input_file, chimera_input_file, output_file):

    header2ind = {}
    header = ""
    hout = open(output_file + ".unsorted", 'w')

    with open(SJ_input_file, 'r') as hin:

        header = hin.readline().rstrip('\n').split('\t')
        for i in range(len(header)):
            header2ind[header[i]] = i

        for line in hin:
            F = line.rstrip('\n').split('\t')
            gene1 = gene_filter(F[header2ind["Gene_1"]].split(';'))
            gene2 = gene_filter(F[header2ind["Gene_2"]].split(';'))
            gene = list(set(gene1) & set(gene2))
            if len(gene) == 0: continue

            splicing_key = F[header2ind["SJ_1"]] + ':' + F[header2ind["SJ_2"]] + '-' + F[header2ind["SJ_3"]]
            print(gene[0] + '\t' + splicing_key + '\t' + F[header2ind["Splicing_Class"]] + '\t' + F[header2ind["Is_Inframe"]] + '\t' + \
                  F[header2ind["SJ_4"]] + '\t' + F[header2ind["SV_Key"]] + '\t' + get_sv_type(F[header2ind["SV_Key"]]), file = hout)


    with open(IR_input_file, 'r') as hin:

        header = hin.readline().rstrip('\n').split('\t')
        for i in range(len(header)):
            header2ind[header[i]] = i

        for line in hin:
            F = line.rstrip('\n').split('\t')
            splicing_key = F[header2ind["Chr"]] + ':' + F[header2ind["Boundary_Pos"]] + '-' + F[header2ind["Boundary_Pos"]]
            splicing_class = "Intron retention"

            print(F[header2ind["Gene_Symbol"]] + '\t' + splicing_key + '\t' + splicing_class + '\t' + '---' + '\t' + \
                  F[header2ind["Read_Count_Vector"]] + '\t' + F[header2ind["SV_Key"]] + '\t' + get_sv_type(F[header2ind["SV_Key"]]), file = hout)

    with open(chimera_input_file, 'r') as hin: 

        header = hin.readline().rstrip('\n').split('\t')
        for i in range(len(header)):
            header2ind[header[i]] = i


        for line in hin: 
            F = line.rstrip('\n').split('\t')

            gene1 = gene_filter(F[header2ind["Gene_1"]].split(';'))
            gene2 = gene_filter(F[header2ind["Gene_2"]].split(';'))
            gene = list(set(gene1) & set(gene2))
            if len(gene) == 0: continue

            # currently, only consider exon_reusage and unspliced_chimera
            if F[header2ind["Chimera_Class"]] not in ["Exon reusage", "Unspliced chimera"]: continue

            splicing_key = ','.join([F[header2ind[x]] for x in ["Chr_1", "Pos_1", "Dir_1", "Chr_2", "Pos_2", "Dir_2", "Inserted_Seq"]])
            print(gene[0] + '\t' + splicing_key + '\t' + F[header2ind["Chimera_Class"]] + '\t' + "---" + '\t' + \
                  F[header2ind["Read_Count_Vector"]] + '\t' + F[header2ind["SV_Key"]] + '\t' + get_sv_type(F[header2ind["SV_Key"]]), file = hout)

    hout.close()

    hout = open(output_file, 'w')
    print('\t'.join(["Gene_Symbol", "Splicing_Key", "Splicing_Class", "Is_Inframe", "Read_Counts", "SV_Key", "SV_Type"]), file = hout)
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

        print("Gene_Symbol" + '\t' + '\t'.join(header), file = hout0)

        for line in hin:
            F = line.rstrip('\n').split('\t')
            genes = F[header2ind["Gene_1"]].split(';') + F[header2ind["Gene_2"]].split(';')
            genes = list(set(genes))

            if "---" in genes: genes.remove("---")
            if len(genes) > 0: 
                genes_nm = list(filter(lambda x: x.find("(NM_") > 0, genes))
                if len(genes_nm) > 0: genes = genes_nm

            gene = genes[0]
            gene = re.sub(r"\(N[MR]_\d+\)", "", gene)

            print(gene + '\t' + '\t'.join(F), file = hout1)

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



def gene_filter(genes):

    genes = list(set(genes))
    if "---" in genes: genes.remove("---")
    if len(genes) > 0: 
        genes_nm = list(filter(lambda x: x.find("(NM_") > 0, genes))
        if len(genes_nm) > 0: genes = genes_nm
   
    if len(genes) > 0: 
        genes = map(lambda x: re.sub(r"\(N[MR]_\d+\)", "", x), genes)
    else:
        genes = [] 

    return(genes)


