#! /usr/bin/env python

import os, glob

def make_savnet_input(output_file, mut_dir, sj_dir, ir_dir, qc_dir):
 
    mut_files = glob.glob(mut_dir + "/*.mutation.txt")
    sj_files = glob.glob(sj_dir + "/*.SJ.out.tab") 
    ir_files = glob.glob(ir_dir + "/*.genomonIR.result.txt")
    qc_files = glob.glob(qc_dir + "/*.Log.final.out")


    def get_weight(log_final_file):

        with open(log_final_file) as hin:
            for line in hin:
                if "Uniquely mapped reads number |" in line:
                    read_num = int(line.replace("Uniquely mapped reads number |", '').strip(' '))
                    weight = float(read_num) / float(10000000)
                    break
        
        return(weight) 


    sample2mut_file = {}
    for mut_file in sorted(mut_files):
        sample = os.path.basename(mut_file).replace(".mutation.txt", "")
        sample = sample[:15]
        sample2mut_file[sample] = os.path.abspath(mut_file)

    sample2sj_file = {}
    for sj_file in sorted(sj_files):
        sample = os.path.basename(sj_file).replace(".SJ.out.tab", "")
        sample = sample[:15]
        sample2sj_file[sample] = os.path.abspath(sj_file)


    sample2weight = {}
    for qc_file in sorted(qc_files):
        sample = os.path.basename(qc_file).replace(".Log.final.out", "")
        sample = sample[:15]
        weight = get_weight(qc_file)
        sample2weight[sample] = weight


    sample2ir_file = {}
    for ir_file in sorted(ir_files):
        sample = os.path.basename(ir_file).replace(".genomonIR.result.txt", '')
        sample = sample[:15]
        sample2ir_file[sample] = os.path.abspath(ir_file)


    with open(output_file, 'w') as hout:
        print >> hout, '\t'.join(["Sample_Name", "Weight", "Mutation_File", "SJ_File", "IR_File"])
        for sample in sorted(sample2mut_file):
            if sample not in sample2sj_file: continue
            if sample not in sample2ir_file: continue
            print >> hout, sample + '\t' + str(round(sample2weight[sample], 4)) + '\t' + sample2mut_file[sample] + '\t' + sample2sj_file[sample] + '\t' + sample2ir_file[sample]



if __name__ == "__main__":

    import sys

    output_file = sys.argv[1]
    mut_dir = sys.argv[2]
    sj_dir = sys.argv[3]
    ir_dir = sys.argv[4]
    qc_dir = sys.argv[5]

    make_savnet_input(output_file, mut_dir, sj_dir, ir_dir, qc_dir)

