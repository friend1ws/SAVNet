#! /usr/bin/env python

from __future__ import print_function

import os, glob

def get_weight(log_final_file):

    with open(log_final_file) as hin:
        for line in hin:
            if "Uniquely mapped reads number |" in line:
                read_num = int(line.replace("Uniquely mapped reads number |", '').strip(' '))
                weight = float(read_num) / float(10000000)
                break

    return(weight)


def make_savnet_input(output_file, mut_dir, sj_dir, ir_dir, qc_dir):
 
    mut_files = glob.glob(mut_dir + "/*.avinput") + glob.glob(mut_dir + "/*.vcf")
    sj_files = glob.glob(sj_dir + "/*.SJ.out.tab") 
    ir_files = glob.glob(ir_dir + "/*.intron_retention.txt")
    qc_files = glob.glob(qc_dir + "/*.Log.final.out")


    sample2mut_file = {}
    for mut_file in sorted(mut_files):
        sample = os.path.basename(mut_file)
        sample = sample.replace(".avinput", "")
        sample = sample.replace(".vcf", "")
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
        sample = os.path.basename(ir_file).replace(".intron_retention.txt", '')
        sample = sample[:15]
        sample2ir_file[sample] = os.path.abspath(ir_file)


    with open(output_file, 'w') as hout:
        print('\t'.join(["Sample_Name", "Weight", "Mutation_File", "SJ_File", "IR_File"]), file = hout)
        for sample in sorted(sample2mut_file):
            if sample not in sample2sj_file: continue
            if sample not in sample2ir_file: continue
            print(sample + '\t' + str(round(sample2weight[sample], 4)) + '\t' + sample2mut_file[sample] + '\t' + sample2sj_file[sample] + '\t' + sample2ir_file[sample], file = hout)


def make_savnet_input_sv(output_file, sv_dir, sj_dir, ir_dir, chimera_dir, qc_dir):

    sv_files = glob.glob(sv_dir + "/*.genomonSV.result.filt.txt")
    sj_files = glob.glob(sj_dir + "/*.SJ.out.tab")
    ir_files = glob.glob(ir_dir + "/*.intron_retention.txt")
    chimera_files = glob.glob(chimera_dir + "/*.chimera.count.txt")
    qc_files = glob.glob(qc_dir + "/*.Log.final.out")

    sample2sv_file = {}
    for sv_file in sorted(sv_files):
        sample = os.path.basename(sv_file).replace(".genomonSV.result.filt.txt", '')
        sample = sample[:15]
        sample2sv_file[sample] = sv_file


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
        sample = os.path.basename(ir_file).replace(".intron_retention.txt", '')
        sample = sample[:15]
        sample2ir_file[sample] = os.path.abspath(ir_file)


    sample2chimera_file = {}
    for chimera_file in sorted(chimera_files):
        sample = os.path.basename(chimera_file).replace(".chimera.count.txt", '')
        sample = sample[:15]
        sample2chimera_file[sample] = chimera_file


    hout = open(output_file, 'w')
    print('\t'.join(["Sample_Name", "Weight", "SV_File", "SJ_File", "IR_File", "Chimera_File"]), file = hout)
    for sample in sorted(sample2sv_file):
        if sample not in sample2sj_file: continue
        if sample not in sample2ir_file: continue
        if sample not in sample2chimera_file: continue
        print(sample + '\t' + str(round(sample2weight[sample], 4)) + '\t' + \
              sample2sv_file[sample] + '\t' + sample2sj_file[sample] + '\t' + \
              sample2ir_file[sample] + '\t' + sample2chimera_file[sample], file = hout) 


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser("make_savnet_input")

    parser.add_argument("--sample_list_file", default = None, required = True, type = str)

    parser.add_argument("--sv", default = False, action = 'store_true',
                        help = "Make SAVNet input for structural variation mode")

    parser.add_argument("--mut_dir", default = None, required = True, type = str,
                        help = "The directory path to mutation or sv files")

    parser.add_argument("--sj_dir", default = None, type = str,
                        help = "The directory path to splicing junction files")
    
    parser.add_argument("--ir_dir", default = None, type = str,
                        help = "The directory path to intron retention files")

    parser.add_argument("--chimera_dir", default = None, type = str,
                        help = "The directory path to chimeric count files")

    parser.add_argument("--qc_dir", default = None, type = str,
                        help = "The directory path to quality control files")

    args = parser.parse_args()

    if args.sv:
        make_savnet_input_sv(args.sample_list_file, args.mut_dir, 
                             args.sj_dir, args.ir_dir, args.chimera_dir, args.qc_dir)
    else:
        make_savnet_input(args.sample_list_file, args.mut_dir, 
                          args.sj_dir, args.ir_dir, args.qc_dir)


