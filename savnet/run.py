#! /usr/bin/env python

from __future__ import print_function
import sys, subprocess, os

from . import preprocess, analysis_network, sample_conf
from .sav import Sav
from .utils import is_tool

from .logger import get_logger
logger = get_logger(__name__)

def savnet_main(args):

    if args.grc == True:
        logger.warning("--grc argument is deprecated and ignored.")


    ##########
    # read sample conf
    sconf = sample_conf.Sample_conf()
    sconf.parse_file(args.sample_list_file, args.sv)

    ##########
    # check if the executables exist
    is_tool("bedtools")
    is_tool("tabix")
    is_tool("bgzip")
    if len(sconf.SJ_files) > 0: is_tool("junc_utils")
    if len(sconf.IR_files) > 0: is_tool("intron_retention_utils")
    if len(sconf.chimera_files) > 0: is_tool("chimera_utils")
    
    ##########
    output_prefix_dir = os.path.dirname(args.output_prefix)
    if output_prefix_dir != "" and not os.path.exists(output_prefix_dir):
       os.makedirs(output_prefix_dir)

    ##########
    logger.info("Merging mutation data.")
    if args.sv == False:
        preprocess.merge_mut2(sconf.mut_files, args.output_prefix + ".mut_merged.vcf", args.reference)
    else:
        preprocess.merge_sv(sconf.sv_files, args.output_prefix + ".sv_merged.txt")

    ##########
    # splicing_junction
    logger.info("Merging splicing junction data.")
    preprocess.merge_SJ2(sconf.SJ_files, args.output_prefix + ".SJ_merged.txt", args.SJ_pooled_control_file, args.SJ_num_thres, args.keep_annotated)

    logger.info("Adding annotation to splicing junction data.")
    annotate_commands = ["junc_utils", "annotate", args.output_prefix + ".SJ_merged.txt", args.output_prefix + ".SJ_merged.annot.txt",
                         "--genome_id", args.genome_id]
    # if args.grc: annotate_commands.append("--grc")
    subprocess.call(annotate_commands)

    logger.info("Checking association betweeen mutation and splicing junction data.")
    if args.sv == False:
        associate_commands = ["junc_utils", "associate", args.output_prefix + ".SJ_merged.annot.txt", args.output_prefix + ".mut_merged.vcf",
                              args.output_prefix + ".SJ_merged.associate.txt",
                              "--donor_size", args.donor_size, "--acceptor_size", args.acceptor_size,
                              "--genome_id", args.genome_id]
        # if args.branchpoint: associate_commands.append("--branchpoint")
        # if args.grc: associate_commands.append("--grc")

    else:
        associate_commands = ["junc_utils", "associate", args.output_prefix + ".SJ_merged.annot.txt", args.output_prefix + ".sv_merged.txt",
                              args.output_prefix + ".SJ_merged.associate.txt", "--sv"]

    subprocess.check_call(associate_commands)
    ##########

    ##########
    # intron_retention
    logger.info("Merging intron retention data.")
    preprocess.merge_intron_retention(sconf.IR_files, args.output_prefix + ".IR_merged.txt", 
                                 args.IR_pooled_control_file, args.IR_ratio_thres, args.IR_num_thres)

    logger.info("Checking association betweeen mutation and intron retention data.")
    if args.sv == False:
        associate_commands = ["intron_retention_utils", "associate", args.output_prefix + ".IR_merged.txt",
                              args.output_prefix + ".mut_merged.vcf", args.output_prefix + ".IR_merged.associate.txt",
                              "--donor_size", args.donor_size, "--acceptor_size", args.acceptor_size]
    else:
        associate_commands = ["intron_retention_utils", "associate", args.output_prefix + ".IR_merged.txt",
                              args.output_prefix + ".sv_merged.txt", args.output_prefix + ".IR_merged.associate.txt", "--sv"]

    subprocess.check_call(associate_commands)
    #########

    #########
    # chimera
    if args.sv:
        logger.info("Merging chimeric junction data.")
        preprocess.merge_chimera(sconf.chimera_files, args.output_prefix + ".chimera_merged.txt", 
                            args.chimera_pooled_control_file, args.chimera_num_thres, args.chimera_overhang_thres)

        logger.info("Checking association betweeen mutation and chimeric junction data.")
        associate_commands = ["chimera_utils", "associate", args.output_prefix + ".chimera_merged.txt",
                              args.output_prefix + ".sv_merged.txt", args.output_prefix + ".chimera_merged.associate.txt", "--genome_id", args.genome_id]
        # if args.grc: associate_commands.append("--grc")

        subprocess.check_call(associate_commands)
    ##########

    ##########
    # organize association
    if args.sv == False:
        logger.info("Organizing splicing association information.")
        preprocess.merge_SJ_IR_files(args.output_prefix + ".SJ_merged.associate.txt", 
                                     args.output_prefix + ".IR_merged.associate.txt",
                                     args.output_prefix + ".splicing.associate.txt")

        logger.info("Creating pickles of splicing association network instances.")
        analysis_network.create_network_list(args.output_prefix + ".splicing.associate.txt",
                                             args.output_prefix + ".splicing_mutation.network.pickles",
                                             args.output_prefix + ".mut_merged.vcf",
                                             sconf.sample_names, sconf.weights)

    else:
        logger.info("Organizing splicing association information.")
        preprocess.merge_SJ_IR_chimera_files_sv(args.output_prefix + ".SJ_merged.associate.txt",
                                                args.output_prefix + ".IR_merged.associate.txt",
                                                args.output_prefix + ".chimera_merged.associate.txt",
                                                args.output_prefix + ".splicing.associate.txt")

        logger.info("Creating pickles of splicing association network instances.")
        analysis_network.create_network_list(args.output_prefix + ".splicing.associate.txt", 
                                             args.output_prefix + ".splicing_mutation.network.pickles",
                                             args.output_prefix + ".sv_merged.txt",
                                             sconf.sample_names, sconf.weights, sv_mode = True)


    logger.info("Extracting splicing associated variants.")
    sav_list_target = analysis_network.extract_sav_list(args.output_prefix + ".splicing_mutation.network.pickles", 
                                                        args.effect_size_thres, 0.5, 0.5, args.log_BF_thres, 1, 
                                                        args.alpha0, args.beta0, args.alpha1, args.beta1, permutation = False)

    logger.info("Extracting of splicing associated variants on permutation pairs to estimate false positive ratios.")
    sav_lists_permutation = []
    for i in range(args.permutation_num):
        temp_sav_list = analysis_network.extract_sav_list(args.output_prefix + ".splicing_mutation.network.pickles", 
                                                          args.effect_size_thres, 0.5, 0.5, args.log_BF_thres, 1,
                                                          args.alpha0, args.beta0, args.alpha1, args.beta1, permutation = True)
        sav_lists_permutation.append(temp_sav_list)

    logger.info("Adding Q-values to splicing associated variants.")
    analysis_network.add_qvalue_to_sav_list(sav_list_target, sav_lists_permutation)

    logger.info("Generating final outputs.")
    with open(args.output_prefix + ".savnet.result.txt", 'w') as hout:
        if args.sv == False:
            print(Sav.print_header_mut, file = hout)
        else:
            print(Sav.print_header_sv, file = hout)
        for sav in sav_list_target:
            print('\n'.join(sav.print_records(sv_mode = args.sv, with_fdr = True)), file = hout)

    with open(args.output_prefix + ".splicing_mutation.count_summary.anno.perm_all.txt", 'w') as hout:
        if args.sv == False:
            print("Permutation_Num" + '\t' + Sav.print_header_mut, file = hout)
        else:
            print("Permutation_Num" + '\t' + Sav.print_header_sv, file = hout)

        for i in range(len(sav_lists_permutation)):
            for sav in sav_lists_permutation[i]:
                print('\n'.join([str(i) + '\t' + x for x in sav.print_records(sv_mode = args.sv, with_fdr = False)]), file = hout)

    if args.debug == False:

        subprocess.call(["rm", "-rf", args.output_prefix + ".mut_merged.vcf"])
        subprocess.call(["rm", "-rf", args.output_prefix + ".sv_merged.txt"])
        subprocess.call(["rm", "-rf", args.output_prefix + ".SJ_merged.txt"])
        subprocess.call(["rm", "-rf", args.output_prefix + ".SJ_merged.annot.txt"])
        subprocess.call(["rm", "-rf", args.output_prefix + ".SJ_merged.associate.txt"])
        subprocess.call(["rm", "-rf", args.output_prefix + ".IR_merged.txt"])
        subprocess.call(["rm", "-rf", args.output_prefix + ".IR_merged.associate.txt"])
        subprocess.call(["rm", "-rf", args.output_prefix + ".chimera_merged.txt"])
        subprocess.call(["rm", "-rf", args.output_prefix + ".chimera_merged.associate.txt"])
        subprocess.call(["rm", "-rf", args.output_prefix + ".splicing.associate.txt"])
        subprocess.call(["rm", "-rf", args.output_prefix + ".splicing_mutation.network.pickles"]) 

