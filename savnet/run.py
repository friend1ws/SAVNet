#! /usr/bin/env python

import sys, subprocess, os
import preprocess, analysis_network, sample_conf
from sav import Sav

def savnet_main(args):

    output_prefix_dir = os.path.dirname(args.output_prefix)
    if output_prefix_dir != "" and not os.path.exists(output_prefix_dir):
       os.makedirs(output_prefix_dir)

    ##########
    # read sample conf
    sconf = sample_conf.Sample_conf()
    sconf.parse_file(args.sample_list_file, args.sv)

    """

    if args.sv == False:
        preprocess.merge_mut(sconf.mut_files, args.output_prefix + ".mut_merged.txt")
    else:
        preprocess.merge_sv(sconf.sv_files, args.output_prefix + ".sv_merged.txt")

    ##########
    # splicing_junction
    preprocess.merge_SJ2(sconf.SJ_files, args.output_prefix + ".SJ_merged.txt", args.SJ_pooled_control_file, args.SJ_num_thres, args.keep_annotated)

    annotate_commands = ["junc_utils", "annotate", args.output_prefix + ".SJ_merged.txt", args.output_prefix + ".SJ_merged.annot.txt",
                         "--genome_id", args.genome_id]
    if args.grc: annotate_commands.append("--grc")
    subprocess.call(annotate_commands)

    if args.sv == False:
        associate_commands = ["junc_utils", "associate", args.output_prefix + ".SJ_merged.annot.txt", args.output_prefix + ".mut_merged.txt",
                              args.output_prefix + ".SJ_merged.associate.txt", "--reference", args.reference_genome,
                              "--mutation_format", "anno", "--donor_size", args.donor_size, "--acceptor_size", args.acceptor_size,
                              "--genome_id", args.genome_id]
        # if args.branchpoint: associate_commands.append("--branchpoint")
        if args.grc: associate_commands.append("--grc")

    else:
        associate_commands = ["junc_utils", "associate", args.output_prefix + ".SJ_merged.annot.txt", args.output_prefix + ".sv_merged.txt",
                              args.output_prefix + ".SJ_merged.associate.txt", "--sv"]

    subprocess.check_call(associate_commands)
    ##########

    ##########
    # intron_retention
    preprocess.merge_intron_retention(sconf.IR_files, args.output_prefix + ".IR_merged.txt", 
                                 args.IR_pooled_control_file, args.IR_ratio_thres, args.IR_num_thres)

    if args.sv == False:
        associate_commands = ["intron_retention_utils", "associate", args.output_prefix + ".IR_merged.txt",
                              args.output_prefix + ".mut_merged.txt", args.output_prefix + ".IR_merged.associate.txt",
                              "--reference", args.reference_genome, "--mutation", "anno",
                              "--donor_size", args.donor_size, "--acceptor_size", args.acceptor_size]
    else:
        associate_commands = ["intron_retention_utils", "associate", args.output_prefix + ".IR_merged.txt",
                              args.output_prefix + ".sv_merged.txt", args.output_prefix + ".IR_merged.associate.txt", "--sv"]

    subprocess.check_call(associate_commands)
    #########

    #########
    # chimera
    if args.sv:
        preprocess.merge_chimera(sconf.chimera_files, args.output_prefix + ".chimera_merged.txt", 
                            args.chimera_pooled_control_file, args.chimera_num_thres, args.chimera_overhang_thres)

        subprocess.call(["chimera_utils", "associate", args.output_prefix + ".chimera_merged.txt",
                         args.output_prefix + ".sv_merged.txt", args.output_prefix + ".chimera_merged.associate.txt"])
    ##########

    ##########
    # organize association
    if args.sv == False:
        preprocess.merge_SJ_IR_files(args.output_prefix + ".SJ_merged.associate.txt", 
                                     args.output_prefix + ".IR_merged.associate.txt",
                                     args.output_prefix + ".splicing.associate.txt")

    else:
        preprocess.merge_SJ_IR_chimera_files_sv(args.output_prefix + ".SJ_merged.associate.txt",
                                                args.output_prefix + ".IR_merged.associate.txt",
                                                args.output_prefix + ".chimera_merged.associate.txt",
                                                args.output_prefix + ".splicing.associate.txt")

    """

    analysis_network.create_network_list(args.output_prefix + ".splicing.associate.txt", 
                                         args.output_prefix + ".splicing_mutatoin.network.pickles",
                                         args.output_prefix + ".mut_merged.txt",
                                         sconf.sample_names, sconf.weights)


    sav_list_target = analysis_network.extract_sav_list(args.output_prefix + ".splicing_mutatoin.network.pickles", args.effect_size_thres, args.log_BF_thres, 
                                                        args.alpha0, args.beta0, args.alpha1, args.beta1, permutation = False)

    sav_lists_permutation = []
    for i in range(args.permutation_num):
        temp_sav_list = analysis_network.extract_sav_list(args.output_prefix + ".splicing_mutatoin.network.pickles", args.effect_size_thres, args.log_BF_thres,
                                                          args.alpha0, args.beta0, args.alpha1, args.beta1, permutation = True)
        sav_lists_permutation.append(temp_sav_list)
 
    analysis_network.add_qvalue_to_sav_list(sav_list_target, sav_lists_permutation)

    with open(args.output_prefix + ".savnet.result.txt", 'w') as hout:
        print >> hout, Sav.print_header_mut
        for sav in sav_list_target:
            print >> hout, '\n'.join(sav.print_records(with_fdr = True))

    with open(args.output_prefix + ".splicing_mutation.count_summary.anno.perm_all.txt", 'w') as hout:
        print >> hout, "Permutation_Num" + '\t' + Sav.print_header_mut
        for i in range(len(sav_lists_permutation)):
            for sav in sav_lists_permutation[i]:
                print >> hout, '\n'.join([str(i) + '\t' + x for x in sav.print_records(with_fdr = False)])
 

    """
        utils.organize_mut_splicing_count(args.output_prefix + ".splicing.associate.txt",
                                    args.output_prefix + ".sv_merged.txt",
                                    args.output_prefix + ".splicing_mutation.count_summary.txt",
                                    args.output_prefix + ".splicing_mutation.mut_info.txt",
                                    args.output_prefix + ".splicing_mutation.splicing_info.txt", True)

    utils.organize_mut_splicing_count2(args.output_prefix + ".splicing.associate.txt",
                                       args.output_prefix + ".mut_merged.txt",
                                       args.output_prefix + ".splicing_mutation.count_summary.txt",
                                       args.output_prefix + ".splicing_mutation.link_info.txt")
    """


    """
    # true combination    
    print >> sys.stderr, "evaluating true combinations"

    utils.convert_pruned_file(args.output_prefix + ".splicing_mutation.count_summary.txt",
                              args.output_prefix + ".splicing_mutation.count_summary.pruned.txt", 
                              sconf.weights, args.effect_size_thres)

    utils.check_significance(args.output_prefix + ".splicing_mutation.count_summary.pruned.txt",
                             args.output_prefix + ".splicing_mutation.count_summary.BIC.txt",
                             sconf.weights, args.log_BF_thres, args.alpha0, args.beta0, args.alpha1, args.beta1)

    utils.add_annotation2(args.output_prefix + ".splicing_mutation.count_summary.BIC.txt",
                         args.output_prefix + ".splicing_mutation.count_summary.anno.txt",
                         sconf.sample_names,
                         args.output_prefix + ".splicing_mutation.link_info.txt")


    # permutation
    for i in range(args.permutation_num):

        print >> sys.stderr, "evaluating permutation " + str(i)
        utils.permute_mut_SJ_pairs(args.output_prefix + ".splicing_mutation.count_summary.txt",
                             args.output_prefix + ".splicing_mutation.count_summary.perm" + str(i) + ".txt")
   
        utils.convert_pruned_file(args.output_prefix + ".splicing_mutation.count_summary.perm" + str(i) + ".txt",
                                  args.output_prefix + ".splicing_mutation.count_summary.pruned.perm" + str(i) + ".txt",
                                  sconf.weights, args.effect_size_thres)

        utils.check_significance(args.output_prefix + ".splicing_mutation.count_summary.pruned.perm" + str(i) + ".txt",
                                 args.output_prefix + ".splicing_mutation.count_summary.BIC.perm" + str(i) + ".txt",
                                 sconf.weights, args.log_BF_thres, args.alpha0, args.beta0, args.alpha1, args.beta1)

        utils.add_annotation2(args.output_prefix + ".splicing_mutation.count_summary.BIC.perm" + str(i) + ".txt",
                              args.output_prefix + ".splicing_mutation.count_summary.anno.perm" + str(i) + ".txt",
                              sconf.sample_names,
                              args.output_prefix + ".splicing_mutation.link_info.txt")

    utils.merge_perm_result(args.output_prefix + ".splicing_mutation.count_summary.anno.perm", 
                            args.output_prefix + ".splicing_mutation.count_summary.anno.perm_all.txt",
                            args.permutation_num)

    utils.calculate_q_value(args.output_prefix + ".splicing_mutation.count_summary.anno.txt",
                            args.output_prefix + ".splicing_mutation.count_summary.anno.perm_all.txt",
                            args.output_prefix + ".savnet.result.txt",
                            args.permutation_num, args.sv)

    if args.debug == False:

        subprocess.call(["rm", "-rf", args.output_prefix + ".mut_merged.txt"])
        subprocess.call(["rm", "-rf", args.output_prefix + ".SJ_merged.txt"])
        subprocess.call(["rm", "-rf", args.output_prefix + ".SJ_merged.annot.txt"])
        subprocess.call(["rm", "-rf", args.output_prefix + ".SJ_merged.associate.txt"])
        subprocess.call(["rm", "-rf", args.output_prefix + ".IR_merged.txt"])
        subprocess.call(["rm", "-rf", args.output_prefix + ".IR_merged.associate.txt"])
        subprocess.call(["rm", "-rf", args.output_prefix + ".splicing_mutation.link_info.txt"])
        subprocess.call(["rm", "-rf", args.output_prefix + ".splicing.associate.txt"])
        subprocess.call(["rm", "-rf", args.output_prefix + ".splicing_mutation.mut_info.txt"])
        subprocess.call(["rm", "-rf", args.output_prefix + ".splicing_mutation.splicing_info.txt"])
        subprocess.call(["rm", "-rf", args.output_prefix + ".splicing_mutation.count_summary.txt"])
        subprocess.call(["rm", "-rf", args.output_prefix + ".splicing_mutation.count_summary.pruned.txt"])
        subprocess.call(["rm", "-rf", args.output_prefix + ".splicing_mutation.count_summary.BIC.txt"])
        subprocess.call(["rm", "-rf", args.output_prefix + ".splicing_mutation.count_summary.anno.txt"]) 

        for i in range(args.permutation_num):
            subprocess.call(["rm", "-rf", args.output_prefix + ".splicing_mutation.count_summary.perm" + str(i) + ".txt"])
            subprocess.call(["rm", "-rf", args.output_prefix + ".splicing_mutation.count_summary.pruned.perm" + str(i) + ".txt"])
            subprocess.call(["rm", "-rf", args.output_prefix + ".splicing_mutation.count_summary.BIC.perm" + str(i) + ".txt"])
            subprocess.call(["rm", "-rf", args.output_prefix + ".splicing_mutation.count_summary.anno.perm" + str(i) + ".txt"])

    """

