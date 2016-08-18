#! /usr/bin/env python

import sys, subprocess, os
import utils, sample_conf

def main(args):

    output_prefix_dir = os.path.dirname(args.output_prefix)
    if output_prefix_dir != "" and not os.path.exists(output_prefix_dir):
       os.makedirs(output_prefix_dir)
    
    ##########
    # read sample conf
    sconf = sample_conf.Sample_conf()
    sconf.parse_file(args.sample_list_file)

    """
    utils.merge_mut(sconf.mut_files, args.output_prefix + ".mut_merged.txt")
    ##########
    # splicing_junction
    utils.merge_SJ2(sconf.SJ_files, args.output_prefix + ".SJ_merged.txt", args.SJ_pooled_control_file, args.SJ_num_thres)

    subprocess.call(["junc_utils", "annotate", args.output_prefix + ".SJ_merged.txt", args.output_prefix + ".SJ_merged.annot.txt", args.resource_dir])

    subprocess.call(["junc_utils", "associate", args.output_prefix + ".SJ_merged.annot.txt", args.output_prefix + ".mut_merged.txt", 
                     args.output_prefix + ".SJ_merged.associate.txt", args.resource_dir, "--reference_genome", args.reference_genome,
                     "--mutation_format", "anno"])

    ##########
    # intron_retention
    utils.merge_intron_retention(sconf.IR_files, args.output_prefix + ".IR_merged.txt", 
                                 args.IR_pooled_control_file, args.IR_ratio_thres, args.IR_num_thres)

    subprocess.call(["genomon_intron_retention", "associate", args.output_prefix + ".IR_merged.txt", 
                     args.output_prefix + ".mut_merged.txt", args.output_prefix + ".IR_merged.associate.txt",
                     "--reference_genome", args.reference_genome, "--mutation", "anno" ])
    #########
    """

    utils.merge_SJ_IR_files(args.output_prefix + ".SJ_merged.associate.txt", 
                            args.output_prefix + ".IR_merged.associate.txt",
                            args.output_prefix + ".splicing.associate.txt")

    utils.organize_mut_splicing_count(args.output_prefix + ".splicing.associate.txt",
                                args.output_prefix + ".mut_merged.txt",
                                args.output_prefix + ".splicing_mutation.count_summary.txt",
                                args.output_prefix + ".splicing_mutation.mut_info.txt", 
                                args.output_prefix + ".splicing_mutation.splicing_info.txt")
 
    # true combination    
    print >> sys.stderr, "evaluating true combinations"

    utils.convert_pruned_file(args.output_prefix + ".splicing_mutation.count_summary.txt",
                              args.output_prefix + ".splicing_mutation.count_summary.pruned.txt", 
                              sconf.weights, args.effect_size_thres)

    utils.check_significance(args.output_prefix + ".splicing_mutation.count_summary.pruned.txt",
                             args.output_prefix + ".splicing_mutation.count_summary.BIC.txt",
                             sconf.weights, args.log_BF_thres, args.alpha0, args.beta0, args.alpha1, args.beta1)

    utils.summarize_result(args.output_prefix + ".splicing_mutation.count_summary.BIC.txt",
                           args.output_prefix + ".genomon_splicing_mutation.result.txt",
                           sconf.sample_names,
                           args.output_prefix + ".splicing_mutation.mut_info.txt",
                           args.output_prefix + ".splicing_mutation.splicing_info.txt")
   
    # permutation
    for i in range(10):

        print >> sys.stderr, "evaluating permutation " + str(i)
    
        utils.permute_mut_SJ_pairs(args.output_prefix + ".splicing_mutation.count_summary.txt",
                             args.output_prefix + ".splicing_mutation.count_summary.perm" + str(i) + ".txt")
    
        utils.convert_pruned_file(args.output_prefix + ".splicing_mutation.count_summary.perm" + str(i) + ".txt",
                                  args.output_prefix + ".splicing_mutation.count_summary.pruned.perm" + str(i) + ".txt",
                                  sconf.weights, args.effect_size_thres)
        
        utils.check_significance(args.output_prefix + ".splicing_mutation.count_summary.pruned.perm" + str(i) + ".txt",
                                 args.output_prefix + ".splicing_mutation.count_summary.BIC.perm" + str(i) + ".txt",
                                 sconf.weights, args.log_BF_thres, args.alpha0, args.beta0, args.alpha1, args.beta1)

        
        utils.summarize_result(args.output_prefix + ".splicing_mutation.count_summary.BIC.perm" + str(i) + ".txt",
                               args.output_prefix + ".genomon_splicing_mutation.result.perm" + str(i) + ".txt",
                               sconf.sample_names,
                               args.output_prefix + ".splicing_mutation.mut_info.txt",
                               args.output_prefix + ".splicing_mutation.splicing_info.txt")

