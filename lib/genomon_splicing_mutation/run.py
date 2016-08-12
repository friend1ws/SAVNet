#! /usr/bin/env python

import sys, subprocess, os
import utils

def main(args):

    output_prefix_dir = os.path.dirname(args.output_prefix)
    if output_prefix_dir != "" and not os.path.exists(output_prefix_dir):
       os.makedirs(output_prefix_dir)

    utils.merge_mut(args.sample_list_file, args.output_prefix + ".mut_merged.txt")
    
    ##########
    # splicing_junction
    utils.merge_SJ2(args.sample_list_file, args.output_prefix + ".SJ_merged.txt", args.SJ_pooled_control_file, 2)

    subprocess.call(["junc_utils", "annotate", args.output_prefix + ".SJ_merged.txt", args.output_prefix + ".SJ_merged.annot.txt", args.resource_dir])

    subprocess.call(["junc_utils", "associate", args.output_prefix + ".mut_merged.txt", args.output_prefix + ".SJ_merged.annot.txt",
                     args.output_prefix, args.resource_dir, "--reference_genome", args.reference_genome, "-f", "anno"])

    # utils.add_gene_symbol(args.output_prefix + ".splicing_mutation.txt", args.output_prefix + ".splicing_mutation.proc.txt")
    ##########

    ########## 
    # intron_retention
    utils.merge_intron_retention(args.sample_list_file, args.output_prefix + ".IR_merged.txt", args.IR_pooled_control_file, 0.05, 3)

    subprocess.call(["genomon_intron_retention", "associate", args.output_prefix + ".IR_merged.txt", 
                     args.output_prefix + ".mut_merged.txt", args.output_prefix + ".IR_merged.associate.txt",
                     "--reference_genome", args.reference_genome, "--mutation", "anno" ])
    #########

    utils.merge_SJ_IR_files(args.output_prefix + ".splicing_mutation.txt", 
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
                              args.output_prefix + ".splicing_mutation.count_summary.pruned.txt", 3.0)

    utils.check_significance(args.output_prefix + ".splicing_mutation.count_summary.pruned.txt",
                             args.output_prefix + ".splicing_mutation.count_summary.BIC.txt")

    utils.summarize_result(args.output_prefix + ".splicing_mutation.count_summary.BIC.txt",
                           args.output_prefix + ".genomon_splicing_mutation.result.txt",
                           args.sample_list_file,
                           args.output_prefix + ".splicing_mutation.mut_info.txt",
                           args.output_prefix + ".splicing_mutation.splicing_info.txt")
   
    # permutation
    for i in range(10):

        print >> sys.stderr, "evaluating permutation " + str(i)
    
        utils.permute_mut_SJ_pairs(args.output_prefix + ".splicing_mutation.count_summary.txt",
                             args.output_prefix + ".splicing_mutation.count_summary.perm" + str(i) + ".txt")
    
        utils.convert_pruned_file(args.output_prefix + ".splicing_mutation.count_summary.perm" + str(i) + ".txt",
                                  args.output_prefix + ".splicing_mutation.count_summary.pruned.perm" + str(i) + ".txt", 3.0)
        
        utils.check_significance(args.output_prefix + ".splicing_mutation.count_summary.pruned.perm" + str(i) + ".txt",
                                 args.output_prefix + ".splicing_mutation.count_summary.BIC.perm" + str(i) + ".txt")

        
        utils.summarize_result(args.output_prefix + ".splicing_mutation.count_summary.BIC.perm" + str(i) + ".txt",
                               args.output_prefix + ".genomon_splicing_mutation.result.perm" + str(i) + ".txt",
                               args.sample_list_file,
                               args.output_prefix + ".splicing_mutation.mut_info.txt",
                               args.output_prefix + ".splicing_mutation.splicing_info.txt")

