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
    sconf.parse_file(args.sample_list_file, args.sv)

    if args.sv == False:

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
        subprocess.call(["intron_retention_utils", "associate", args.output_prefix + ".IR_merged.txt", 
                         args.output_prefix + ".mut_merged.txt", args.output_prefix + ".IR_merged.associate.txt",
                         "--reference_genome", args.reference_genome, "--mutation", "anno"])
        #########
        utils.merge_SJ_IR_files(args.output_prefix + ".SJ_merged.associate.txt", 
                                args.output_prefix + ".IR_merged.associate.txt",
                                args.output_prefix + ".splicing.associate.txt")
        
        utils.organize_mut_splicing_count(args.output_prefix + ".splicing.associate.txt",
                                    args.output_prefix + ".mut_merged.txt",
                                    args.output_prefix + ".splicing_mutation.count_summary.txt",
                                    args.output_prefix + ".splicing_mutation.mut_info.txt", 
                                    args.output_prefix + ".splicing_mutation.splicing_info.txt")

    else:

        utils.merge_sv(sconf.sv_files, args.output_prefix + ".sv_merged.txt")

        ##########
        # splicing_junction
        """
        utils.merge_SJ2(sconf.SJ_files, args.output_prefix + ".SJ_merged.txt", args.SJ_pooled_control_file, args.SJ_num_thres)

        subprocess.call(["junc_utils", "annotate", args.output_prefix + ".SJ_merged.txt", args.output_prefix + ".SJ_merged.annot.txt", args.resource_dir])
        """
        subprocess.call(["junc_utils", "associate", args.output_prefix + ".SJ_merged.annot.txt", args.output_prefix + ".sv_merged.txt",
                         args.output_prefix + ".SJ_merged.associate.txt", args.resource_dir, "--sv"])
        ##########
        # intron_retention
        utils.merge_intron_retention(sconf.IR_files, args.output_prefix + ".IR_merged.txt",
                                     args.IR_pooled_control_file, args.IR_ratio_thres, args.IR_num_thres)

        subprocess.call(["intron_retention_utils", "associate", args.output_prefix + ".IR_merged.txt",
                         args.output_prefix + ".sv_merged.txt", args.output_prefix + ".IR_merged.associate.txt",
                         "--sv"])
        ##########
        # chimera
        utils.merge_chimera(sconf.chimera_files, args.output_prefix + ".chimera_merged.txt", 
                            args.chimera_pooled_control_file, args.chimera_num_thres, args.chimera_overhang_thres)
         
        subprocess.call(["chimera_utils", "associate", "--is_grc", args.output_prefix + ".chimera_merged.txt",
                         args.output_prefix + ".sv_merged.txt", args.output_prefix + ".chimera_merged.associate.txt"])
        ##########

        utils.merge_SJ_IR_chimera_files_sv(args.output_prefix + ".SJ_merged.associate.txt",
                                           args.output_prefix + ".IR_merged.associate.txt",
                                           args.output_prefix + ".chimera_merged.associate.txt",
                                           args.output_prefix + ".splicing.associate.txt")

        utils.organize_mut_splicing_count(args.output_prefix + ".splicing.associate.txt",
                                    args.output_prefix + ".sv_merged.txt",
                                    args.output_prefix + ".splicing_mutation.count_summary.txt",
                                    args.output_prefix + ".splicing_mutation.mut_info.txt",
                                    args.output_prefix + ".splicing_mutation.splicing_info.txt", True)


    # true combination    
    print >> sys.stderr, "evaluating true combinations"

    utils.convert_pruned_file(args.output_prefix + ".splicing_mutation.count_summary.txt",
                              args.output_prefix + ".splicing_mutation.count_summary.pruned.txt", 
                              sconf.weights, args.effect_size_thres)

    utils.check_significance(args.output_prefix + ".splicing_mutation.count_summary.pruned.txt",
                             args.output_prefix + ".splicing_mutation.count_summary.BIC.txt",
                             sconf.weights, args.log_BF_thres, args.alpha0, args.beta0, args.alpha1, args.beta1)

    utils.add_annotation(args.output_prefix + ".splicing_mutation.count_summary.BIC.txt",
                         args.output_prefix + ".splicing_mutation.count_summary.anno.txt",
                         sconf.sample_names,
                         args.output_prefix + ".splicing_mutation.mut_info.txt",
                         args.output_prefix + ".splicing_mutation.splicing_info.txt", args.sv)
   
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

        
        utils.add_annotation(args.output_prefix + ".splicing_mutation.count_summary.BIC.perm" + str(i) + ".txt",
                             args.output_prefix + ".splicing_mutation.count_summary.anno.perm" + str(i) + ".txt",
                             sconf.sample_names,
                             args.output_prefix + ".splicing_mutation.mut_info.txt",
                             args.output_prefix + ".splicing_mutation.splicing_info.txt", args.sv)

    utils.calculate_q_value(args.output_prefix + ".splicing_mutation.count_summary.anno.txt",
                            args.output_prefix + ".splicing_mutation.count_summary.anno.perm",
                            args.output_prefix + ".genomon_splicing_mutation.result.txt",
                            args.permutation_num, args.sv)

    if args.debug == False:
        """
        subprocess.call(["rm", "-rf", args.output_prefix + ".splicing.associate.txt"])
        subprocess.call(["rm", "-rf", args.output_prefix + ".splicing_mutation.mut_info.txt"])
        subprocess.call(["rm", "-rf", args.output_prefix + ".splicing_mutation.splicing_info.txt"])
        subprocess.call(["rm", "-rf", args.output_prefix + ".splicing_mutation.count_summary.txt"])
        subprocess.call(["rm", "-rf", args.output_prefix + ".splicing_mutation.count_summary.pruned.txt"])
        subprocess.call(["rm", "-rf", args.output_prefix + ".splicing_mutation.count_summary.BIC.txt"])
        subprocess.call(["rm", "-rf", args.output_prefix + ".splicing_mutation.count_summary.anno.txt"]) 
        """
        for i in range(args.permutation_num):
            subprocess.call(["rm", "-rf", args.output_prefix + ".splicing_mutation.count_summary.perm" + str(i) + ".txt"])
            subprocess.call(["rm", "-rf", args.output_prefix + ".splicing_mutation.count_summary.pruned.perm" + str(i) + ".txt"])
            subprocess.call(["rm", "-rf", args.output_prefix + ".splicing_mutation.count_summary.BIC.perm" + str(i) + ".txt"])
            subprocess.call(["rm", "-rf", args.output_prefix + ".splicing_mutation.count_summary.anno.perm" + str(i) + ".txt"])


