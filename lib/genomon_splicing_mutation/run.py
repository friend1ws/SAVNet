#! /usr/bin/env python

import subprocess, os
import utils

def main(args):

    output_prefix_dir = os.path.dirname(args.output_prefix)
    if output_prefix_dir != "" and not os.path.exists(output_prefix_dir):
       os.makedirs(output_prefix_dir)

    utils.merge_SJ(args.sample_list_file, args.output_prefix + ".SJ_merged.txt", args.pooled_control_file, 2)

    utils.merge_mut(args.sample_list_file, args.output_prefix + ".mut_merged.txt")

    subprocess.call(["junc_utils", "annotate", args.output_prefix + ".SJ_merged.txt", args.output_prefix + ".SJ_merged.annot.txt", args.resource_dir])

    subprocess.call(["junc_utils", "associate", args.output_prefix + ".mut_merged.txt", args.output_prefix + ".SJ_merged.annot.txt",
                     args.output_prefix, args.resource_dir, "--reference_genome", args.reference_genome, "-f", "anno"])

    utils.add_gene_symbol(args.output_prefix + ".splicing_mutation.txt", args.output_prefix + ".splicing_mutation.proc.txt")

    utils.organize_mut_SJ_count(args.output_prefix + ".splicing_mutation.proc.txt",
                                args.output_prefix + ".mut_merged.txt",
                                args.output_prefix + ".splicing_mutation.proc.count_summary.txt",
                                args.output_prefix + ".splicing_mutation.proc.mut_info.txt", 
                                args.output_prefix + ".splicing_mutation.proc.SJ_info.txt")

    utils.check_significance(args.output_prefix + ".splicing_mutation.proc.count_summary.txt",
                             args.output_prefix + ".splicing_mutation.proc.count_summary.BIC.txt")


    utils.summarize_result(args.output_prefix + ".splicing_mutation.proc.count_summary.BIC.txt",
                           args.output_prefix + ".genomon_splicing_mutation.result.txt",
                           args.sample_list_file,
                           args.output_prefix + ".splicing_mutation.proc.mut_info.txt",
                           args.output_prefix + ".splicing_mutation.proc.SJ_info.txt")


