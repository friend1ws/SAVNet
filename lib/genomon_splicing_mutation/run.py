#! /usr/bin/env python

import os
import utils

def main(args):

    output_prefix_dir = os.path.dirname(args.output_prefix)
    if output_prefix_dir != "" and not os.path.exists(output_prefix_dir):
       os.makedirs(output_prefix_dir)

    # utils.merge_SJ(args.sample_list_file, args.output_prefix + ".SJ_merged.txt", args.pooled_control_file, 2)

    utils.merge_mut(args.sample_list_file, args.output_prefix + ".mut_merged.txt")

