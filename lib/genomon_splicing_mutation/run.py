#! /usr/bin/env python

import utils

def main(args):

    utils.merge_SJ(args.sample_list_file, args.output_prefix + ".SJ_merged.txt", args.pooled_control_file, 2)

 
