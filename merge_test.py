#! /usr/bin/env python

import sys
import junc_utils.utils

SJ_list_file = sys.argv[1]
output_file = sys.argv[2]
control_file = sys.argv[3]

file_list = []
with open(SJ_list_file, 'r') as hin:
    for line in hin:
        in_file = line.rstrip('\n')
        print >> sys.stderr, "processing: " + in_file
        junc_utils.utils.proc_star_junction(in_file, output_file + ".tmp", control_file,
                                            2, 10, True, False)


