#! /usr/bin/env python

import sys, glob, pysam

SJ_list_file = sys.argv[1]
control_file = sys.argv[2]

control_db = pysam.TabixFile(control_file)

file_list = [] 
with open(SJ_list_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        file_list.append(F[2])


junc2list = {}
for file in file_list:
    print >> sys.stderr, "proceccing: " + file
    with open(file, 'r') as hin:
        for line in hin:

            F = line.rstrip('\n').split('\t')
            if F[5] != "0": continue
            if int(F[6]) < 2: continue
            key = F[0] + '\t' + F[1] + '\t' + F[2]
            if key not in junc2list: junc2list[key] = 1


junc2count = {}
junc2control = {}
ind = 0
for file in file_list:
    print >> sys.stderr, "processing: " + file
    with open(file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            key = F[0] + '\t' + F[1] + '\t' + F[2]

            if F[1] == "10003574" and F[2] ==  "10009695":
                pass

            if key not in junc2list: continue

            if key in junc2control:
                if junc2control[key] == 1: continue
            else:

                ##########
                # remove control files
                tabixErrorFlag = 0
                try:
                    records = control_db.fetch(F[0], int(F[1]) - 5, int(F[1]) + 5)
                except Exception as inst:
                    print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
                    tabixErrorMsg = str(inst.args)
                    tabixErrorFlag = 1

                control_flag = 0;
                if tabixErrorFlag == 0:
                    for record_line in records:
                        record = record_line.split('\t')
                        if F[0] == record[0] and F[1] == record[1] and F[2] == record[2]:
                            control_flag = 1

                if control_flag == 1: 
                    junc2control[key] = 1
                    continue
                else: junc2control[key] = 0

            if key not in junc2count: junc2count[key] = ["0"] * len(file_list)
            junc2count[key][ind] = F[6]

    ind = ind + 1

for junc in sorted(junc2count):
    print junc + '\t' + ','.join(junc2count[junc])


