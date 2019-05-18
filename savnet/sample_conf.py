#! /usr/bin/env python

from __future__ import print_function

import sys, os

class Sample_conf(object):

    def __init__(self):

        self.sample_names = []
        self.mut_files = []
        self.sv_files = [] 
        self.SJ_files = [] 
        self.IR_files = [] 
        self.chimera_files = []
        self.weights = [] 

    def parse_file(self, file_path, sv_mode):

        header2ind = {}
        with open(file_path) as hin:

            # read header
            header = hin.readline().rstrip('\n').split('\t')
            for (i, cname) in enumerate(header):
    
                if sv_mode == False:
                    if cname not in ["Sample_Name", "Mutation_File", "SJ_File", "IR_File", "Weight"]:
                        print("header should any of Sample_Name, Mutation_File, SJ_File, IR_File or Weight", file = sys.stderr)
                        sys.exit(1)
                else:
                    if cname not in ["Sample_Name", "SV_File", "SJ_File", "IR_File", "Chimera_File", "Weight"]:
                        print("header should any of Sample_Name, SV_File, SJ_File, IR_File, Chimera_File or Weight", file = sys.stderr)
                        sys.exit(1)
                header2ind[cname] = i
        
            # check header
            if "Sample_Name" not in header2ind:
                print("Sample_Name column should be included in the sample config file", file = sys.stderr)
                sys.exit(1)

            if sv_mode == False and "Mutation_File" not in header2ind:
                print("Mutation_File column should be included in the sample config file", file = sys.stderr)
                sys.exit(1)

            if sv_mode == True and "SV_File" not in header2ind:
                print("SV_File column should be included in the sample config file", file = sys.stderr)
                sys.exit(1)

            if sv_mode == False:
                if "SJ_File" not in header2ind and "IR_File" not in header2ind:
                    print("At leaset one of SJ_File column or IR_File column should be included inthe smaple config file", file = sys.stderr)
                    sys.exit(1)
                if "SJ_File" not in header2ind and "IR_File" not in header2ind and "Chimera_File" not in header2ind:
                    print("At leaset one of SJ_File, IR_File or Chimera_File column should be included inthe smaple config file", file = sys.stderr)
                    sys.exit(1)

        
            for line in hin:
                F = line.rstrip('\n').split('\t')

                self.sample_names.append(F[header2ind["Sample_Name"]])

                if sv_mode == False:
                    if not os.path.exists(F[header2ind["Mutation_File"]]):
                        print(F[header2ind["Mutation_File"]] + " does not exist", file = sys.stderr)
                        sys.exit(1)
                    else:
                        self.mut_files.append(F[header2ind["Mutation_File"]])
                else:
                    if not os.path.exists(F[header2ind["SV_File"]]):
                        print(F[header2ind["SV_File"]] + " does not exist", file = sys.stderr)
                        sys.exit(1)
                    else:
                        self.sv_files.append(F[header2ind["SV_File"]])

                if "SJ_File" in header2ind:
                    if not os.path.exists(F[header2ind["SJ_File"]]):
                        print(F[header2ind["SJ_File"]] + " does not exist", file = sys.stderr)
                        sys.exit(1)
                    else:
                        self.SJ_files.append(F[header2ind["SJ_File"]])

                if "IR_File" in header2ind:
                    if not os.path.exists(F[header2ind["IR_File"]]):
                        print(F[header2ind["IR_File"]] + " does not exist", file = sys.stderr)
                        sys.exit(1)
                    else:
                        self.IR_files.append(F[header2ind["IR_File"]])

                if "Chimera_File" in header2ind:
                    if not os.path.exists(F[header2ind["Chimera_File"]]):
                        print(F[header2ind["Chimera_File"]] + " does not exist", file = sys.stderr)
                        sys.exit(1)
                    else:
                        self.chimera_files.append(F[header2ind["Chimera_File"]])

                if "Weight" in header2ind:
                    self.weights.append(float(F[header2ind["Weight"]]))
                    # self.weights.append(0.10)
                else:
                    self.weights.append(1.0)
