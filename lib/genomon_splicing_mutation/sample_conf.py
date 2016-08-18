#! /usr/bin/env python

import sys, os

class Sample_conf(object):

    def __init__(self):

        self.sample_names = []
        self.mut_files = [] 
        self.SJ_files = [] 
        self.IR_files = [] 
        self.weights = [] 

    def parse_file(self, file_path):

        header2ind = {}
        with open(file_path) as hin:

            # read header
            header = hin.readline().rstrip('\n').split('\t')
            for (i, cname) in enumerate(header):
                if cname not in ["Sample_Name", "Mutation_File", "SJ_File", "IR_File", "Weight"]:
                    print >> sys.stderr, "header should any of Sample_Name, Mutation_File, SJ_File, IR_File or Weight"
                    sys.exit(1)
                header2ind[cname] = i
        
            # check header
            if "Sample_Name" not in header2ind:
                print >> sys.stderr, "Sample_Name column should be included in the sample config file"
                sys.exit(1)

            if "Mutation_File" not in header2ind:
                print >> sys.stderr, "Mutation_File column should be included in the sample config file"
                sys.exit(1)

            if "SJ_File" not in header2ind and "IR_File" not in header2ind:
                print >> sys.stderr, "At leaset one of SJ_File column or IR_File column should be included inthe smaple config file"
                sys.exit(1)

        
            for line in hin:
                F = line.rstrip('\n').split('\t')

                self.sample_names.append(F[header2ind["Sample_Name"]])

                if not os.path.exists(F[header2ind["Mutation_File"]]):
                    print >> sys.stderr, F[header2ind["Mutation_File"]] + " does not exist"
                    sys.exit(1)
                else:
                    self.mut_files.append(F[header2ind["Mutation_File"]])

                if "SJ_File" in header2ind:
                    if not os.path.exists(F[header2ind["SJ_File"]]):
                        print >> sys.stderr, F[header2ind["SJ_File"]] + " does not exist" 
                        sys.exit(1)
                    else:
                        self.SJ_files.append(F[header2ind["SJ_File"]])

                if "IR_File" in header2ind:
                    if not os.path.exists(F[header2ind["IR_File"]]):
                        print >> sys.stderr, F[header2ind["IR_File"]] + " does not exist"
                        sys.exit(1)
                    else:
                        self.IR_files.append(F[header2ind["IR_File"]])

                if "Weight" in header2ind:
                    self.weights.append(float(F[header2ind["Weight"]]))
                else:
                    self.weights.append(1.0)
