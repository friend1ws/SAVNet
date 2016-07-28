#! /usr/bin/env python

import sys, re, subprocess

input_file = sys.argv[1]
output_file = sys.argv[2]

header2ind = {}
header = ""
hout0 = open(output_file + ".tmp0", 'w')
hout1 = open(output_file + ".tmp1", 'w')
with open(input_file, 'r') as hin:

    header = hin.readline().rstrip('\n').split('\t')
    for i in range(len(header)):
        header2ind[header[i]] = i

    print >> hout0, "Gene_Symbol" + '\t' + '\t'.join(header)

    for line in hin:
        F = line.rstrip('\n').split('\t')
        genes = F[header2ind["Gene_1"]].split(';') + F[header2ind["Gene_2"]].split(';')
        genes = list(set(genes))

        if "---" in genes: genes.remove("---")
        if len(genes) > 0: 
            genes_nm = filter(lambda x: x.find("(NM_") > 0, genes)
            if len(genes_nm) > 0: genes = genes_nm

        gene = genes[0]
        gene = re.sub(r"\(N[MR]_\d+\)", "", gene)

        print >> hout1, gene + '\t' + '\t'.join(F)

hout0.close()
hout1.close()

hout2 = open(output_file + ".tmp2", 'w')
subprocess.call(["sort", "-k1", output_file + ".tmp1"], stdout = hout2)
hout2.close()

hout3 = open(output_file, 'w')
subprocess.call(["cat", output_file + ".tmp0", output_file + ".tmp2"], stdout = hout3)
hout3.close()

subprocess.call(["rm", "-rf", output_file + ".tmp0"])
subprocess.call(["rm", "-rf", output_file + ".tmp1"])
subprocess.call(["rm", "-rf", output_file + ".tmp2"])

