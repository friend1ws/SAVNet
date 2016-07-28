#! /bin/bash

<<_COM_

echo "python generate_omega_mut_list.py /home/omega3/omega_project/genomon2_2_0_alpha/DLBC/mutation /home/omega3/omega_rna/star/DLBC/ 01 > DLBC_mut_SJ_list.txt"
python generate_omega_mut_list.py /home/omega3/omega_project/genomon2_2_0_alpha/DLBC/mutation /home/omega3/omega_rna/star/DLBC/ 01 > DLBC_mut_SJ_list.txt

echo "python merge_SJ.py DLBC_mut_SJ_list.txt /home/yshira/project/inframe_junc/output/control.bed.gz > DLBC_SJ_merged.txt "
python merge_SJ.py DLBC_mut_SJ_list.txt /home/yshira/project/inframe_junc/output/control.bed.gz > DLBC_SJ_merged.txt 


echo "python merge_mut.py DLBC_mut_list.txt DLBC_mut2sample.txt"
python merge_mut.py DLBC_mut_list.txt DLBC_mut2sample.txt

echo "junc_utils annotate DLBC_SJ_merged.txt DLBC_SJ_merged.annot.txt /home/yshira/mysoftware/junc_utils/resource/"
junc_utils annotate DLBC_SJ_merged.txt DLBC_SJ_merged.annot.txt /home/yshira/mysoftware/junc_utils/resource/

echo "junc_utils associate DLBC_mut2sample.txt DLBC_SJ_merged.annot.txt DLBC /home/yshira/mysoftware/junc_utils/resource --reference_genome /home/w3varann/database/GRCh37/GRCh37.fa -f anno"
junc_utils associate DLBC_mut2sample.txt DLBC_SJ_merged.annot.txt DLBC /home/yshira/mysoftware/junc_utils/resource --reference_genome /home/w3varann/database/GRCh37/GRCh37.fa -f anno

_COM_

echo "python proc_splicing_mutation.py DLBC.splicing_mutation.txt DLBC.splicing_mutation.proc.txt"
python proc_splicing_mutation.py DLBC.splicing_mutation.txt DLBC.splicing_mutation.proc.txt


echo "python summarize_mut_SJ.py DLBC.splicing_mutation.proc.txt DLBC_mut2sample.txt DLBC.splicing_mutation.rawdata.txt DLBC.splicing_mutation.mutation.txt DLBC.splicing_mutation.splicing.txt"
python summarize_mut_SJ.py DLBC.splicing_mutation.proc.txt DLBC_mut2sample.txt DLBC.splicing_mutation.rawdata.txt DLBC.splicing_mutation.mutation.txt DLBC.splicing_mutation.splicing.txt

echo "python calc_BIC.py DLBC.splicing_mutation.rawdata.txt > DLBC.splicing_mutation.bic.txt"
python calc_BIC.py DLBC.splicing_mutation.rawdata.txt > DLBC.splicing_mutation.bic.txt

python organize_result.py DLBC.splicing_mutation.bic.txt DLBC_mut_SJ_list.txt DLBC.splicing_mutation.mutation.txt DLBC.splicing_mutation.splicing.txt > DLBC.splicing_mutation.final_result.txt



