# SAVNet

## Introduction

 software for extracting splicing associated variants (SAVs) from somatic mutation, splicing junction and intron retention data. 
 This software has been used for large-scale exome-transcriptome sequence analysis (see our [preprint](https://www.biorxiv.org/content/early/2017/09/28/162560)).

## Dependency

### Python

Python (>= 2.7), `pysam (>= 0.8.1)`, [`junc_utils`](https://github.com/friend1ws/junc_utils), [`intron retention_utils`](https://github.com/friend1ws/intron_retention_utils) packages.

### Software

[hstlib](http://www.htslib.org)

## Install 
```
git clone  https://github.com/friend1ws/SAVNet.git
cd SAVNet
python setup.py build install
```


## Input files

### How to prepare sample config file

The sample config file should be tab-delimited and start with headers.
Currently, there are 5 columns, **Sample_Name**, **Mutation_File**, **SJ_File**, **IR_File** and **Weight**.
Splicing junctions in **SJ_File** furnish evidences for *Exon-skipping*, *Alternative 5'SS (splice site)* and *Alternative 3'SS*,
whereas intron retentions in **IR_File** constitute evidences for *Intron retention*
**Sample_Name** and **Mutation_File** is required column.
Either **SJ_File** or **IR_File** need to be specified
(When, e.g., **IR_File** is not specified, SAVNet does not consider mutations associated with *intron retention*).
**Weight** is optional.

#### Sample_Name

The field for sample names. This is used for the output file to show which samples have the identified SAVs.

#### Mutation_File

The pathes to mutation calling data for each sample (on hg19, hg38 or mm10 reference genome). 
Currently, SAVNet accepts only [Annovar input file](http://annovar.openbioinformatics.org/en/latest/user-guide/input/) format
(We will soon update the software so that it can accept the vcf format).
Only the first five columns (Chr, Start, End, Ref, Alt) are used, and others are ignored.
When the chromosome names are not chr-prefix (when using Genome Reference Consortium (GRC) genomes such as GRCh37), 
add `--grc` option in execution of `savnet` command below.
Lines starting with # are skipped (as comments).


#### SJ_File

The pathes to the splicing junction files (SJ.out.tab) generated by the [`STAR`](https://github.com/alexdobin/STAR) alignment.
The splicing junctions appearing in SJ.out tab files greatly changes by STAR parameters  such as `--outSJfilter`.
We currently recommend to use the following parameters (together with other settings such as sorting and thread numbers):
```
-outSJfilterCountUniqueMin 1 1 1 1 --outSJfilterCountTotalMin 1 1 1 1 \
--outSJfilterOverhangMin 12 12 12 12 --outSJfilterDistToOtherSJmin 0 0 0 0 --alignSJstitchMismatchNmax -1 -1 -1 -1 
```

#### IR_File

The pathes to the intron retention files generated by [`intron_retention_utils`](https://github.com/friend1ws/intron_retention_utils)
`simple_count` command. We currenty recommed to use the fowwloing parameters:
```
--intron_retention_check_size 10 --mapping_qual_thres 20
```

#### Weight

This is used for negate the diversity of the numbers of total RNA-seq reads among samples.
Although this field is optional, we think you can obtain reasonable results without setting this.
However, we currently recommend to set this as the number of uniquely aligned read pairs (derived from the STAR Log.final.out files) divided by 10^7.


### How to prepare pooled control files

When either one or both of pooled control files for splicing junctions and intron retentions are specified, 
SAVNet removes the splicing variations registered in these files before associating mutations and splicings.
This will greatly help improving the accuracy of SAV calls, 
focusing on splicing variations that cannot be observed in normal control samples, as well as reducing the computational cost.
We recommend to specify these parameters using as many samples as possible (hopefully at least >= 10 control samples).

#### Splicing junction control files

The pooled control file for splicing junction can be generated by [`junc_utils`](https://github.com/friend1ws/junc_utils) 
`merge_control` command. We currently recommend to use the following parameters:
```
junc_utils merge_control --read_num_thres 2 --keep_annotated --sample_num_thres 1 ${input_list} ${output_file}
```
The value of `--sample_num_thres` can be tuned for large number of control samples.


#### Intron retention control files


## Commands

```
 savnet [-h] [--version] [--grc] [--genome_id {hg19,hg38,mm10}] [--sv]
              [--branchpoint] [--donor_size donor_size]
              [--acceptor_size acceptor_size]
              [--SJ_pooled_control_file SJ_POOLED_CONTROL_FILE]
              [--IR_pooled_control_file IR_POOLED_CONTROL_FILE]
              [--chimera_pooled_control_file CHIMERA_POOLED_CONTROL_FILE]
              [--SJ_num_thres SJ_NUM_THRES] [--keep_annotated]
              [--IR_num_thres IR_NUM_THRES] [--IR_ratio_thres IR_RATIO_THRES]
              [--chimera_num_thres CHIMERA_NUM_THRES]
              [--chimera_overhang_thres CHIMERA_OVERHANG_THRES]
              [--permutation_num PERMUTATION_NUM] [--alpha0 ALPHA0]
              [--beta0 BETA0] [--alpha1 ALPHA1] [--beta1 BETA1]
              [--log_BF_thres LOG_BF_THRES]
              [--effect_size_thres EFFECT_SIZE_THRES] [--debug]
              sample_list.txt output_prefix reference.fa
```


## About result
The following columns are added to the input files:

* **Gene_Symbol**: Gene symbol affected by the SAV
* **Sample_Name**: Samples having the SAV
* **Mutation_Key**: {chromosome},{position},{reference allele},{alternative allele}
* **Motif_Pos**: The position of disrupted or newly created splicing motif ({chromosome}:{start}-{end})
* **Mutation_Type**: splicing donor dirsuption, splicing acceptor disruption, splicing donor creation or splicing acceptor creation
