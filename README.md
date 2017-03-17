# SAVNet

## Introduction

 software for extracting splicing associated variants from somatic mutation, splicing junction and intron retention data.

## Dependency

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


