Quick Start
===========

Here, we use somatic mutations, splicing junctions and intron retentions
collected from 26 lung cell-line whole genome and transcriptome sequencing data 
by `Suzuki et al. (NAR, 2014) <https://doi.org/10.1093/nar/gku885>`_.

1. Change the directory to the workspace for this quick start

.. code-block:: bash

  % mkdir -p savnet_quick_start/resource
  % cd savnet_quick_start


2. Download the necessary resources processed for SAVNet execution.

.. code-block:: bash

  # Somatic mutation
  % wget https://storage.googleapis.com/friend1ws_package_data/savnet/mutation.tar.gz

  # Splicing junction
  % wget https://storage.googleapis.com/friend1ws_package_data/savnet/junction.tar.gz

  # Intron retention
  % wget https://storage.googleapis.com/friend1ws_package_data/savnet/intron_retention.tar.gz

  # Reference genome
  % wget https://storage.googleapis.com/friend1ws_package_data/common/GRCh37.fa
  
  
3. Download the script to create input list file, and execute it.

.. code-block:: bash

  # Download the script for creating input list file
  % wget https://storage.googleapis.com/friend1ws_package_data/savnet/make_savnet_input.py
 
  # Generate the input list file by running the script
  % python make_savnet_input.py --sample_list_file savnet.input.txt --mut_dir mutation --sj_dir junction --ir_dir intron_retention --qc_dir qc

  # Confirm the input list file
  % cat savnet.input.txt

  # Following keys will be shown
  Sample_Name     Weight  Mutation_File   SJ_File IR_File
  A427    5.2935  /home/friend1ws/savnet_quick_start/resource/mutation/A427.mutation.txt  /home/friend1ws/savnet_quick_start/resource/junction/A427.SJ.out.tab    /home/friend1ws/savnet_quick_start/resource/intron_retention/A427.genomonIR.result.txt
  A549    2.5843  /home/friend1ws/savnet_quick_start/resource/mutation/A549.mutation.txt  /home/friend1ws/savnet_quick_start/resource/junction/A549.SJ.out.tab    /home/friend1ws/savnet_quick_start/resource/intron_retention/A549.genomonIR.result.txt
  ABC-1   4.6718  /home/friend1ws/savnet_quick_start/resource/mutation/ABC-1.mutation.txt /home/friend1ws/savnet_quick_start/resource/junction/ABC-1.SJ.out.tab   /home/friend1ws/savnet_quick_start/resource/intron_retention/ABC-1.genomonIR.result.txt


4. Execute Savnet

.. code-block:: bash

  % cd ../

  % savnet resource/savnet_input.txt lung_cellline/lung_cellline resource/GRCh37.fa



5. Confirm the output file

.. code-block:: bash

  % cat lung_cellline/lung_cellline.savnet.result.txt
