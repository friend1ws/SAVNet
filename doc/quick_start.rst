Quick Start
===========

Here, we use somatic mutations, splicing junctions and intron retentions
collected from 26 lung cell-line whole genome and transcriptome sequencing data
by `Suzuki et al. (NAR, 2014) <https://doi.org/10.1093/nar/gku885>`_.


1. Change the directory to the workspace for this quick start.

.. code-block:: none

  % mkdir -p savnet_quick_start/resource
  % cd savnet_quick_start/resource


2. Download the necessary resources processed for SAVNet execution.

.. code-block:: none

  # Somatic mutation
  % wget https://zenodo.org/api/files/6dd7ee0c-4dce-45fe-9459-478fcd0e8102/vcf.tar.gz
  % tar -zxvf vcf.tar.gz

  # Splicing junction
  % wget https://zenodo.org/api/files/6dd7ee0c-4dce-45fe-9459-478fcd0e8102/junction.tar.gz
  % tar -zxvf junction.tar.gz

  # Intron retention
  % wget https://zenodo.org/api/files/6dd7ee0c-4dce-45fe-9459-478fcd0e8102/intron_retention.tar.gz
  % tar -zxvf intron_retention.tar.gz

  # Quality check
  % wget https://zenodo.org/api/files/6dd7ee0c-4dce-45fe-9459-478fcd0e8102/qc.tar.gz
  % tar -zxvf qc.tar.gz


When using Annovar input file format file (optional):

.. code-block:: none

  # Somatic mutation in Annovar input file format
  % wget https://zenodo.org/api/files/6dd7ee0c-4dce-45fe-9459-478fcd0e8102/annovar.tar.gz
  % tar -zxvf annovar.tar.gz

  # Reference genome
  % wget ftp://ftp.ncbi.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh37/special_requests/GRCh37-lite.fa.gz


  
  
3. Download the script to create input list file, and execute it.

.. code-block:: none

  # Download the script for creating input list file
  % wget https://zenodo.org/api/files/6dd7ee0c-4dce-45fe-9459-478fcd0e8102/make_savnet_input.py

  # Generate the input list file by running the script
  % python make_savnet_input.py --sample_list_file savnet.input.txt --mut_dir vcf --sj_dir junction --ir_dir intron_retention --qc_dir qc

  # Confirm the input list file
  % cat savnet.input.txt

  #
  # Following keys will be shown
  #

  Sample_Name     Weight  Mutation_File   SJ_File IR_File
  A427    5.2935  /home/friend1ws/savnet_quick_start/resource/vcf/A427.vcf  /home/friend1ws/savnet_quick_start/resource/junction/A427.SJ.out.tab    /home/friend1ws/savnet_quick_start/resource/intron_retention/A427.intron_retention.txt
  A549    2.5843  /home/friend1ws/savnet_quick_start/resource/vcf/A549.vcf  /home/friend1ws/savnet_quick_start/resource/junction/A549.SJ.out.tab    /home/friend1ws/savnet_quick_start/resource/intron_retention/A549.intron_retention.txt
  ABC-1   4.6718  /home/friend1ws/savnet_quick_start/resource/vcf/ABC-1.vcf /home/friend1ws/savnet_quick_start/resource/junction/ABC-1.SJ.out.tab   /home/friend1ws/savnet_quick_start/resource/intron_retention/ABC-1.intron_retention.txt


When using Annovar input file format file (optional):

.. code-block:: none

  % python make_savnet_input.py --sample_list_file savnet.input.anno.txt --mut_dir annovar --sj_dir junction --ir_dir intron_retention --qc_dir qc

  #
  # Following keys will be shown
  #

  Sample_Name     Weight  Mutation_File   SJ_File IR_File
  A427    5.2935  /home/friend1ws/savnet_quick_start/resource/annovar/A427.avinput  /home/friend1ws/savnet_quick_start/resource/junction/A427.SJ.out.tab    /home/friend1ws/savnet_quick_start/resource/intron_retention/A427.intron_retention.txt
  A549    2.5843  /home/friend1ws/savnet_quick_start/resource/annovar/A549.avinput  /home/friend1ws/savnet_quick_start/resource/junction/A549.SJ.out.tab    /home/friend1ws/savnet_quick_start/resource/intron_retention/A549.intron_retention.txt
  ABC-1   4.6718  /home/friend1ws/savnet_quick_start/resource/annovar/ABC-1.avinput /home/friend1ws/savnet_quick_start/resource/junction/ABC-1.SJ.out.tab   /home/friend1ws/savnet_quick_start/resource/intron_retention/ABC-1.intron_retention.txt



4. Execute SAVNet

.. code-block:: none

  % cd ../
  % savnet resource/savnet.input.txt lung_cellline/lung_cellline

It will take 10 to 20 minutes for completing the calculation.

When using Annovar input file format file (optional):

.. code-block:: none

  % cd ../
  % savnet resource/savnet.input.anno.txt lung_cellline/lung_cellline --reference resource/GRCh37.fa


5. Confirm the output file.

.. code-block:: none

  % cat lung_cellline/lung_cellline.savnet.result.txt

  #
  # Following keys will be shown
  #

  Gene_Symbol     Sample_Name     Mutation_Key    Motif_Pos       Mutation_Type   Is_Canonical    Splicing_Key    Splicing_Class  Is_Inframe      Supporting_Read_Num     Score   Q_Value
  ABCC9   RERF-LC-Ad1     12,21981996,T,A 12:21981994-21982000,-  Acceptor disruption     Canonical       12:21981983-21991011    Alternative 3'SS        In-frame        34      100.2136        0.0247
  ABCD4   RERF-LC-Ad1     14,74754513,C,G 14:74754507-74754515,-  Donor disruption        Non-canonical   14:74753520-74754520    Alternative 5'SS        ---     102     300.2769        0.02
  ABCD4   RERF-LC-Ad1     14,74754513,C,G 14:74754507-74754515,-  Donor disruption        Non-canonical   14:74753520-74754909    Exon skipping   ---     3       300.2769        0.02
  ABLIM3  H1648   5,148630908,T,A 5:148630904-148630910,+ Acceptor creation       Canonical       5:148630068-148630908   Alternative 3'SS        In-frame        5       11.323  0.051
