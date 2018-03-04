Quick Start
===========

Here, we use somatic mutations, splicing junctions and intron retentions collected from 26 lung cellline 

1. Change the directory to the workspace for this quick start

.. code-block:: bash

  % mkdir -p savnet_quick_start/resource
  % cd savnet_quick_start
  

2. Download the necessary resources for savnet. 

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

  # Script for creating input list file
  % wget https://storage.googleapis.com/friend1ws_package_data/savnet/make_savnet_input.py
  
  
4. Execute savnet

.. code-block:: bash

  % savnet resource/savnet_input.txt lung_cellline/lung_cellline resource/GRCh37.fa
  
  
5. Confirm the output file

.. code-block:: bash

  % cat lung_cellline/lung_cellline.savnet.result.txt
  
  
