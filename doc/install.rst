Install
##########

To use savnet, ``htslib`` and ``bedtools`` need to be installed and added to the PATH.

.. code-block:: bash

  % wget https://github.com/samtools/htslib/releases/download/1.7/htslib-1.7.tar.bz2
  % tar jxvf htslib-1.7.tar.bz2
  % cd htslib-1.7 && make && export PATH=$PATH:$PWD && cd ..

  % wget https://github.com/arq5x/bedtools2/releases/download/v2.27.0/bedtools-2.27.0.tar.gz
  % tar zxvf bedtools-2.27.0.tar.gz
  % cd bedtools2 && make && export PATH=$PATH:$PWD/bin && cd ..


You may want to check whether they are installed correctly
through typing following commands:

.. code-block:: bash

  % tabix

  # Following messages will appear
  Version: 1.7
  Usage:   tabix [OPTIONS] [FILE] [REGION [...]]

  Indexing Options:

  % bedtools

  # Following messages will appear
  bedtools is a powerful toolset for genome arithmetic.

  Version:   v2.27.0
  About:     developed in the quinlanlab.org and by many contributors worldwide.


If help messages appear, then installation seems to be successful.

Next, install the savnet package using ``pip``

.. code-block:: bash

  % pip install savnet

Then, other dependent python packages
(``pysam``, ``annot_utils``, ``junc_utils``, ``intron_retention_utils``, ``chimera_utils``)
will also be installed.
If you are using shared computational environment where you do not have root account,
you may want to add the option ``--user``.

To confirm successful installation, please type:

.. code-block:: bash

  % savnet --version

  # Following messages will appear
  savnet 0.3.0


Alternatively, if docker is installed in your environment,
you can directly use the savnet through docker.

.. code-block:: bash

  % docker run friend1ws/savnet:0.3.0 --version
