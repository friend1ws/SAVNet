FROM ubuntu:16.04
MAINTAINER Yuichi Shiraishi <friend1ws@gmail.com> 

RUN apt-get update && apt-get install -y \
    git \
    wget \
    bzip2 \
    make \
    gcc \
    zlib1g-dev \
    python \
    python-pip \
    bedtools=2.25.0-1

RUN wget https://github.com/samtools/htslib/releases/download/1.3.2/htslib-1.3.2.tar.bz2 && \
    tar jxvf htslib-1.3.2.tar.bz2 && \
    cd htslib-1.3.2 && \
    make && \
    make install


RUN pip install --upgrade pip
RUN pip install --upgrade setuptools

RUN pip install pysam==0.11.2.1 numpy==1.13.3

RUN wget https://github.com/friend1ws/annot_utils/archive/v0.1.0.tar.gz && \
    tar xzvf v0.1.0.tar.gz && \
    cd annot_utils-0.1.0/resource && \
    bash prep_data.sh && \
    cd .. && \
    python setup.py build install

RUN wget https://github.com/friend1ws/junc_utils/archive/v0.2.0.tar.gz && \
    tar xzvf v0.2.0.tar.gz && \
    rm -rf v0.2.0.tar.gz && \
    cd junc_utils-0.2.0 && \
    python setup.py build install

RUN wget https://github.com/friend1ws/intron_retention_utils/archive/v0.3.0.tar.gz && \
    tar xzvf v0.3.0.tar.gz && \
    cd intron_retention_utils-0.3.0 && \
    python setup.py build install

RUN wget https://github.com/friend1ws/SAVNet/archive/v0.2.0.tar.gz && \
    tar xzvf v0.2.0.tar.gz && \
    cd SAVNet-0.2.0 && \
    python setup.py build install


