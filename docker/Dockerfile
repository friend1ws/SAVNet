FROM ubuntu:16.04
MAINTAINER Yuichi Shiraishi <friend1ws@gmail.com> 

RUN apt-get update && apt-get install -y \
    git \
    wget \
    bzip2 \
    make \
    gcc \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    python \
    python-pip

RUN wget https://github.com/samtools/htslib/releases/download/1.7/htslib-1.7.tar.bz2 && \
    tar -jxvf htslib-1.7.tar.bz2 && \
    cd htslib-1.7 && make && make install 

RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.27.0/bedtools-2.27.0.tar.gz && \
    tar -zxvf bedtools-2.27.0.tar.gz && \
    cd bedtools2 && make && make install

RUN pip install --upgrade pip
RUN pip install --upgrade setuptools

RUN pip install pysam==0.13
RUN pip install annot_utils==0.2.1
RUN pip install junc_utils==0.4.1
RUN pip install intron_retention_utils==0.5.1 
RUN pip install chimera_utils==0.5.1
RUN pip install savnet==0.3.2

