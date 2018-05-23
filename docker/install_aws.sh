#!/bin/bash

# Run this as root on Amazon AWS to install all dependencies

# update package repository
echo "deb http://cran.uni-muenster.de/bin/linux/ubuntu trusty/" >> /etc/apt/sources.list
gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9
gpg -a --export E084DAB9 | sudo apt-key add -

# install required packages
sudo apt-get update && sudo apt-get install -y \
    build-essential \
    g++ \
    git \
    cmake \
    zlib1g-dev \
    ant \
    libbz2-dev \
    liblzma-dev \
    libboost-date-time-dev \
    libboost-program-options-dev \
    libboost-system-dev \
    libboost-filesystem-dev \
    libboost-iostreams-dev \
    python \
    python-dev \
    python-pip \
    gfortran \
    libncurses-dev \
    datamash \
    libncurses5-dev \
    hdf5-tools \
    libhdf5-dev \
    wget \
    openjdk-8-jdk \
    x11-apps \
    emacs \
    vcftools \
    r-base \
    feh \
    firefox \
    libcurl4-openssl-dev \
    libxml2-dev \
    libcanberra-gtk-module \
    gthumb \
    xpdf \
    && apt-get clean

# set environment
export BOOST_ROOT=/usr
export SEQTK_ROOT=/opt/htslib
export IDIR=/opt
export bindir=/usr/local/bin/

# install delly & tutorial tools
cd ${IDIR} \
    && wget 'https://github.com/samtools/htslib/releases/download/1.6/htslib-1.6.tar.bz2' \
    && tar -xjf htslib-1.6.tar.bz2 \
    && cd htslib-1.6 \
    && make \
    && make install
cd ${IDIR} \
    && wget 'https://github.com/samtools/bcftools/releases/download/1.6/bcftools-1.6.tar.bz2' \
    && tar -xjf bcftools-1.6.tar.bz2 \
    && cd bcftools-1.6 \
    && make all \
    && make install
cd ${IDIR} \
    && wget 'https://github.com/samtools/samtools/releases/download/1.6/samtools-1.6.tar.bz2' \
    && tar -xjf samtools-1.6.tar.bz2 \
    && cd samtools-1.6 \
    && make all \
    && make install
cd ${IDIR} \
    && wget https://github.com/igvteam/igv/archive/v.2.3.77.tar.gz \
    && tar -xzf v.2.3.77.tar.gz \
    && cd igv-v.2.3.77/ \
    && ant
cd ${IDIR} \
    && git clone https://github.com/dellytools/delly.git \
    && cd /opt/delly/ \
    && make all \
    && make install
cd ${IDIR} \
    && git clone https://github.com/dellytools/svprops.git \
    && cd /opt/svprops/ \
    && make all \
    && make install
cd ${IDIR} \
    && git clone https://github.com/tobiasrausch/alfred.git \
    && cd /opt/alfred/ \
    && make all \
    && make install
cd ${IDIR} \
    && wget 'https://github.com/arq5x/bedtools2/releases/download/v2.26.0/bedtools-2.26.0.tar.gz' \
    && tar -xzf bedtools-2.26.0.tar.gz \
    && cd bedtools2/ \
    && make \
    && make install
cd ${IDIR} \
    && wget 'https://github.com/lh3/bwa/releases/download/v0.7.16/bwa-0.7.16a.tar.bz2' \
    && tar -xjf bwa-0.7.16a.tar.bz2 \
    && cd bwa-0.7.16a/ \
    && make \
    && install -p /opt/bwa-0.7.16a/bwa /usr/local/bin/
cd ${IDIR} \
    && wget 'https://github.com/gt1/libmaus2/archive/2.0.403-release-20171019153510.tar.gz' \
    && tar -xzf 2.0.403-release-20171019153510.tar.gz \
    && cd libmaus2* \
    && ./configure \
    && make \
    && make install
cd ${IDIR} \
    && wget 'https://github.com/gt1/biobambam2/archive/2.0.79-release-20171006114010.tar.gz' \
    && tar -xzf 2.0.79-release-20171006114010.tar.gz \
    && cd biobambam2* \
    && ./configure \
    && make \
    && make install
cd ${IDIR} \
    && git clone --recursive git://github.com/ekg/freebayes.git \
    && cd freebayes \
    && git checkout tags/v1.1.0 \
    && git submodule update --init --recursive \
    && make \
    && make install

# install python variant filtering dependencies
pip install -r requirements.txt

# install bioconductor packages
Rscript pkginstall.R

# Create data directory
mkdir -p /data && chmod a+rwx /data
