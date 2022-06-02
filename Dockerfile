# use the ubuntu base image
FROM ubuntu:16.04

# maintainer
MAINTAINER Tobias Rausch rausch@embl.de

# install required packages
RUN apt-get update && apt-get install -y \
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
    imagemagick \
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
    xvfb \
    pkg-config \
    && apt-get clean

# set environment
ENV BOOST_ROOT /usr
ENV SEQTK_ROOT /opt/htslib
ENV EBROOTHTSLIB /opt/htslib
ENV bindir /usr/local/bin/

# install delly & tutorial tools
RUN cd /opt \
    && git clone https://github.com/samtools/htslib.git \
    && cd /opt/htslib \
    && make \
    && make lib-static \
    && make install
RUN cd /opt \
    && git clone https://github.com/samtools/bcftools.git \
    && cd /opt/bcftools \
    && make all \
    && make install
RUN cd /opt \
    && git clone https://github.com/samtools/samtools.git \
    && cd /opt/samtools \
    && make all \
    && make install
RUN cd /opt \
    && wget 'https://github.com/arq5x/bedtools2/releases/download/v2.26.0/bedtools-2.26.0.tar.gz' \
    && tar -xzf bedtools-2.26.0.tar.gz \
    && rm bedtools-2.26.0.tar.gz \
    && cd bedtools2/ \
    && make \
    && make install
RUN cd /opt \
    && wget 'https://github.com/lh3/bwa/releases/download/v0.7.16/bwa-0.7.16a.tar.bz2' \
    && tar -xjf bwa-0.7.16a.tar.bz2 \
    && rm bwa-0.7.16a.tar.bz2 \
    && cd bwa-0.7.16a/ \
    && make \
    && install -p /opt/bwa-0.7.16a/bwa /usr/local/bin/
RUN cd /opt \
    && wget 'https://github.com/gt1/libmaus2/archive/2.0.403-release-20171019153510.tar.gz' \
    && tar -xzf 2.0.403-release-20171019153510.tar.gz \
    && rm 2.0.403-release-20171019153510.tar.gz \
    && cd libmaus2* \
    && ./configure \
    && make \
    && make install
RUN cd /opt \
    && wget 'https://github.com/gt1/biobambam2/archive/2.0.79-release-20171006114010.tar.gz' \
    && tar -xzf 2.0.79-release-20171006114010.tar.gz \
    && rm 2.0.79-release-20171006114010.tar.gz \
    && cd biobambam2* \
    && ./configure \
    && make \
    && make install
RUN cd /opt \
    && git clone --recursive git://github.com/ekg/freebayes.git \
    && cd freebayes \
    && git checkout tags/v1.1.0 \
    && git submodule update --init --recursive \
    && make \
    && make install
RUN cd /opt \
    && git clone https://github.com/tobiasrausch/alfred.git \
    && cd /opt/alfred/ \
    && make all \
    && make install
RUN cd /opt \
    && git clone https://github.com/tobiasrausch/coral.git \
    && cd /opt/coral/ \
    && make all \
    && make install
RUN cd /opt \
    && git clone https://github.com/dellytools/delly.git \
    && cd /opt/delly/ \
    && make all \
    && make install
RUN cd /opt \
    && git clone https://github.com/dellytools/svprops.git \
    && cd /opt/svprops/ \
    && make all \
    && make install
RUN cd /opt \
    && wget https://github.com/igvteam/igv/archive/v.2.3.77.tar.gz \
    && tar -xzf v.2.3.77.tar.gz \
    && rm v.2.3.77.tar.gz \
    && cd igv-v.2.3.77/ \
    && ant
RUN cd /opt \
    && git clone https://github.com/tobiasrausch/variant-calling.git

# install python variant filtering dependencies
RUN pip install -r /opt/variant-calling/docker/requirements.txt

# install bioconductor packages
RUN Rscript /opt/variant-calling/docker/pkginstall.R

# ldconfig
RUN ldconfig

# Create data directory
RUN mkdir -p /data && chmod a+rwx /data

# add user
RUN groupadd -r -g 1000 ubuntu && useradd -r -g ubuntu -u 1000 -m ubuntu
USER ubuntu

# Set user environment
ENV LD_LIBRARY_PATH /usr/local/lib
RUN echo alias open=feh >> ~/.bashrc

# copy data as local user to /data directory
# docker run -i -t -v /opt/dev/variant-calling/data/:/data/tmp variant-calling /bin/bash

# by default /bin/bash is executed
CMD ["/bin/bash"]
