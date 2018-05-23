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
    libcanberra-gtk-module \
    gthumb \
    xpdf \
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
    && git clone https://github.com/tobiasrausch/alfred.git \
    && cd /opt/alfred/ \
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

# Create data directory
RUN mkdir -p /data && chmod a+rwx /data

# add user
RUN groupadd -r -g 1000 ubuntu && useradd -r -g ubuntu -u 1000 -m ubuntu
USER ubuntu
RUN cd /data \
    && wget ftp://ftp-exchange.embl-heidelberg.de/pub/exchange/rausch/outgoing/rd.tar.gz \
    && tar -xzf rd.tar.gz \
    && rm rd.tar.gz
RUN cd /data \
    && wget ftp://ftp-exchange.embl-heidelberg.de/pub/exchange/rausch/outgoing/sv.tar.gz \
    && tar -xzf sv.tar.gz \
    && rm sv.tar.gz

# by default /bin/bash is executed
CMD ["/bin/bash"]
