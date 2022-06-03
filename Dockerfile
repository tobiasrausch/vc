# use the ubuntu base image
FROM ubuntu:20.04

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

# install tutorial tools
RUN cd /opt \
    && git clone --recursive https://github.com/tobiasrausch/vc.git \
    && cd /opt/vc \
    && make all

# by default /bin/bash is executed
CMD ["/bin/bash"]
