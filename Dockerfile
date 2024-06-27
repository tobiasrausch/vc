# use the ubuntu base image
FROM ubuntu:22.04

# maintainer
MAINTAINER Tobias Rausch rausch@embl.de

# install required apt packages
RUN apt-get update && apt-get install -y \
    autoconf \
    build-essential \
    curl \
    g++ \
    git \
    cmake \
    gthumb \
    libcurl4-gnutls-dev \
    libboost-date-time-dev \
    libboost-program-options-dev \
    libboost-system-dev \
    libboost-filesystem-dev \
    libboost-iostreams-dev \
    libbz2-dev \
    libdeflate-dev \
    libhdf5-dev \
    libncurses-dev \
    liblzma-dev \    
    python3 \
    datamash \
    gthumb \
    pkg-config \
    wget \
    zlib1g-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# install tutorial tools
RUN cd /opt \
    && git clone --recursive https://github.com/tobiasrausch/vc.git \
    && cd /opt/vc \
    && make all

# by default /bin/bash is executed
CMD ["/bin/bash"]
