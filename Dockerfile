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
    libbz2-dev \
    liblzma-dev \
    libboost-date-time-dev \
    libboost-program-options-dev \
    libboost-system-dev \
    libboost-filesystem-dev \
    libboost-iostreams-dev \
    python \
    datamash \
    wget \
    openjdk-8-jdk \
    x11-apps \
    emacs \
    r-base \
    firefox \
    libcurl4-openssl-dev \
    libxml2-dev \
    gthumb \
    pkg-config \
    && apt-get clean

# install tutorial tools
RUN cd /opt \
    && git clone --recursive https://github.com/tobiasrausch/vc.git \
    && cd /opt/vc \
    && make all

# by default /bin/bash is executed
CMD ["/bin/bash"]
