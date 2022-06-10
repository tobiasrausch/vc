# use the ubuntu base image
FROM ubuntu:18.04

# maintainer
MAINTAINER Tobias Rausch rausch@embl.de

# install required apt packages
RUN apt-get update && apt-get install -yq \
    build-essential \
    g++ \
    git \
    cmake \
    python \
    datamash \
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
