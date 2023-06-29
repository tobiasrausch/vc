SHELL := /bin/bash

# Targets
TARGETS = .mamba .tools .rstats .pcks .check
PBASE=$(shell pwd)

all: ${TARGETS}

.mamba:
	curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(shell uname)-$(shell uname -m).sh" && bash Mambaforge-$(shell uname)-$(shell uname -m).sh -b -p mamba && rm "Mambaforge-$(shell uname)-$(shell uname -m).sh" && touch .mamba

.tools: .mamba
	export PATH=${PBASE}/mamba/bin:${PATH} && mamba install -y -c conda-forge -c bioconda datamash samtools bcftools bedtools htslib bwa delly alfred freebayes igv wally minimap2 && touch .tools

.rstats: .mamba .tools
	export PATH=${PBASE}/mamba/bin:${PATH} && mamba install -y -c conda-forge -c bioconda bioconductor-genomicfeatures r-ggplot2 r-reshape2 bioconductor-dnacopy && touch .rstats

.pcks: .mamba .tools .rstats
	export PATH=${PBASE}/mamba/bin:${PATH} && mamba install -y -c conda-forge -c bioconda cyvcf2 numpy pysam && pip install gdown && touch .pcks

.check: .mamba .tools .rstats .pcks
	export PATH=${PBASE}/mamba/bin:${PATH} && delly --version && touch .check

download: .mamba .tools .rstats .pcks
	export PATH=${PBASE}/mamba/bin:${PATH} && cd data/ && gdown ${FILE} && tar -xzf sv.tar.gz && rm sv.tar.gz

clean:
	rm -rf $(TARGETS) $(TARGETS:=.o) mamba/
