# Structural variant calling tutorial using [delly](https://github.com/dellytools/delly).

## Installation

`make all`

## Load the conda environment with all required tools

```bash
export PATH=`pwd`/conda/bin:${PATH}
source activate base
```

## SV Calling

### Discovery of chromothripsis in a cancer sample

In this practical we will analyze germline and somatic structural variants (SVs) of a Chromothripsis sample from this [study](https://www.ncbi.nlm.nih.gov/pubmed/22265402). The anonomyzed data has been filtered for chr2 only to speed up all analyses hereafter. The alignment file of the tumor genome is called `tumor.bam` and the alignment file of the control genome is called `control.bam`.

## Structural Variant Alignment QC

Paired-end methods can be affected by a skewed insert size distribution, read-depth methods by non-uniform coverage and split-read methods suffer from high sequencing error rates that cause mis-mappings. Prior to any structural variant discovery you should therefore evaluate the quality of the data such as the percentage of mapped reads, singletons, duplicates, properly paired reads and the insert size & coverage distributions. [Picard](http://broadinstitute.github.io/picard/), [SAMtools](http://www.htslib.org), [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [Alfred](https://github.com/tobiasrausch/alfred) compute some of these alignment statistics as shown below for the tumor sample.

```bash
cd data/sv/
samtools flagstat tumor.bam
alfred qc -r chr2.fa -o qc.tsv.gz -j qc.json.gz tumor.bam
zcat qc.tsv.gz | grep ^ME | datamash transpose | column -t
```

Instead of parsing the tab-delimited file, you can also upload the JSON file `qc.json.gz` to the [Alfred web interface](https://www.gear-genomics.com/alfred/).

As you can see in the QC results, I heavily downsampled the data to 7x coverage to speed-up all analyses but this implies that some SVs will have low support. Regarding the QC interpretation, there are some general things to watch out for such as mapping percentages below 70%, >20% duplicates or multiple peaks in the insert size distribution. Be aware that many alignment statistics vary largely by protocol and hence, it's usually best to compare multiple different sequencing runs using the same protocol (DNA-seq, RNA-seq, ChIP-seq, paired-end, single-end or mate-pair) against each other, which then highlights the outliers.

***Exercises***

* What is the median coverage of the data set?
* Given the insert size distribution, what would be a suitable cutoff to define deletion supporting paired-ends?


## Germline Structural Variants

Before diving into SV calling, let's get an idea of how structural variants (SVs) look like in short read, next-generation sequencing data. I have prepared a simple [BED](https://bedtools.readthedocs.io/) file with some "simple" germline structural variants such as deletions and some more complex examples.

```bash
cat svs.bed
```

Using [IGV](http://software.broadinstitute.org/software/igv/) we can then browse these SVs interactively.

```bash
./igv.sh -g chr2.fa
```

Once IGV has started use 'File' and 'Load from File' to load the `tumor.bam` and `control.bam` alignment file. Then import the file `svs.bed` from your working directory using 'Regions' and 'Import Regions'. The structural variants can then be browsed easily using 'Regions' and 'Region Navigator'. Select a structural variant in the Region Navigator and click 'View', which will center the IGV alignment view on the selected structural variant. It's usually best to zoom out once then by clicking on the '-' sign in the toolbar at the top, so you can view all supporting abnormal paired-ends. To highlight the abnormal paired-ends please right click in IGV on the BAM file and activate 'View as pairs'. In the same menu, please open 'Color alignments by' and then switch to "pair orientation' for inversions and duplications. For deletions, you want to color the alignments by "insert size". 

Most [1000 Genomes SVs](https://www.nature.com/articles/nature15394) are deletions because these are easier to detect in low coverage sequencing data. Deletions cause a drop in read-depth and can be detected by spanning paired-ends of abnormally large insert size. In addition, split-reads have a prefix alignment before the deletion and a suffix alignment after the deletion. Other SV types such as inversions are much harder to detect because these are balanced rearrangements that do not cause a read-depth change. Besides the simple SV types (deletions, duplications, inversions) one can also find more complex SVs in the germline. Three example regions for that are in the `svs.bed` file.

```bash
cat svs.bed | grep "complex"
```

As part of the SV consortium we validated some of the above complex SVs using PacBio. The reads are in a separate FASTA file called `pacbio.sv1.fa` and `pacbio.sv2.fa`. We need the subsequence of the reference to create a pairwise dotplot of the PacBio read against the reference. [SAMtools](http://www.htslib.org) is a convenient tool to extract such subsequences of a FASTA file.

```bash
head pacbio.sv1.fa
samtools faidx chr2.fa chr2:18905691-18907969
head pacbio.sv2.fa
samtools faidx chr2.fa chr2:96210505-96212783
```

Please align the above genomic reference subsequences against the respective PacBio read using [Maze](https://www.gear-genomics.com/alfred/) available on [gear-genomics.com](https://www.gear-genomics.com/). 

***Exercises***

* What kind of SV is present in the region chr2:18905691-18907969 ?
* What kind of SV is present in the region chr2:96210505-96212783 ?


## Structural Variant Calling

As a first step we will discover structural variants using [Delly](https://github.com/dellytools/delly). [Delly](https://www.ncbi.nlm.nih.gov/pubmed/22962449) calls structural variants jointly on the tumor and normal genome and outputs a [BCF](https://samtools.github.io/hts-specs) file, the binary encoding of [VCF](https://samtools.github.io/hts-specs). You can also provide a text file with regions to exclude from the analysis of structural variants. The default exclude map of Delly removes the telomeric and centromeric regions of all human chromosomes since these regions cannot be accurately analyzed with short-read data.

```bash
cd /data/sv
delly call -n -q 20 -g chr2.fa -x hg19.ex -o sv.bcf tumor.bam control.bam
```

VCF was originally designed for short variants and that's why all SV callers heavily use the INFO fields to encode additional information about the SV such as the structural variant end (INFO:END) and the SV type (INFO:SVTYPE). You can look at the header of the BCF file using grep, '-A 2' includes the first two structural variant records in the file:

```bash
bcftools view sv.bcf | grep "^#" -A 2
```

[Delly](https://github.com/dellytools/delly) uses the VCF:INFO fields for structural variant site information, such as how confident the structural variant prediction is and how accurate the breakpoints are. The genotype fields contain the actual sample genotype, its genotype quality and genotype likelihoods and various count fields for the variant and reference supporting reads and spanning pairs. If you browse through the VCF file you will notice that a subset of the Delly structural variant predictions have been refined using split-reads. These precise variants are flagged in the VCF info field with the tag 'PRECISE', all others are listed as 'IMPRECISE'. Please note that this BCF file contains germline and somatic structural variants but also false positives caused by repeat-induced mis-mappings or incomplete reference sequences. [SVprops](https://github.com/dellytools/svprops) is a simple program that converts Delly's BCF output to a tab-delimited SV site list.

```bash
svprops sv.bcf | head
```

Every record of the BCF file is converted to one tab-delimited row with all kinds of summary statistics. [SVprops](https://github.com/dellytools/svprops) also provides a 'column-view' listing summary statistics for all samples present in the BCF file.

```bash
sampleprops sv.bcf
```

This initial SV calling cannot differentiate somatic and germline structural variants. For instance, the 2 complex variants we looked at before are still present in the Delly output. A proximal duplication causes 2 paired-end signatures (deletion-type and duplication-type):

```bash
bcftools view sv.bcf chr2:18905691-18907969 | awk '$2>=18905691 && $2<=18907969'
```

***Exercises***

* What is the fraction of deletions that has been called precisely (at single nucleotide resolution) by Delly?
* Is the germline inverted duplication present in Delly's output files?
* What type of SVs delineate a proximal inverted duplication?
* Some of the SVs have nucleotide resolution and the alternative haplotype is present in INFO:CONSENSUS. [Blat](https://genome.ucsc.edu/cgi-bin/hgBlat) the consensus sequence of some of these deletions. What genomic element has been deleted for DEL00001919 (bcftools view sv.bcf | grep "DEL00001919")? Is it a known variant?


## Somatic Filtering

Delly ships with a basic somatic filtering subcommand that uses the matched control and possibly additional control samples from unrelated individuals. The somatic filtering requires a sample file listing tumor and control sample names from the VCF file.

```bash
cat spl.tsv
```

There are many parameters available to tune the somatic structural variant filtering. Below we require a minimum variant allele frequency of 25%, no support in the matched normal and an overall confident structural variant site prediction with the VCF filter field being equal to PASS.

```bash
delly filter -p -f somatic -o somatic.bcf -a 0.25 -s spl.tsv sv.bcf
sampleprops somatic.bcf
```

***Exercises***

* What is the average size of the somatic SVs?
* Is there any SV type that is enriched among the somatic SVs?


## Structural Variant Visualization using IGV

IGV is excellent for inspecting small variants but for these long-range complex re-arrangements its less powerful. You cannot view the entire event in IGV but you can still visualize the breakpoint regions separately, denoted as left breakpoint (L) and right breakpoint (R) below.

```bash
svprops somatic.bcf | tail -n +2 | awk '{print $1"\t"($2-500)"\t"($2+500)"\t"$5"L";}' > somatic.bp.bed
svprops somatic.bcf | tail -n +2 | awk '{print $3"\t"($4-500)"\t"($4+500)"\t"$5"R";}' >> somatic.bp.bed
sort -k4,4 somatic.bp.bed
```

We can then use [IGV](http://software.broadinstitute.org/software/igv/) to browse the variants interactively using the chr2 reference sequence.

```bash
./igv.sh -g chr2.fa
```

As previously, once IGV has started use 'File' and 'Load from File' to load the tumor and normal bam alignment file. Then import 'somatic.bp.bed' from your working directory using 'Regions' and 'Import Regions'. The somatic structural variants can then be browsed easily using 'Regions' and 'Region Navigator' and you probably want to switch on 'View as pairs' and 'Color alignments by pair orientation'. You may also want to play around with the split screen view, where you can right click a read and select 'View mate region in split screen'.


## Visualizing Complex Rearrangements

To reveal any complex re-arrangement patterns it's worthwhile to create a read-depth plot, overlay the long-range structural variants and visualize the B-allele frequency of germline SNPs to highlight loss-of-heterozygosity (LOH) events.

```bash
coral call -m chr2.map.fa -g chr2.fa -v chr2.snps.bcf -l control.bam -s id -o out tumor.bam
zcat out.adaptive.cov.gz | head
```

A simple read-depth plot can then be generated using

```bash
Rscript cnBafSV.R out.adaptive.cov.gz
open cov.png
```

Next, we will try to visualize the somatic SV calls on top of the log2 read-depth ratio plot.

```bash
svprops somatic.bcf > svs.tsv
head svs.tsv
Rscript cnBafSV.R out.adaptive.cov.gz svs.tsv
open cov.png
```

Lastly, we can augment the plots with the allelic depth of SNPs. For the amplifications, we would expect the B-allele frequency to branch, similar to LOH patterns for somatic deletions.

```bash
zcat out.baf.gz | head
Rscript cnBafSV.R out.adaptive.cov.gz svs.tsv out.baf.gz
open cov.png
```
