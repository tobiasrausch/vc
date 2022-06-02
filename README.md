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

In this practical we will analyze germline and somatic structural variants (SVs) of a chromothripsis sample from a recent [cancer study](https://www.ncbi.nlm.nih.gov/pubmed/22265402). The anonymized data was filtered for chr2 to speed up all subsequent analysis. The tumor genome alignment file is named `tumor.bam` and the control genome alignment file is named `control.bam`.

### Structural variant alignment quality control

Before each discovery of structural variants, you should assess the quality of the data,
as, for example, paired-end methods are hampered by skewed insert size distributions, read-depth methods by uneven coverage, and split-read methods by high sequencing error rates. Common quality criteria are e.g. the percentage of reads mapped, number of singletons and duplicates, number of properly paired reads and the shape of the insert size and coverage distributions. 
[Picard](http://broadinstitute.github.io/picard/), [SAMtools](http://www.htslib.org), [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [Alfred](https://github.com/tobiasrausch/alfred) are commonly used quality control tools that compute some of these alignment statistics as shown below for the tumor sample.

```bash
cd data/sv/
samtools flagstat tumor.bam
alfred qc -r chr2.fa -o qc.tsv.gz -j qc.json.gz tumor.bam
zcat qc.tsv.gz | grep ^ME | datamash transpose | column -t
```

Instead of parsing the tab-delimited file, you can also upload the JSON file `qc.json.gz` to the [Alfred web application](https://www.gear-genomics.com/alfred/).


As you can see from the QC results, the data has been downsampled to 7x coverage to speed up all analyses.
This implies that some SVs will have only weak support due to low coverage. In terms of QC interpretation, there are some general things to look out for, such as mapping percentages below 70%, >20% duplicates, or multiple peaks in the insert size distribution. Notice that many alignment statistics vary greatly depending on the protocol used, so it's usually best to compare several different sequencing runs from the same protocol (DNA-seq, RNA-seq, ChIP-seq, paired-end, single-end, or mate-pair) to highlight outliers.
 
#### Exercises

* What is the median coverage of the data set?
* Given the insert size distribution, what would be a suitable cutoff to define deletion supporting paired-ends?


### Germline Structural Variants

Before we dive into SV calling, let's get an idea of what structural variants (SVs) look like in short-read sequencing data. For this I have prepared a [BED file](https://bedtools.readthedocs.io/) with some "simple" germline structure variants like deletions and some more complex examples.

```bash
cat svs.bed
```

Using [IGV](http://software.broadinstitute.org/software/igv/) we can browse these SVs interactively.

```bash
igv -g chr2.fa
```

Once IGV has started use 'File' and 'Load from File' to load the `tumor.bam` and `control.bam` alignment file. Then import the `svs.bed` file from your working directory using 'Regions' and 'Import Regions'.
You can then easily navigate to the structural variants with 'Regions' and 'Region Navigator'.
Select a structural vaariant in the region navigator and click 'View', which will center the IGV alignment view on the selected structural variant.
You can zoom in and out using the '+' and '-' signs in the toolbar at the top.
To highlight the abnormal paired-ends please right click in IGV on the BAM file and activate 'View as pairs'. In the same menu, please open 'Color alignments by' and then switch to "pair orientation' for inversions and duplications. For deletions, you want to color the alignments by "insert size". 

### Plotting structural variants

IGV is excellent for interactive browsing but for large numbers of SVs you can use command-line tools such as [wally](https://github.com/tobiasrausch/wally) to plot multiple SVs in batch.

```bash
wally region -R svs.bed -cp -g chr2.fa tumor.bam control.bam
```

### Complex structural variants

Even in germline genomes we can observe complex SVs and two example regions are in the `svs.bed` file.

```bash
cat svs.bed | grep "complex"
```

As part of the [1000 Genomes SV consortium](https://www.nature.com/articles/nature15394) we validated some of the above complex SVs using PacBio. The reads are in a separate FASTA file called `pacbio.sv1.fa` and `pacbio.sv2.fa`. We need the subsequence of the reference to create a pairwise dotplot of the PacBio read against the reference. [SAMtools](http://www.htslib.org) is a convenient tool to extract such subsequences of a FASTA file.

```bash
samtools faidx chr2.fa chr2:18905691-18907969 > sv1.fa
samtools faidx chr2.fa chr2:96210505-96212783 > sv2.fa
```

Please align the above genomic reference subsequences `sv1.fa` and `sv2.fa` against the respective PacBio read `pacbio.sv1.fa` and `pacbio.sv2.fa` using [Maze](https://www.gear-genomics.com/maze/) available on [gear-genomics.com](https://www.gear-genomics.com/). 

***Exercises***

* What kind of SV is present in the region chr2:18905691-18907969 ?
* What kind of SV is present in the region chr2:96210505-96212783 ?


### Somatic structural variant calling

[Delly](https://github.com/dellytools/delly) is a method for detecting structural variants. 
Using the tumor and normal genome alignment, delly calculates structural variants and outputs them as a BCF file, the binary encoding of [VCF](https://samtools.github.io/hts-specs).
You can also provide a text file with regions to exclude from structural variant analysis. Delly's default exclusion map removes the telomeric and centromeric regions of all human chromosomes because these repetitive regions cannot be analyzed with short-read data.

```bash
delly call -q 20 -g chr2.fa -x hg19.ex -o sv.bcf tumor.bam control.bam
```

#### VCF encoding of structural variants

VCF was originally designed for short variants and that's why all SV callers heavily use the VCF INFO fields to encode additional information about the SV such as the structural variant end (INFO:END) and the SV type (INFO:SVTYPE). You can look at the header of the BCF file using grep where '-A 2' includes the first two structural variant records after the header in the file:

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
