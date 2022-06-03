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
zcat qc.tsv.gz | grep ^ME | datamash transpose
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


### Delly structural variant calling

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

[Delly](https://github.com/dellytools/delly) uses the VCF:INFO fields for structural variant site information, such as how confident the structural variant prediction is and how accurate the breakpoints are. The genotype fields contain the actual sample genotype, its genotype quality and genotype likelihoods and various count fields for the variant and reference supporting reads and spanning pairs. If you browse through the VCF file you will notice that a subset of the Delly structural variant predictions have been refined using split-reads. These precise variants are flagged in the VCF info field with the tag 'PRECISE', all others are listed as 'IMPRECISE'. Please note that this BCF file contains germline and somatic structural variants but also false positives caused by repeat-induced mis-mappings or incomplete reference sequences.

#### Querying VCF files

[Bcftools](https://github.com/samtools/bcftools) offers many possibilities to query and reformat SV calls. For instance, to output a table with the chromosome, start, end, identifier and genotype of each SV we can use:

```bash
bcftools query -f "%CHROM\t%POS\t%INFO/END\t%ID[\t%GT]\n" sv.bcf | head
```

This initial SV calling cannot differentiate somatic and germline structural variants. For instance, the 2 complex variants we looked at before are still present in the delly output. A proximal duplication causes 2 paired-end signatures (deletion-type and duplication-type):

```bash
bcftools view sv.bcf chr2:18905691-18907969 | awk '$2>=18905691 && $2<=18907969'
```

#### Exercises

* What is the fraction of deletions that has been called precisely (at single nucleotide resolution) by Delly?
* Is the other complex structural variant still present in Delly's output file?

### Somatic structural variant filtering

Delly's somatic filtering requires a sample file listing tumor and control sample names from the VCF file.

```bash
cat spl.tsv
```

There are many parameters available to tune the somatic structural variant filtering. Below we require a minimum variant allele frequency of 25%, no support in the matched normal and an overall confident structural variant site prediction with the VCF filter field being equal to PASS.

```bash
delly filter -p -f somatic -o somatic.bcf -a 0.25 -s spl.tsv sv.bcf
```

As expected, the somatic SVs have a homozygous reference genotype in the control sample.

```bash
bcftools query -f "%CHROM\t%POS\t%INFO/END\t%ID[\t%GT]\n" somatic.bcf
```

#### Exercises

* What is the average size of the somatic SVs?
* Is there any SV type that is enriched among the somatic SVs?


### Complex structural variant visualization

IGV and [wally](https://github.com/tobiasrausch/wally) support so-called split views to visualize the breakpoints of long-range SVs greater than 10,000kbp.

```bash
bcftools query -f "%CHROM\t%POS\t%INFO/END\t%ID\n" somatic.bcf | awk '$3-$2>10000 {print $1"\t"($2-500)"\t"($2+500)"\t"$4"L\n"$1"\t"($3-500)"\t"($3+500)"\t"$4"R";}' > somatic.bp.bed
wally region -R somatic.bp.bed -s 2 -cp -g chr2.fa tumor.bam control.bam
```

To reveal any higher-order SV class such as chromothripsis we need to integrate read-depth with structural variant predictions. Let's first create a simple read-depth plot.

```bash
delly cnv -u -z 10000 -o cnv.bcf -c cnv.cov.gz -g chr2.fa -m chr2.map.fa tumor.bam 
Rscript cnBafSV.R cnv.cov.gz
```

Now we can overlay the somatic structural variants on top of the read-depth information.

```bash
bcftools query -f "%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%ID\n" somatic.bcf > svs.tsv
Rscript cnBafSV.R cnv.cov.gz svs.tsv
```

Lastly, we can augment the plots with the allelic depth of SNPs. For amplifications, we would expect the variant allele frequencies to deviate from the expected 50% for heterozygous variants. To speed up things, we only compute SNPs in the region 1 - 50Mbp.

```bash
bcftools mpileup -d 50 -r chr2:1-50000000 -a FORMAT/AD -f chr2.fa tumor.bam | bcftools call -mv -Ob -o calls.bcf
bcftools view -g het -m2 -M2 -v snps calls.bcf | bcftools query -f "%POS,[%AD\n]" - | awk 'BEGIN {FS=","} $2+$3>10 {print $1"\t"$3/($2+$3);}' > baf.tsv
Rscript cnBafSV.R cnv.cov.gz svs.tsv baf.tsv 
```
