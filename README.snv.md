# SNV Calling

## Introduction

### Disovering a Causative Mutation in a Rare Disease Patient

This practical provides an introduction into variant discovery and genotyping. We will cover single-nucleotide variants, short insertions and deletions (InDels) and large structural variants. All data of this practical has been anonomyzed and subsampled to speed up the analyses. We will start with a [rare disease](https://www.ncbi.nlm.nih.gov/pubmed/23999272) case that we analyzed in 2012. This infant of consanguineous parents suffered from severe combined immunodeficiency and the details are described in this [publication](https://www.ncbi.nlm.nih.gov/pubmed/23561803).

## Reference Indices

We will first map the data to the Human reference genome using [BWA](https://github.com/lh3/bwa). To speed up the mapping the reference genome needs to be indexed. BWA uses an FM-Index which is built around the [Burrows-Wheeler transform](https://de.wikipedia.org/wiki/Burrows-Wheeler-Transformation).

```bash
cd /data/rd/
bwa index chr7.fa
ls -rt1 chr7.fa*
```

It is also useful to build an index of the FASTA reference file using [SAMtools](http://www.htslib.org) to allow a quick extraction of subsequences from the reference genome.

```bash
samtools faidx chr7.fa
```

We can now, for instance, extract 50bp from position 10017.

```bash
samtools faidx chr7.fa chr7:10017-10067
```

[bedtools](http://bedtools.readthedocs.io/en/latest/) can be used to create a simple BED file with the start and end of each chromosome. We could also use the command makewindows to tile each chromosome into 10kbp windows.


```bash
bedtools makewindows -g <(cut -f 1,2 chr7.fa.fai) -n 1 > chr7.bed
bedtools nuc -fi chr7.fa -bed chr7.bed
```

***Exercises***

* What is the length of chr7?
* What is the GC-content of chr7?
* What is the proportion of Ns in chr7?

## Alignment

Once the index has been built we can map the paired-end [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) data against the reference and convert it to [BAM](http://www.htslib.org).

```bash
bwa mem chr7.fa read1.fq.gz read2.fq.gz | samtools view -bT chr7.fa - > rd.bam
```

Using [samtools](http://www.htslib.org) we can have a look at the header and the first few alignment records:

```bash
samtools view -H rd.bam
samtools view rd.bam | head
```

Please familiarize yourself with the [BAM](http://www.htslib.org) format, the required fields present in every bam alignment record are explained below:


| Col   | Field    | Description                              |
|-------|----------|------------------------------------------|
|  1    |   QNAME  |    Query template NAME                   |
|  2    |   FLAG   |    bitwise FLAG                          |
|  3    |   RNAME  |    Reference sequence NAME               |
|  4    |   POS    |    1-based leftmost mapping POSition     |
|  5    |   MAPQ   |    MAPping Quality                       |
|  6    |   CIGAR  |    CIGAR string                          |
|  7    |   RNEXT  |    Ref. name of the mate/next read       |
|  8    |   PNEXT  |    Position of the mate/next read        |
|  9    |   TLEN   |    observed Template LENgth              |
|  10   |   SEQ    |    segment SEQuence                      |
|  11   |   QUAL   |    ASCII of Phred-scaled base QUALity+33 |

The bitwise FLAG can be decoded using the [explain flag tool](https://broadinstitute.github.io/picard/explain-flags.html) from the picard distribution.

We need to sort the alignments and can then also built an index to allow a random extraction of alignments.

```bash
samtools sort -o rd.srt.bam rd.bam
samtools index rd.srt.bam
```

## Mark Duplicates and Alignment Quality Control

Unless you are using a PCR-free library, PCR duplicates are common in DNA-sequencing and should be flagged prior to variant calling.

```bash
bammarkduplicates I=rd.srt.bam O=rd.rmdup.bam M=rd.metrics.tsv index=1 rmdup=0
```

SAMtools flagstat computes some basic alignment statistics such as the number of properly paired reads and singletons.

```bash
samtools flagstat rd.rmdup.bam
```

[Alfred](https://github.com/tobiasrausch/alfred) can be used to compute the insert size distribution, the coverage distribution and alignment error rates.

```bash
alfred qc -r chr7.fa rd.rmdup.bam
```

The output file has multiple sections. Most are data matrices for quality control plots but there is also one section with summary alignment metrics.

```bash
zcat qc.tsv.gz | grep ^ME | datamash transpose | column -t
```

To create the QC plots we run the Rscript that is part of [Alfred](https://github.com/tobiasrausch/alfred).

```bash
Rscript /opt/alfred/R/stats.R qc.tsv.gz
```

The output PDF file can then be converted to PNGs.

```bash
convert qc.tsv.gz.pdf qc.png
# Base content distribution
open qc-0.png
# Base quality distribution
open qc-1.png
# Coverage
open qc-4.png
# Insert size
open qc-5.png
```

For calling exonic variants we are primarily interested in the coverage distribution across exons. One could, for instance, use [R Statistics](https://www.r-project.org/) to download exon coordinates for hg19 or download coding regions from UCSC or Ensembl. The file exons.bed.gz contains CCDS coding regions for hg19. We have to subset this bed file to the subsequence of chr7 that we are using in this practical.

```bash
zcat exons.bed.gz | head
bedtools intersect -a <(zcat exons.bed.gz) -b chr7.bed | gzip -c > exons.chr7.bed.gz
zcat exons.chr7.bed.gz | head
```

With the exonic coordinates, [Alfred](https://github.com/tobiasrausch/alfred) can be used to compute the avg. coverage per target region. In our case the targets are CCDS exons but the same method can be used to compute on-target rates for exome capture data sets.

```bash
alfred qc -r chr7.fa -b exons.chr7.bed.gz rd.rmdup.bam
Rscript /opt/alfred/R/stats.R qc.tsv.gz
convert qc.tsv.gz.pdf qc.png
# Exon coverage distribution
open qc-10.png
```



***Exercises***

* What is the median coverage of the data set?
* What is the meaning of the different library layouts (F+, F-, R+, R-)?
* What is the duplicate fraction in the library?
* Would it make sense to sequence this library deeper to achieve 30x coverage?


## Variant Calling

Once the alignment is sorted and duplicates are marked we can run a variant caller such as [FreeBayes](https://github.com/ekg/freebayes) to scan the alignments for differences compared to the reference.

```bash
freebayes --fasta-reference chr7.fa -b rd.rmdup.bam -v snv.vcf
```

Compressing and indexing of the output VCF file will again speed up random access to the file.

```bash
bgzip snv.vcf
tabix snv.vcf.gz
```

The [VCF](https://samtools.github.io/hts-specs) format has multiple header lines starting with the hash # sign. Below the header lines is one record for each variant. The record format is described in the below table:

| Col | Field  | Description         |
|-----|--------|---------------------|
| 1   | CHROM  | Chromosome name |
| 2   | POS    | 1-based position. For an indel, this is the position preceding the indel. |
| 3   | ID     | Variant identifier. Usually the dbSNP rsID. |
| 4   | REF    | Reference sequence at POS involved in the variant. For a SNP, it is a single base. |
| 5   | ALT    | Comma delimited list of alternative sequence(s). |
| 6   | QUAL   | Phred-scaled probability of all samples being homozygous reference. |
| 7   | FILTER | Semicolon delimited list of filters that the variant fails to pass. |
| 8   | INFO   | Semicolon delimited list of variant information. |
| 9   | FORMAT | Colon delimited list of the format of individual genotypes in the following fields. |
| 10+ | Samples| Individual genotype information defined by FORMAT. |

You can look at the header of the VCF file using grep, '-A 1' includes the first variant record in the file:

```bash
bcftools view snv.vcf.gz | grep "^#" -A 1
```

Using [BCFtools](https://samtools.github.io/bcftools/bcftools.html) we can generate some useful summary statistics such as the [transition/transversion ratio](https://en.wikipedia.org/wiki/Transversion).

```bash
bcftools stats snv.vcf.gz | grep "TSTV"
```

***Exercises***

* How many SNPs have been called (hint: bcftools stats, SN tag)?
* How many InDels have been called (hint: bcftools stats, SN tag)?
* How many C>T mutations have been called (hint: bcftools stats, ST tag)?


## Filtering Variants

In most applications researchers use external ground truth data to calibrate a variant calling pipeline. In our case we do not know the ground truth so we will illustrate some filtering options based on summary statistics such as the transition/transversion ratio. In most species, transitions are far more likely than transversions and for humans we would expect a transition/transversion ratio of approximately 2.

```bash
bcftools stats snv.vcf.gz | grep "TSTV"
bcftools filter -i '%QUAL>20' snv.vcf.gz  | bcftools stats | grep "TSTV"
bcftools filter -e '%QUAL<=20 || %QUAL/INFO/AO<=2 || SAF<=2 || SAR<=2' snv.vcf.gz  | bcftools stats | grep "TSTV"
```

Another useful bulk metric is the length of indels in exons because most InDel polymorphisms should be in-frame. If you perform variant calling on a large population cohort with hundreds of samples of different ancestry then [heterozygosity](https://en.wikipedia.org/wiki/Zygosity) is another metric that could be useful. For our single sample case study we move on with a simple threshold based filtering strategy to subset the VCF to exonic variants.

```bash
bcftools filter -O z -o exon.vcf.gz -R <(zcat exons.bed.gz) -e '%QUAL<=20 || %QUAL/INFO/AO<=2 || SAF<=2 || SAR<=2' snv.vcf.gz
bcftools stats exon.vcf.gz | egrep "^SN|TSTV"
```

[SAMtools](http://www.htslib.org) also includes a basic alignment viewer called tview that is useful to spot-check variants in the raw alignment data. For instance, to view the alignment data for the first two exonic variants:


```bash
bcftools view exon.vcf.gz | grep "^#" -A 2
samtools tview -d t -p chr7:299825 rd.rmdup.bam chr7.fa
```


***Exercises***

* Are the first two exonic variants homozygous or heterozygous?
* What is the genotype and the allelic depth for both variants that FreeBayes emits?
* Spot-check some heterozygous variants using samtools tview.
* Plot the InDel length distribution of all called InDels (hint: bcftools stats, IDD tag).


## Variant Annotation

Variant annotation and classification is a challenging process.

* you can use transcript annotations from Ensembl, UCSC or RefSeq
* there is a long list of mutation damaging prediction tools such as PolyPhen, MutationTaster or Sift
* you can annotate variants with allele frequency information from variation archives such as 1000 Genomes, ExAC or gnomAD
* you can check the expression of genes in your studied tissue using GTEx
* you can prioritize mutations in genes that interact with known candidate genes of the disease
* you can categorize known mutations into benign and pathogenic using ClinVar


In the recent years a number of convenient pipelines have been developed that ease the annotation of variants with some of the above information. In this practical we will use [VEP](http://www.ensembl.org/info/docs/tools/vep/index.html) because it can be run directly online. We first dump all SNPs in a VEP compliant format which you can
then copy and paste into the VEP application. Make sure you use the hg19/GRCh37 version available [here](http://grch37.ensembl.org/Homo_sapiens/Tools/VEP).

```bash
bcftools query -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\n" exon.vcf.gz
```

A working variant annotation and classification pipeline can easily reduce an initial call set of several thousands of exonic variants to a handful mutation candidates. In a rare disease setting additional power can be gained by taking advantage of the suspected inheritance model (autosomal recessive, autosomal dominant, etc.). 

```bash
bcftools query -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%GT]\n" exon.vcf.gz
```

***Exercises***

* What would be a useful additional filter in our case given that the index patient has consanguineous parents?
* How could we use variation archives to further filter the list of exonic variants?
* Which annotation features could be used to rank the mutation list for a clinician?


## Variant Validation

Once a putative causative variant has been identified these are usually validated in the index patient using PCR and Sanger sequencing. If a specific inheritance model is suspected the parents are also tested. For the likely causative variant we will design primers using [Primer3Plus](https://www.ncbi.nlm.nih.gov/pubmed/17485472). These primers are locally unique and have the appropriate Tm but they are not necessarily unique in the entire genome. We can use [Silica](https://www.gear-genomics.com/silica), a tool for [In-silico PCR](https://en.wikipedia.org/wiki/In_silico_PCR), to check genome-wide uniqueness. Both methods, Primer3Plus and Silica, are combined in [Verdin](https://www.gear-genomics.com/verdin) to automatically design primers for short variants and large structural variants. The example shows how SNVs, InDels and SVs are encoded and then you can design primers for the likely causative SNV. We do not have the time to run the actual PCR experiment and sequence the breakpoint mutation but the Sanger validation files of the original study are in the data folder.

```bash
ls *.ab1
```

You can analyze these trace files using [Indigo](https://www.gear-genomics.com/indigo). Indigo is primarily for discovering InDels in Sanger traces but it also aligns the trace file to the reference genome so we can use it to compare the alignments and traces of the index patient to her parents. The files are patient.ab1, mother.ab1 and father.ab1. For your convenience, a screenshot of the sanger traces is provided in sanger.png.

```bash
open sanger.png
```

You can of course also spot-check this variant in the raw alignment.

```bash
samtools tview -d t -p chr7:2954850 rd.rmdup.bam chr7.fa
```


***Exercises***

* Why should we not put the primers directly next to the mutation?
* Why did we not select primers more than 1000bp away from the mutation?
* Is the gene of interest on the forward or reverse strand?
* Does the candidate gene make sense for a patient with severe immunodeficiency?
* Does the candidate gene interact with NFKB1?
* Can you spot the mutation in the traces and the alignment against the reference?
* What is the validated genotype by Sanger sequencing of the mother, father and patient for the mutation?
