## Structural Variant Calling

As a first step we will compute structural variants using [Delly](https://github.com/dellytools/delly). [Delly](https://www.ncbi.nlm.nih.gov/pubmed/22962449) calls structural variants jointly on the tumor and normal genome and outputs a [BCF](https://samtools.github.io/hts-specs) file, the binary encoding of [VCF](https://samtools.github.io/hts-specs). You can also provide a text file with regions to exclude from the analysis of structural variants. The default exclude map of Delly includes the telomeric and centromeric regions of all human chromosomes since these regions cannot be accurately analyzed with short-read data.

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

## Interactive IGV

With an X11 display we can start [IGV](http://software.broadinstitute.org/software/igv/) and browse the variants interactively using the chr2 reference sequence.

```bash
./igv.sh -g chr2.fa
```

As previously, once IGV has started use 'File' and 'Load from File' to load the tumor and normal bam alignment file. Then import 'somatic.bp.bed' from your working directory using 'Regions' and 'Import Regions'. The somatic structural variants can then be browsed easily using 'Regions' and 'Region Navigator' and you probably want to switch on 'View as pairs' and 'Color alignments by pair orientation'. You may also want to play around with the split screen view, where you can right click a read and select 'View mate region in split screen'.


## IGV in batch mode

Without an X11 display we can still use IGV to plot snapshots of the putative SVs.

```bash
./runIGVBatch.sh somatic.bp.bed
# For instance, an inversion (left breakpoint)
open INV00000238L.png
# and the right breakpoint
open INV00000238R.png
```


## Visualizing Complex Rearrangements

To reveal any complex re-arrangement patterns it's worthwhile to create a so-called read-depth ratio plot using [R statistics](https://www.r-project.org) and overlay the somatic structural variants. We therefore first calculate fragement counts in non-overlapping 10kbp windows using [Alfred](https://github.com/tobiasrausch/alfred).

```bash
alfred count_dna -o tumor.cov.gz tumor.bam
zcat tumor.cov.gz | head
alfred count_dna -o control.cov.gz control.bam
zcat control.cov.gz | head
```

[Alfred](https://github.com/tobiasrausch/alfred) has a simple Rscript to plot these raw fragment counts across the genome. In our case, 

```bash
Rscript /opt/alfred/R/rd.R tumor.cov.gz
convert tumor.wholegenome.pdf tumor.wholegenome.png
open tumor.wholegenome.png
Rscript /opt/alfred/R/rd.R control.cov.gz
convert control.wholegenome.pdf control.wholegenome.png
open control.wholegenome.png
```

For a matched tumor-control sample it is common practice to normalize the tumor against the control (log2 read-depth ratio plot). This method is implement in the Rscript rdbaf.R.

```bash
paste <(zcat tumor.cov.gz) <(zcat control.cov.gz | cut -f 5) | gzip -c > cov.gz
Rscript rdbaf.R cov.gz
open cov.png
```

Next, we will try to visualize the somatic SV calls on top of the log2 read-depth ratio plot.

```bash
svprops somatic.bcf > svs.tsv
Rscript rdbaf.R cov.gz svs.tsv
open cov.png
```

## Allelic Depth for SNPs

We can also augment the plots with SNP variant allele frequencies. The file 'snps.pos' contains 1000 Genomes SNP locations on chr2 where the control sample is heterozygous.

```bash
head snps.pos
```

Using [SAMtools](http://www.htslib.org) and [BCFtools](http://www.htslib.org) we can annotate the allelic depth for every SNP.

```bash
samtools mpileup -go mp.bcf -f chr2.fa -t AD -l snps.pos tumor.bam
bcftools call -vmAO z -o mp.vcf.gz mp.bcf
zcat mp.vcf.gz | grep "^#" -A 1
```

Let us first filter for bi-allelic SNPs at coverage >=5 and then extract AD (allelic depth). In the last step we calculate the so-called B-allele frequency.

```bash
bcftools filter -i 'DP>=10' mp.vcf.gz | bcftools view -m2 -M2 -v snps - | bcftools query -f '%CHROM\t%POS[\t%AD]\n' > ad.tsv
cat ad.tsv | tr ',' '\t' | awk '{print $2"\t"($4)/($3+$4);}' > baf.tsv
```

Next, we integrate the B-allele frequency plot with the read-depth plot and the structural variants.

```bash
Rscript rdbaf.R cov.gz svs.tsv baf.tsv
open cov.png
```
