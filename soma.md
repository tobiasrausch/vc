## Structural Variant Calling

As a first step we will compute structural variants using [Delly](https://github.com/tobiasrausch/delly). [Delly](https://www.ncbi.nlm.nih.gov/pubmed/22962449) calls structural variants jointly on the tumor and normal genome and outputs a [BCF](https://samtools.github.io/hts-specs) file, the binary encoding of [VCF](https://samtools.github.io/hts-specs). You can also provide a text file with regions to exclude from the analysis of structural variants. The default exclude map of Delly includes the telomeric and centromeric regions of all human chromosomes since these regions cannot be accurately analyzed with short-read data.

```shell
delly call -t DEL -q 20 -g chr2.fa -x hg19.ex -o del.bcf tumor.bam control.bam
delly call -t DUP -q 20 -g chr2.fa -x hg19.ex -o dup.bcf tumor.bam control.bam
delly call -t INV -q 20 -g chr2.fa -x hg19.ex -o inv.bcf tumor.bam control.bam
```

VCF was originally designed for short variants and that's why all SV callers heavily use the INFO fields to encode additional information about the SV such as the structural variant end (INFO:END) and the SV type (INFO:SVTYPE). You can look at the header of the BCF file using grep, '-A 1' includes the first structural variant record in the file:

```shell
bcftools view del.bcf | grep "^#" -A 1
```

[Delly](https://github.com/tobiasrausch/delly) uses the VCF:INFO fields for structural variant site information, such as how confident the structural variant prediction is and how accurate the breakpoints are. The genotype fields contain the actual sample genotype, its genotype quality and genotype likelihoods and various count fields for the variant and reference supporting reads and spanning pairs. If you browse through the VCF file you will notice that a subset of the Delly structural variant predictions have been refined using split-reads. These precise variants are flagged in the VCF info field with the tag 'PRECISE', all others are listed as 'IMPRECISE'. Please note that this BCF file contains germline and somatic structural variants but also false positives caused by repeat-induced mis-mappings or incomplete reference sequences. [SVprops](https://github.com/tobiasrausch/svprops) is a simple program that converts Delly's BCF output to a tab-delimited SV site list.

```shell
svprops del.bcf | head
```

Every record of the BCF file is converted to one tab-delimited row with all kinds of summary statistics. [SVprops](https://github.com/tobiasrausch/svprops) also provides a 'column-view' listing summary statistics for all samples present in the BCF file.

```shell
sampleprops del.bcf
```

This initial SV calling cannot differentiate somatic and germline structural variants. For instance, the 2 complex variants we looked at before are still present in the Delly output. A proximal duplication causes 2 paired-end signatures (deletion-type and duplication-type):

```shell
bcftools view del.bcf chr2:18905691-18907969 | awk '$2>=18905691 && $2<=18907969'
bcftools view dup.bcf chr2:18905691-18907969 | awk '$2>=18905691 && $2<=18907969'
```

***Exercises***

* What is the fraction of deletions that has been called precisely (at single nucleotide resolution) by Delly?
* Is the germline inverted duplication present in Delly's output files?
* What type of SVs delineate a proximal inverted duplication?
* [Blat](https://genome.ucsc.edu/cgi-bin/hgBlat) the consensus sequence of the first precise deletion (bcftools view del.bcf | grep "^#" -A 1). Is it a known variant?


## Somatic Filtering

Delly ships with a basic somatic filtering subcommand that uses the matched control and possibly additional control samples from unrelated individuals. The somatic filtering requires a sample file listing tumor and control sample names from the VCF file.

```shell
cat spl.tsv
```

There are many parameters available to tune the somatic structural variant filtering. Below we require a minimum variant allele frequency of 25%, no support in the matched normal and an overall confident structural variant site prediction with the VCF filter field being equal to PASS.

```shell
delly filter -pt DEL -f somatic -o sdel.bcf -a 0.25 -s spl.tsv del.bcf
delly filter -pt DUP -f somatic -o sdup.bcf -a 0.25 -s spl.tsv dup.bcf
delly filter -pt INV -f somatic -o sinv.bcf -a 0.25 -s spl.tsv inv.bcf
```

Using [BCFtools](http://www.htslib.org) we can merge all somatic structural variants together in a single BCF file.

```shell
bcftools concat -a -O b -o somatic.bcf sdel.bcf sdup.bcf sinv.bcf
sampleprops somatic.bcf
```

***Exercises***

* What is the average size of the somatic SVs?
* Is there any SV type that is enriched among the somatic SVs?

## Structural Variant Visualization using IGV

IGV is excellent for inspecting small variants but for these long-range complex re-arrangements its less powerful. You cannot view the entire event in IGV but you can still visualize the breakpoint regions separately, denoted as left breakpoint (L) and right breakpoint (R) below.

```shell
svprops somatic.bcf | tail -n +2 | awk '{print $1"\t"($2-500)"\t"($2+500)"\t"$5"L";}' > somatic.bp.bed
svprops somatic.bcf | tail -n +2 | awk '{print $3"\t"($4-500)"\t"($4+500)"\t"$5"R";}' >> somatic.bp.bed
./igv.sh -g chr2.fa
```

As previously, once IGV has started use 'File' and 'Load from File' to load the tumor and normal bam alignment file. Then import 'somatic.bp.bed' from your working directory using 'Regions' and 'Import Regions'. The somatic structural variants can then be browsed easily using 'Regions' and 'Region Navigator' and you probably want to switch on 'View as pairs' and 'Color alignments by pair orientation'. You may also want to play around with the split screen view, where you can right click a read and select 'View mate region in split screen'.


## Visualizing Complex Rearrangements

To reveal any complex re-arrangement patterns it's worthwhile to create a so-called read-depth ratio plot using [R statistics](https://www.r-project.org) and overlay the somatic structural variants. We therefore first calculate read counts in 10kbp windows using Delly's cov tool.

```shell
cov -g chr2.fa -q 20 -f tn.cov.gz tumor.bam control.bam
R
```

In R we first load and then plot the coverage data.

```R
library(ggplot2)
library(reshape2)
library(scales)
cov = read.table("tn.cov.gz", header=T)
cov$pos = cov$start + (cov$end - cov$start) / 2
cov = cov[cov$control!=0,c("pos", "tumor", "control")]
covNorm = median(cov$control) / median(cov$tumor)
cov$lg2 = log2(covNorm*cov$tumor/cov$control)
cov = melt(cov, id.vars=c("pos"))

p=ggplot(data=cov, aes(x=pos, y=value)) + geom_point(alpha=1/3)
p=p + xlab("chr2") + ylab("log2 | coverage | coverage")
p=p + scale_x_continuous(labels = comma)
p=p + facet_grid(variable ~ ., scales="free")
p
ggsave("cov.pdf", width=9, height=6)
```

When you quit R please save the workspace so we can continue amending this initial coverage plot. Next, we will try to visualize the somatic SV calls on top of the log2 read-depth ratio plot.

```shell
svprops somatic.bcf > svs.tsv
```

In R we then augment the read-depth plot.

```R
library(ggplot2)
library(reshape2)
library(scales)
cov = cov[cov$variable == "lg2",]
p=ggplot(data=cov, aes(x=pos, y=value)) + geom_point(alpha=1/3)
p=p + ylab("Log2 read-depth ratio") + xlab("chr2")
p=p + scale_x_continuous(labels = comma)
p
ggsave("cov1.pdf", width=9, height=6)

sv = read.table("svs.tsv", header=T)
sv$type = factor(paste0(sv$svtype, "-", sv$ct))
sv = sv[,c("chr","start","end","type","id")]
p2 = p + geom_curve(data=sv, aes(x=start, xend=end, col=type), y=2, yend=2, curvature=-0.5)
p2 = p2 + scale_y_continuous(limits=c(-2,3))
p2 = p2 + labs(colour="SV type")
p2
ggsave("cov2.pdf", width=9, height=6)
```

In addition, you might also want to segment the read-depth ratios using the [DNAcopy](https://bioconductor.org/packages/release/bioc/html/DNAcopy.html) R package.

```R
library(ggplot2)
library(reshape2)
library(scales)
library(DNAcopy)
seg=segments.summary(segment(smooth.CNA(CNA(cov$value, rep("chr2", nrow(cov)), cov$pos, data.type="logratio", sampleid="tumor"))))
p3=p2 + geom_segment(data=seg, aes(x=loc.start, y=seg.median, xend=loc.end, yend=seg.median), colour="darkorange")
p3=p3 + theme(legend.position="bottom")
p3
ggsave("cov3.pdf", width=9, height=6)
```

For later retrival save again the workspace when you now quit R.


## Allelic Depth for SNPs

We can also augment the plots with SNP variant allele frequencies. The file 'snps.pos' contains 1000 Genomes SNP locations on chr2 where the control sample is heterozygous.

```shell
head snps.pos
```

Using [SAMtools](http://www.htslib.org) and [BCFtools](http://www.htslib.org) we can annotate the allelic depth for every SNP.

```shell
samtools mpileup -go mp.bcf -f chr2.fa -t AD -l snps.pos tumor.bam
bcftools call -vmAO z -o mp.vcf.gz mp.bcf
zcat mp.vcf.gz | grep "^#" -A 1
```

Let us first filter for bi-allelic SNPs at coverage >=5 and then extract AD (allelic depth). In the last step we calculate the so-called B-allele frequency.

```shell
bcftools filter -i 'DP>=10' mp.vcf.gz | bcftools view -m2 -M2 -v snps - | bcftools query -f '%CHROM\t%POS[\t%AD]\n' > ad.tsv
cat ad.tsv | tr ',' '\t' | awk '{print $2"\t"($4)/($3+$4);}' > baf.tsv
```

Let's look at the allelic counts using R.

```R
library(ggplot2)
library(reshape2)
library(scales)
co = read.table("baf.tsv", header=F)
colnames(co) = c("pos", "baf")
q=ggplot(data=co, aes(x=pos, y=baf))
q=q + geom_jitter(alpha=1/8, size=0.5, width=0, height=0.1)
q=q + xlab("chr2") + ylab("Variant Allele Frequency") 
q=q + scale_x_continuous(labels = comma)
q
```

Of course, this plot is more informative if we can also see the read-depth.

```R
library(grid)
grid.newpage()
pushViewport(viewport(layout=grid.layout(6,1)))
print(q, vp = viewport(layout.pos.row=1:2, layout.pos.col=1))
print(p3, vp = viewport(layout.pos.row=3:6, layout.pos.col=1))
```
