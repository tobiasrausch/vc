## Germline Structural Variants

To get an idea of how structural variants look like in short read, next-generation sequencing data let's download the [1000 Genomes SV polymorphism catalogue](https://www.ncbi.nlm.nih.gov/pubmed/26432246).

```shell
wget 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz'
wget 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz.tbi'
```

To increase the chance that our sample of interest is a carrier of one of these germline SVs we will subset to common SVs with a population allele frequency above 80% on chr2.

```shell
bcftools view ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz 2 | bcftools filter -O z -o common.vcf.gz -i 'EUR_AF>=0.8' -
/g/solexa/home/rausch/scripts/cpp/svprops/src/svprops common.vcf.gz  | tail -n +2 | cut -f 1,2,4,5 > svs.bed
```

We then start IGV using the chr2 reference sequence.

```shell
./igv.sh -g chr2.fa
```

Once IGV has started use 'File' and 'Load from File' to load the control.bam alignment file. Then import 'svs.bed' from your working directory using 'Regions' and 'Import Regions'. The structural variants can then be browsed easily using 'Regions' and 'Region Navigator'. Select a structural variant in the Region Navigator and click 'View', which will center the IGV alignment view on the selected structural variant. It's usually best to zoom out once then by clicking on the '-' sign in the toolbar at the top, so you can view all supporting abnormal paired-ends. To highlight the abnormal paired-ends please right click in IGV on the bam file and activate 'View as pairs'. In the same menu, please open 'Color alignments by' and then switch to "pair orientation' for inversions and duplications. For deletions, you want to color the alignments by "insert size". 


As you will notice most 1000 Genomes SVs are deletions because these are easier to detect in low coverage sequencing data. Deletions cause a drop in read-depth and can be detected by spanning paired-ends of abnormally large insert size. In addition, split-reads have a prefix alignment before the deletion and a suffix alignment after the deletion. Other SV types such as inversions are much harder to detect because these are balanced rearrangements that do not cause a read-depth change. Besides the simple SV types (deletions, duplications, inversions) one can also find more complex SVs in the germline. Please have a look at the below regions in IGV.

```shell
samtools faidx chr2.fa chr2:18905691-18907969
samtools faidx chr2.fa chr2:72439663-72441941
samtools faidx chr2.fa chr2:96210505-96212783
```

As part of the SV consortium we validated the above complex SVs using PacBio. The reads are in a separate FASTA file called pacbio.fa. Please align the above genomic reference subsequences against these PacBio reads using [gear.embl.de](gear.embl.de). 

***Exercises***

* What kind of SV is happening in the region chr2:18905691-18907969 ?
* What kind of SV is happening in the region chr2:72439663-72441941 ?
* What kind of SV is happening in the region chr2:96210505-96212783 ?


## Structural Variant Alignment QC

Paired-end methods can be affected by a skewed insert size distribution, read-depth methods by non-uniform coverage and split-read methods suffer from high sequencing error rates that cause mis-mappings. Prior to any structural variant discovery you should therefore evaluate the quality of the data such as the percentage of mapped reads, singletons, duplicates, properly paired reads and the insert size & coverage distributions. [Picard](http://broadinstitute.github.io/picard/), [SAMtools](http://www.htslib.org), [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [Alfred](https://github.com/tobiasrausch/alfred) compute some of these alignment statistics as shown below for the tumor sample.

```shell
samtools flagstat tumor.bam
alfred qc -r chr2.fa -o stats tumor.bam
cat stats.metrics.tsv | datamash transpose | column -t
```

There are some general things to watch out for such as mapping percentages below 70%, >20% duplicates or multiple peaks in the insert size distribution. However, these statistics vary largely by protocol and hence, it's usually best to compare multiple different sequencing runs using the same protocol (DNA-seq, RNA-seq, ChIP-seq, paired-end, single-end or mate-pair) against each other, which then highlights the outliers.

All distribution files from [Alfred](https://github.com/tobiasrausch/alfred) are simple tab-delimited text files that can be easily visualized in [R](https://www.r-project.org/).


```R
library(ggplot2)
cov=read.table("stats.coverage.tsv", header=T)
p=ggplot(data=cov, aes(x=Coverage, y=Count))
p=p + geom_line()
p=p + coord_cartesian(xlim=c(0,50))
p
isize=read.table("stats.isize.tsv", header=T)
q=ggplot(data=isize, aes(x=InsertSize, y=Count))
q=q + geom_line(aes(group=Layout, color=Layout))
q=q + coord_cartesian(xlim=c(0,600))
q
quit()
```

***Exercises***

* What is the median coverage of the data set?
* What is the meaning of the different library layouts (F+, F-, R+, R-)?
* What is the duplicate fraction in the library?

