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

***Exercises***

* What is the length of chr7?
* What is the GC-content of chr7? (hint: bedtools nuc)
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
alfred qc -r chr7.fa -o stats rd.rmdup.bam
cat stats.metrics.tsv | datamash transpose | column -t
```

All distribution files are simple tab-delimited text files that can be easily visualized in [R](https://www.r-project.org/).


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
* Would it make sense to sequence this library deeper to achieve 30x coverage?


