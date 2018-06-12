## Germline Deletions

To get an idea of how structural variants (SVs) look like in short read, next-generation sequencing data one can use the [1000 Genomes SV polymorphism catalogue](https://www.ncbi.nlm.nih.gov/pubmed/26432246) because these are germline SVs and a subset of those will be present in our sample. To increase the chance that our sample of interest is a carrier of one of these germline SVs I have already subsetted the SV list to common SVs with a population allele frequency above 50% on chr2. The below [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) file shows a couple of germline SVs.

```bash
cd /data/sv/
cat svs.bed
```

Using [IGV](http://software.broadinstitute.org/software/igv/) we can then browse these SVs interactively if an X11 display is available. If you only have a terminal we can still use IGV in batch mode and create snapshots of the putative SV loci.

## Interactive IGV 

With an X11 display we can start [IGV](http://software.broadinstitute.org/software/igv/) and browse the variants interactively using the chr2 reference sequence.

```bash
./igv.sh -g chr2.fa
```

Once IGV has started use 'File' and 'Load from File' to load the control.bam alignment file. Then import 'svs.bed' from your working directory using 'Regions' and 'Import Regions'. The structural variants can then be browsed easily using 'Regions' and 'Region Navigator'. Select a structural variant in the Region Navigator and click 'View', which will center the IGV alignment view on the selected structural variant. It's usually best to zoom out once then by clicking on the '-' sign in the toolbar at the top, so you can view all supporting abnormal paired-ends. To highlight the abnormal paired-ends please right click in IGV on the bam file and activate 'View as pairs'. In the same menu, please open 'Color alignments by' and then switch to "pair orientation' for inversions and duplications. For deletions, you want to color the alignments by "insert size". 


## IGV in batch mode

Without an X11 display we can still use IGV to plot snapshots of the putative SVs.

```bash
./runIGVBatch.sh svs.bed
```

Once IGV is finished we can simply browse these snapshots. Some clear germline deletions are for instance.

```bash
open DEL11.png
open DEL13.png
```


## Complex Germline Structural Variants

Most 1000 Genomes SVs are deletions because these are easier to detect in low coverage sequencing data. Deletions cause a drop in read-depth and can be detected by spanning paired-ends of abnormally large insert size. In addition, split-reads have a prefix alignment before the deletion and a suffix alignment after the deletion. Other SV types such as inversions are much harder to detect because these are balanced rearrangements that do not cause a read-depth change. Besides the simple SV types (deletions, duplications, inversions) one can also find more complex SVs in the germline. Two example regions for that are in the file complex.bed.

```bash
./runIGVBatch.sh complex.bed
open complex1.png
open complex2.png
```

As part of the SV consortium we validated the above complex SVs using PacBio. The reads are in a separate FASTA file called pacbio.sv1.fa and pacbio.sv2.fa. We need the subsequence of the reference to create a pairwise dotplot of the PacBio read against the reference. [SAMtools](http://www.htslib.org) is a convenient tool to extract such subsequences of a FASTA file.

```bash
cat complex.bed
samtools faidx chr2.fa chr2:18905691-18907969
samtools faidx chr2.fa chr2:96210505-96212783
```


Please align the above genomic reference subsequences against the respective PacBio read using [Maze](https://gear.embl.de/maze/) available on [gear.embl.de](https://gear.embl.de). 

***Exercises***

* What kind of SV is present in the region chr2:18905691-18907969 ?
* What kind of SV is present in the region chr2:96210505-96212783 ?


## Structural Variant Alignment QC

Paired-end methods can be affected by a skewed insert size distribution, read-depth methods by non-uniform coverage and split-read methods suffer from high sequencing error rates that cause mis-mappings. Prior to any structural variant discovery you should therefore evaluate the quality of the data such as the percentage of mapped reads, singletons, duplicates, properly paired reads and the insert size & coverage distributions. [Picard](http://broadinstitute.github.io/picard/), [SAMtools](http://www.htslib.org), [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [Alfred](https://github.com/tobiasrausch/alfred) compute some of these alignment statistics as shown below for the tumor sample.

```bash
samtools flagstat tumor.bam
alfred qc -r chr2.fa tumor.bam
zcat qc.tsv.gz | grep ^ME | datamash transpose | column -t
Rscript /opt/alfred/R/stats.R qc.tsv.gz
convert qc.tsv.gz.pdf qc.png
# Coverage distribution
open qc-4.png
# Insert size distribution
open qc-5.png
```

There are some general things to watch out for such as mapping percentages below 70%, >20% duplicates or multiple peaks in the insert size distribution. However, these statistics vary largely by protocol and hence, it's usually best to compare multiple different sequencing runs using the same protocol (DNA-seq, RNA-seq, ChIP-seq, paired-end, single-end or mate-pair) against each other, which then highlights the outliers.

***Exercises***

* What is the median coverage of the data set?
* Given the insert size distribution, what would be a suitable cutoff to define deletion supporting paired-ends?

