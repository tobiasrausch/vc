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

```shell
bcftools view snv.vcf.gz | grep "^#" -A 1
```

Using [BCFtools](https://samtools.github.io/bcftools/bcftools.html) we can generate some useful summary statistics such as the [transition/transversion ratio](https://en.wikipedia.org/wiki/Transversion).

```shell
bcftools stats snv.vcf.gz | grep "TSTV"
```

***Exercises***

* How many SNPs have been called (hint: bcftools stats, SN tag)?
* How many InDels have been called (hint: bcftools stats, SN tag)?
* How many C>T mutations have been called (hint: bcftools stats, ST tag)?


## Filtering Variants

In most applications researchers use external ground truth data to calibrate a variant calling pipeline. In our case we do not know the ground truth so we will illustrate some filtering options based on summary statistics such as the transition/transversion ratio. In most species, transitions are far more likely than transversions and for humans we would expect a transition/transversion ratio of approximately 2.

```shell
bcftools stats snv.vcf.gz | grep "TSTV"
bcftools filter -i '%QUAL>20' snv.vcf.gz  | bcftools stats | grep "TSTV"
bcftools filter -e '%QUAL<=20 || %QUAL/AO<=2 || SAF<=2 || SAR<=2' snv.vcf.gz  | bcftools stats | grep "TSTV"
```

Another useful bulk metric is the length of indels in exons because most InDel polymorphisms should be in-frame. In order to check this we first need to get exon coordinates.

```R
library(GenomicFeatures)
db=makeTxDbFromUCSC(genome="hg19", tablename="ccdsGene")
ex=keepStandardChromosomes(reduce(exons(db), ignore.strand=T))
df=data.frame(chr=seqnames(ex), start=start(ex), end=end(ex))
write.table(df, "exons.bed", quote=F, row.names=F, sep="\t")
```

The exon file needs to be compressed and indexed to calculate in-frame and out-frame InDel statistics.

```shell
bgzip exons.bed
tabix -S 1 -s1 -b2 -e3 exons.bed.gz
bcftools stats -E exons.bed.gz snv.vcf.gz | grep "FS"
```

As you can see, we have only one in-frame exonic InDel because our data is downsampled. This metric is most useful for a large population VCF file with hundreds of samples. Similarly, [heterozygosity](https://en.wikipedia.org/wiki/Zygosity) is a useful metric for such population cohorts, which tends to be higher in Africans compared to Europeans or East Asians. In our case we move on with our simple threshold based filtering and subset the VCF to exonic variants.

```shell
bcftools filter -O z -o exon.vcf.gz -R <(zcat exons.bed.gz | tail -n +2) -e '%QUAL<=20 || %QUAL/AO<=2 || SAF<=2 || SAR<=2' snv.vcf.gz
bcftools stats exon.vcf.gz | egrep "SN|TSTV"
```

***Exercises***

* Plot the InDel length distribution of all called InDels (hint: bcftools stats, IDD tag).

