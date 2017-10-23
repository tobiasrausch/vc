## Variant Annotation

Variant annotation and classification is a challenging process. You can

* use transcript annotations from Ensembl, UCSC or RefSeq
* there is a long list of mutation prediction tools such as PolyPhen, MutationTaster or Sift
* you can annotate variants with allele frequency information from variation archives such as dbSNP, ExAC or gnomAD
* you can check the expression of genes in your studied tissue using GTEx.

In the recent years a number of convenient pipelines have been developed that ease the annotation of variants with some of the above information. In this practical we will use [VEP](http://www.ensembl.org/info/docs/tools/vep/index.html) because it can be run directly online. We first dump all SNPs in a VEP compliant format which you can
then copy and paste into the VEP application. Make sure you use the hg19/GRCh37 version available [here](http://grch37.ensembl.org/Homo_sapiens/Tools/VEP).

```shell
bcftools query -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\n" exon.vcf.gz
```

A working variant annotation and classification pipeline can easily reduce an initial call set of several thousands of exonic variants to a handful mutation candidates. In a rare disease setting additional power can be gained by taking advantage of the suspected inheritance model (autosomal recessive, autosomal dominant, etc.). 

```shell
bcftools query -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%GT]\n" exon.vcf.gz
```

***Exercises***

* What would be a useful additional filter in our case given that the index patient has consanguineous parents?
* How could we use variation archives to further filter the list of exonic variants?
* Which annotation features could be used to rank the mutation list for a clinician?


## Variant Validation

Once a putative causative variant has been identified these are usually validated in the index patient using PCR and Sanger sequencing. If a specific inheritance model is suspected the parents are also tested. For the likely causative variant we will design primers using [Primer3Plus](https://www.ncbi.nlm.nih.gov/pubmed/17485472) available on [gear.embl.de](https://gear.embl.de). We will first select a left primer in the preceeding sequence of the mutation and then a right primer in the suceeding sequence. Don't forget to use hg19/GRChr37 as the reference genome.

```shell
samtools faidx chr7.fa chr7:2954300-2954850
samtools faidx chr7.fa chr7:2954900-2955450
```

The primers are locally unique and have the appropriate Tm but they are not necessarily unique in the entire genome. We can use [Silica](https://gear.embl.de/silica), a tool for [In-silico PCR](https://en.wikipedia.org/wiki/In_silico_PCR), to check genome-wide uniqueness. Try different combinations of left and right primers and possibly change parameters in Primer3Plus to generate further candidates until you found a good pair of primers to validate the mutation.

We do not have the time to run the PCR experiment and sequence the breakpoint mutation but the Sanger validation files of the original study are in the data folder:

You can analyze these trace files using [Indigo](https://gear.embl.de/indigo). Indigo is primarily for discovering InDels in Sanger traces but it also aligns the trace file to the reference genome so we can use it to compare the alignments and traces of the index patient to her parents. The files are patient.ab1, mother.ab1 and father.ab1.

***Exercises***

* Why should we not put the primers directly next to the mutation?
* Why did we not select primers more than 1000bp away from the mutation?
* Is the gene of interest on the forward or reverse strand?
* Can you spot the mutation in the traces and the alignment against the reference?
* What is the validated genotype by Sanger sequencing of the mother, father and patient for the given mutation?
