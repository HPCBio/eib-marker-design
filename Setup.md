Setup for marker design
================

This document contains setup steps that you should only have to do once.

## Packages to install

Certain Bioconductor packages must be installed for these tutorials to
work. They can be installed with the following code. It is best to run
this installation code as soon as R is opened, before doing any other
work.

``` r
BiocManager::install(c("VariantAnnotation", "Rsamtools"))
```

If you are asked to update packages, I recommend answering “a” for “all”
just to make sure everything is up-to-date.

## Data

For convenience, this repository has a `data` folder where you can put
your VCF, reference FASTA, and any other input files. You don’t have to
put them here, but if they are somewhere else you will need to specify
the full path to them.

In this case I have two VCFs for two different projects, but you might
just have one.

``` r
vcf1 <- "data/331_new_data_vcf_IBRC.vcf"
vcf2 <- "data/unique_173_clones.recode.vcf"
genfasta <- "data/TDr96_F1_v2_PseudoChromosome.rev07.fasta"
```

## Building indices

An index is a type of file that enables bioinformatics software to
quickly find any genomic region of interest within a file. For example,
if you needed the sequence from a region of chromosome 5, an index of
your reference FASTA file will help the software to jump directly to
that region, rather than having to read all of chromosomes 1 through 4
first.

``` r
library(VariantAnnotation)
library(Rsamtools)
```

First we’ll index the reference genome.

``` r
indexFa(genfasta)
```

You should now have a file called
data/TDr96\_F1\_v2\_PseudoChromosome.rev07.fasta.fai. If this worked,
you’ll be able to load the reference genome into R:

``` r
refgenome <- FaFile(genfasta)
seqinfo(refgenome)
```

    ## Seqinfo object with 20 sequences from an unspecified genome:
    ##   seqnames seqlengths isCircular genome
    ##   chrom_01   30583384       <NA>   <NA>
    ##   chrom_02   23804961       <NA>   <NA>
    ##   chrom_03   19070176       <NA>   <NA>
    ##   chrom_04   22296860       <NA>   <NA>
    ##   chrom_05   32752956       <NA>   <NA>
    ##   ...             ...        ...    ...
    ##   chrom_16   24109376       <NA>   <NA>
    ##   chrom_17   22358743       <NA>   <NA>
    ##   chrom_18   24280954       <NA>   <NA>
    ##   chrom_19   31611005       <NA>   <NA>
    ##   chrom_20   33023525       <NA>   <NA>

Next, we’ll index the VCFs. We will also compress them first into a
special format for Bioconductor/Samtools.

``` r
bg1 <- bgzip(vcf1)
bg2 <- bgzip(vcf2)
indexTabix(bg1, format = "vcf")
indexTabix(bg2, format = "vcf")
```

The code below gets those file names back without redoing the
compression.

``` r
bg1 <- paste0(vcf1, ".bgz")
bg2 <- paste0(vcf2, ".bgz")
```

You should see “.bgz” and “.bgz.tbi” files in your data folder now. If
your compressing and indexing worked, you should be able to preview the
VCF now.

``` r
scanVcfHeader(bg1)
```

    ## class: VCFHeader 
    ## samples(331): TDr2946A TDr1489A ... TDr0900280 TDr8700211
    ## meta(2): fileformat Tassel
    ## fixed(1): FILTER
    ## info(3): NS DP AF
    ## geno(5): GT AD DP GQ PL

``` r
scanVcfHeader(bg2)
```

    ## class: VCFHeader 
    ## samples(173): /home/ibrcuser/work/plant2/2019_Yam/1200.PCflye_ac334_ns167_align/PCrev6_all_bams/DRS_002.all.rd.bam
    ##   /home/ibrcuser/work/plant2/2019_Yam/1200.PCflye_ac334_ns167_align/PCrev6_all_bams/DRS_003.all.rd.bam ...
    ##   /home/ibrcuser/work/plant2/2019_Yam/1200.PCflye_ac334_ns167_align/PCrev6_all_bams/TDr_199.all.rd.bam
    ##   /home/ibrcuser/work/plant2/2019_Yam/1200.PCflye_ac334_ns167_align/PCrev6_all_bams/TDr_200.all.rd.bam
    ## meta(5): fileformat reference samtoolsCommand samtoolsVersion contig
    ## fixed(2): FILTER ALT
    ## info(17): INDEL IDV ... DP4 MQ
    ## geno(7): GT PL ... GQ GP
