Setup for marker design
================
Lindsay Clark, HPCBio, Roy J. Carver Biotechnology Center, University of
Illinois, Urbana-Champaign

This document contains setup steps that you should only have to do once.
If you haven’t already, you’ll need to clone this repository to your
computer and open it as an RStudio project as described in the
[README](README.md).

## Packages to install

Certain Bioconductor and CRAN packages must be installed for these
tutorials to work. They can be installed with the following code. It is
best to run this installation code as soon as R is opened, before doing
any other work.

``` r
install.packages(c("tidyverse", "ape", "adegenet"))
BiocManager::install(c("VariantAnnotation", "Rsamtools"))
```

If you are asked to update packages, I recommend answering “a” for “all”
just to make sure everything is up-to-date.

## Data

For convenience, this repository has a `data` folder where you can put
your VCF, reference FASTA, and any other input files. You don’t have to
put them here, but if they are somewhere else you will need to specify
the full path to them.

In this case I have three VCFs for two different projects, but you might
just have one.

``` r
vcf1 <- "data/331_new_data_vcf_IBRC.vcf"
vcf2 <- "data/unique_173_clones.recode.vcf"
vcf3 <- "data/yam336.m2M2vsnps_missing0.9.recode.vcf"
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
library(dplyr)
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
hdr1 <- scanVcfHeader(bg1)
hdr1
```

    ## class: VCFHeader 
    ## samples(331): TDr2946A TDr1489A ... TDr0900280 TDr8700211
    ## meta(2): fileformat Tassel
    ## fixed(1): FILTER
    ## info(3): NS DP AF
    ## geno(5): GT AD DP GQ PL

``` r
hdr2 <- scanVcfHeader(bg2)
hdr2
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

## Importing SNP metadata only from a larger VCF

In some cases, you may have a much larger, unfiltered VCF. Although you
don’t plan to mine markers from the VCF, you still need to know the
location and alleles of those SNPs so they can be annotated in the
flanking sequence that you will use for KASP design. The script
`snp_positions_from_vcf.R` reads the VCF and saves the SNP metadata to
an R object. The script is designed to be run on a bioinformatics server
or cluster from the terminal like so:

``` bash
Rscript snp_positions_from_vcf.R --args yam336.m2M2vsnps_missing0.9.recode.vcf yam336.m2M2vsnps_missing0.9.recode.rds
```

After successfully running it, you can download the `.rds` file to your
`data` directory.

## Quality control

There are a few things you might want to check over in your VCF. We’ll
run through a similar set of checks for the two files, but they will
differ slightly depending on how the files were generated.

### TASSEL

Let’s look at the first 50 sample names. These should match the samples
that were in your study.

``` r
samples(hdr1) %>% head(50)
```

    ##  [1] "TDr2946A" "TDr1489A" "TDr2284A" "TDr1499A" "TDr1509A" "TDr1510A" "TDr3782A" "TDr1858C" "TDr1576A" "TDr1585A" "TDr1585C" "TDr1598A" "TDr1622A" "TDr1628A"
    ## [15] "TDr1631C" "TDr1649A" "TDr1650B" "TDr1653A" "TDr1655A" "TDr1663A" "TDr1686A" "TDr1707A" "TDr1709A" "TDr1711A" "TDr3872A" "TDr1732A" "TDr1735A" "TDr2029A"
    ## [29] "TDr1760A" "TDr1763C" "TDr1804A" "TDr1775A" "TDr1798A" "TDr1805A" "TDr1807A" "TDr1829A" "TDr1850A" "TDr1899A" "TDr1922C" "TDr1935A" "TDr2608A" "TDr2041B"
    ## [43] "TDr2121A" "TDr2155A" "TDr2159A" "TDr2161C" "TDr2167A" "TDr2207A" "TDr2210A" "TDr3311B"

You can also take a look at the variants themselves. Here we’ll read the
SNP metadata without reading the genotypes.

``` r
vcf1 <- readVcf(bg1, genome = genfasta,
               param = ScanVcfParam(geno = NA))
rowRanges(vcf1)
```

    ## GRanges object with 136429 ranges and 5 metadata columns:
    ##                     seqnames    ranges strand | paramRangeID            REF                ALT      QUAL      FILTER
    ##                        <Rle> <IRanges>  <Rle> |     <factor> <DNAStringSet> <DNAStringSetList> <numeric> <character>
    ##      chrom_01_32840    OM_01     32840      * |           NA              G                  A        NA        PASS
    ##      chrom_01_45700    OM_01     45700      * |           NA              A                  T        NA        PASS
    ##      chrom_01_58956    OM_01     58956      * |           NA              T                  C        NA        PASS
    ##      chrom_01_62865    OM_01     62865      * |           NA              A                  T        NA        PASS
    ##      chrom_01_65124    OM_01     65124      * |           NA              G                  T        NA        PASS
    ##                 ...      ...       ...    ... .          ...            ...                ...       ...         ...
    ##   chrom_20_33017477    OM_20  33017477      * |           NA              A                  G        NA        PASS
    ##   chrom_20_33017823    OM_20  33017823      * |           NA              C                  T        NA        PASS
    ##   chrom_20_33017839    OM_20  33017839      * |           NA              T                  C        NA        PASS
    ##   chrom_20_33017917    OM_20  33017917      * |           NA              T                  C        NA        PASS
    ##   chrom_20_33018263    OM_20  33018263      * |           NA              A                  G        NA        PASS
    ##   -------
    ##   seqinfo: 20 sequences from data/TDr96_F1_v2_PseudoChromosome.rev07.fasta genome; no seqlengths

In the `seqnames` column, we see chromosome names that have been
shortened by TASSEL. These don’t match the names in
`seqinfo(refgenome)`, so we’ll have to keep that in mind. There are
functions in the `src` directory to help with name conversion. If your
VCF is different and you need to do some conversion other than “OM” to
“chrom”, you might modify the `fixTasselNames` function.

``` r
source("src/marker_stats.R")

ranges1a <- rowRanges_correctedSeqnames(vcf1)
ranges1a
```

    ## GRanges object with 136429 ranges and 5 metadata columns:
    ##                     seqnames    ranges strand | paramRangeID            REF                ALT      QUAL      FILTER
    ##                        <Rle> <IRanges>  <Rle> |     <factor> <DNAStringSet> <DNAStringSetList> <numeric> <character>
    ##      chrom_01_32840 chrom_01     32840      * |           NA              G                  A        NA        PASS
    ##      chrom_01_45700 chrom_01     45700      * |           NA              A                  T        NA        PASS
    ##      chrom_01_58956 chrom_01     58956      * |           NA              T                  C        NA        PASS
    ##      chrom_01_62865 chrom_01     62865      * |           NA              A                  T        NA        PASS
    ##      chrom_01_65124 chrom_01     65124      * |           NA              G                  T        NA        PASS
    ##                 ...      ...       ...    ... .          ...            ...                ...       ...         ...
    ##   chrom_20_33017477 chrom_20  33017477      * |           NA              A                  G        NA        PASS
    ##   chrom_20_33017823 chrom_20  33017823      * |           NA              C                  T        NA        PASS
    ##   chrom_20_33017839 chrom_20  33017839      * |           NA              T                  C        NA        PASS
    ##   chrom_20_33017917 chrom_20  33017917      * |           NA              T                  C        NA        PASS
    ##   chrom_20_33018263 chrom_20  33018263      * |           NA              A                  G        NA        PASS
    ##   -------
    ##   seqinfo: 20 sequences from an unspecified genome; no seqlengths

Next, we will use the DNA sequence to make sure that the reference
genome matches the VCF. First, we will retrieve the nucleotide at each
SNP position from the reference genome.

``` r
refcheck1 <- scanFa(refgenome, ranges1a)
refcheck1
```

    ## DNAStringSet object of length 136429:
    ##          width seq                                                                                                                               names               
    ##      [1]     1 G                                                                                                                                 chrom_01
    ##      [2]     1 T                                                                                                                                 chrom_01
    ##      [3]     1 T                                                                                                                                 chrom_01
    ##      [4]     1 A                                                                                                                                 chrom_01
    ##      [5]     1 T                                                                                                                                 chrom_01
    ##      ...   ... ...
    ## [136425]     1 A                                                                                                                                 chrom_20
    ## [136426]     1 C                                                                                                                                 chrom_20
    ## [136427]     1 T                                                                                                                                 chrom_20
    ## [136428]     1 T                                                                                                                                 chrom_20
    ## [136429]     1 A                                                                                                                                 chrom_20

Does that match the reference allele in the VCF?

``` r
mean(refcheck1 == rowRanges(vcf1)$REF)
```

    ## [1] 0.8517471

It only matches 85% of the time. However, that is because TASSEL
sometimes lists the common allele as the reference allele, rather than
the allele that matches the reference. So, we’ll also check the
alternative allele.

``` r
mean(rowRanges(vcf1)$REF == refcheck1 | unlist(rowRanges(vcf1)$ALT) == refcheck1)
```

    ## [1] 1

Now we see a 100% match, so everything looks good.

### Samtools

Let’s look at the first 10 sample names. In this case they contain the
full path to the BAM file, so we might want to keep that in mind if we
need to look up genotypes by sample later.

``` r
samples(hdr2) %>% head(10)
```

    ##  [1] "/home/ibrcuser/work/plant2/2019_Yam/1200.PCflye_ac334_ns167_align/PCrev6_all_bams/DRS_002.all.rd.bam"
    ##  [2] "/home/ibrcuser/work/plant2/2019_Yam/1200.PCflye_ac334_ns167_align/PCrev6_all_bams/DRS_003.all.rd.bam"
    ##  [3] "/home/ibrcuser/work/plant2/2019_Yam/1200.PCflye_ac334_ns167_align/PCrev6_all_bams/DRS_004.all.rd.bam"
    ##  [4] "/home/ibrcuser/work/plant2/2019_Yam/1200.PCflye_ac334_ns167_align/PCrev6_all_bams/DRS_006.all.rd.bam"
    ##  [5] "/home/ibrcuser/work/plant2/2019_Yam/1200.PCflye_ac334_ns167_align/PCrev6_all_bams/DRS_007.all.rd.bam"
    ##  [6] "/home/ibrcuser/work/plant2/2019_Yam/1200.PCflye_ac334_ns167_align/PCrev6_all_bams/DRS_009.all.rd.bam"
    ##  [7] "/home/ibrcuser/work/plant2/2019_Yam/1200.PCflye_ac334_ns167_align/PCrev6_all_bams/DRS_012.all.rd.bam"
    ##  [8] "/home/ibrcuser/work/plant2/2019_Yam/1200.PCflye_ac334_ns167_align/PCrev6_all_bams/DRS_016.all.rd.bam"
    ##  [9] "/home/ibrcuser/work/plant2/2019_Yam/1200.PCflye_ac334_ns167_align/PCrev6_all_bams/DRS_018.all.rd.bam"
    ## [10] "/home/ibrcuser/work/plant2/2019_Yam/1200.PCflye_ac334_ns167_align/PCrev6_all_bams/DRS_022.all.rd.bam"

You can also take a look at the variants themselves. Here we’ll read the
SNP metadata without reading the genotypes.

``` r
vcf2 <- readVcf(bg2, genome = genfasta,
               param = ScanVcfParam(geno = NA))
rowRanges(vcf2)
```

    ## GRanges object with 136429 ranges and 5 metadata columns:
    ##                     seqnames    ranges strand | paramRangeID            REF                ALT      QUAL      FILTER
    ##                        <Rle> <IRanges>  <Rle> |     <factor> <DNAStringSet> <DNAStringSetList> <numeric> <character>
    ##      chrom_01_32840 chrom_01     32840      * |           NA              G                  A       999           .
    ##      chrom_01_45700 chrom_01     45700      * |           NA              T                  A       999           .
    ##      chrom_01_58956 chrom_01     58956      * |           NA              T                  C       999           .
    ##      chrom_01_62865 chrom_01     62865      * |           NA              A                  T       999           .
    ##      chrom_01_65124 chrom_01     65124      * |           NA              T                  G       999           .
    ##                 ...      ...       ...    ... .          ...            ...                ...       ...         ...
    ##   chrom_20_33017477 chrom_20  33017477      * |           NA              A                  G       999           .
    ##   chrom_20_33017823 chrom_20  33017823      * |           NA              C                  T       999           .
    ##   chrom_20_33017839 chrom_20  33017839      * |           NA              T                  C       999           .
    ##   chrom_20_33017917 chrom_20  33017917      * |           NA              T                  C       999           .
    ##   chrom_20_33018263 chrom_20  33018263      * |           NA              A                  G       999           .
    ##   -------
    ##   seqinfo: 2253 sequences from data/TDr96_F1_v2_PseudoChromosome.rev07.fasta genome

We can see if the reference allele is what we would expect, given the
sequence from our reference genome.

``` r
refcheck2 <- scanFa(refgenome, rowRanges(vcf2))
mean(rowRanges(vcf2)$REF == refcheck2) # this should return 1
```

    ## [1] 1

The reference allele matches the reference genome 100% of the time.

We can check with our larger VCF as well.

``` r
vcf3 <- readRDS(sub("vcf$", "rds", vcf3))

refcheck3 <- scanFa(refgenome, rowRanges(vcf3))
mean(rowRanges(vcf3)$REF == refcheck3) # this should return 1
```

    ## [1] 1
