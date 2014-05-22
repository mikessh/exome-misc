==========
EXOME-misc
==========

Miscellaneous tools for exome-seq: annotating SNP, counting molecules on start, etc.

To see the list of executable scripts run 

>$java -jar exome-tools.jar

To execute specific script run 

>$java -cp exome-tools.jar script_name

Description
===========

1. MolCountExome and MolCountSNP
================================

Estimates the number of initial molecules using the number of read (pairs) with unique mapping offsets aka distinct reads.
Those estimates are used to calculate 'molecular coverage' of exome-sequencing data and individual SNPs. 

Usage:

>$groovy MolCountExome input(sam or bam) exome_ref(refseq format) output [read length, default = 100]

>$groovy MolCountSNP input(sam or bam) varinats(VCF format) output [read length, default = 100]

Note:

Exome table should have the following columns (for whole exome just download RefSeq genes UCSC GB table):

```name``` ```chrom``` ```strand``` ```txStart``` ```txEnd``` ```cdsStart``` ```cdsEnd``` ```exonStarts``` ```exonEnds``` ```name2```

2. VcfAnnot and MutectAnnot
===========================

These scripts perform a comprehensive annotation of VCF and MuTect format entries:

- Search for matching dbSNP and COMSIC entries
- Search for parent transcript and segment (exon, intron or UTR)
- Translate mutation-related codons, find missense mutations and frameshift indels

Usage:

>$VcfAnnot.groovy [options] file1.vcf[,file2.vcf,...]|folder/with/vcfs [path to folder with refGene, COSMIC and dbSNP data if not in script dir]

>$MutectAnnot.groovy [options] mutect_output1[,mutect_output2,...] [path to folder with refGene, COSMIC and dbSNP data if not in script dir]

Note:

A bundle (134 Mb) with annotations should be downloaded from https://www.dropbox.com/s/q42rm9cdns9o5fx/vcf_annot_bundle.zip

It should be either placed in script parent directory or specified during execution 