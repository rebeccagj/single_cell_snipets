# Adding GFP to a 10X Genomics `cellranger` reference genome

Use the pre-complied 10X mm10 to acquire the genome and gtf; it's a faster download (since it's tarball) than doing them individually. You still have to remake the genome with the files you modify using `cellranger mkref`.

```
wget http://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-mm10-3.0.0.tar.gz
tar -xzvf refdata-cellranger-mm10-3.0.0.tar.gz
```


## Edit fasta

Let's see what the top of the first file, the fasta genome sequence looks like.

```
cd refdata-cellranger-mm10-3.0.0/
cat fasta/genome.fa | head
>1 dna:chromosome chromosome:GRCm38:1:1:195471971:1 REF
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
```
The start of the file looks at the first chromosome. `N` means this file has repeat regions hard masked (read more about this [here](https://uswest.ensembl.org/info/genome/genebuild/assembly_repeats.html) if you care).

Now let's look at how each chromosome is encoded:

```
grep "dna:" genome.fa
>1 dna:chromosome chromosome:GRCm38:1:1:195471971:1 REF
>10 dna:chromosome chromosome:GRCm38:10:1:130694993:1 REF
...
>Y dna:chromosome chromosome:GRCm38:Y:1:91744698:1 REF
>JH584299.1 dna:scaffold scaffold:GRCm38:JH584299.1:1:953012:1 REF
>GL456233.1 dna:scaffold scaffold:GRCm38:GL456233.1:1:336933:1 REF
```

This is the kind of info we need to mimic for our GFP. Lets make a GFP chromosome with the eGFP sequence.

```
touch fasta/eGFP.fa
emacs fasta/eGFP.fa
```

This can be copy and pasted in (replace sequence `ATG...` and sequence length `1:725` with your transcript values):

```
>GFP dna:chromosome chromosome:GRCm38:GFP:1:725:1 REF
ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAA
GTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGC
TGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAG
CAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTA
CAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGG
ACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAAC
GGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACAC
CCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACG
```

(These are spaces in the `>GFP dna:...`, not tabs, which is a clarification you will appreciate later.)

Append to mm10 fasta & double check it.

```
cat fasta/eGFP.fa >> fasta/genome.fa
cat fasta/genome.fa | tail 
```

## Edit GTF

For reference: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references#gtf

We can build our gtf using this table from the 10x link:

|Column|Name|Description|
|---|---|---|
|1|	Chromosome|	Must refer to a chromosome/contig in the genome fasta.|
|2|	Source|	Unused.|
|3|	Feature|	cellranger count only uses rows where this line is exon.|
|4|	Start|	Start position on the reference (1-based inclusive).|
|5|	End	End| position on the reference (1-based inclusive).|
|6|	Score|	Unused.|
|7|	Strand|	Strandedness of this feature on the reference: + or -.|
|8|	Frame|	Unused.|
|9|	Attributes|	A semicolon-delimited list of key-value pairs of the form key "value". The attribute keys transcript_id and gene_id are required; gene_name is optional and may be non-unique, but if present will be preferentially displayed in reports.|

Let's see what this looks like in the pre-built genome by looking at the very first gene

```
sed -n '6p' genes/genes.gtf
1       ensembl_havana  gene    3205901 3671498 .       -       .       gene_id "ENSMUSG00000051951"; gene_version "5"; gene_name "Xkr4"; gene_source "ensembl_havana"; gene_biotype "protein_coding"

```

We can mimic this for GFP as well.

```
touch genes/eGFP.gtf
emacs genes/eGFP.gtf
```

And paste in the following:

```
GFP	annot	exon	1	725	.	+	.	gene_id "GFP"; transcript_id "eGFP"; gene_name "GFP"
```

Notes:

* Note that `GFP` chromosome name exactly matches the `eGFP.fa` name. 
* I made up the `Source` value. 
* Feature column **must** be `exon` for cellranger STAR parameters. 
* I think `transcript_id` and `gene_id` should be unique. 
* **IMPORTANT**: Column 1:9 should be separated by tabs, not spaces!! Edit this in your emacs/vim if necessary.


Append to mm10 GTF & double check it.

```
cat genes/eGFP.gtf >> genes/genes.gtf
cat genes/genes.gtf | tail
```

## Make New Reference

Install `cellranger` if not already installed to create the custom mm10 + GFP

```
mv fasta/genome.fa ~/mm10_eGFP_fasta.fa
mv genes/genes.gtf ~/mm10_eGFP_annot.gtf
cellranger mkref --genome=mm10gfp --fasta=~/mm10_eGFP_fasta.fa --genes=~/mm10_eGFP_annot.gtf
```