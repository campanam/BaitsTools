# BaitsTools Tutorial

Michael G. Campana, 2017  
Smithsonian Conservation Biology Institute  
Contact: campanam@si.edu  

This series of tutorials covers basic operation of BaitsTools. Each tutorial covers a subcommand and any relevant methods specific to that subcommand. The first two tutorials also demonstrate the basics of bait quality control and filtration. Each tutorial assumes that the previous tutorial has been completed.  

## Table of Contents  
1. [tilebaits and Basic Sequence Filtration](#1-tilebaits-and-basic-sequence-filtration)  
2. [aln2baits and Sequence Filtration by GC Content and Melting Temperature](#2-aln2baits-and-sequence-filtration-by-gc-content-and-melting-temperature)  
3. [bed2baits](#3-bed2baits)  
4. [annot2baits](#4-annot2baits)  
5. [checkbaits](#5-checkbaits)  
6. [vcf2baits](#6-vcf2baits)  
7. [stacks2baits](#7-stacks2baits)  
8. [References](#8-references)

## 1. tilebaits and Basic Sequence Filtration 
tilebaits generates baits from a list of DNA sequences in FASTA or FASTQ format. In this tutorial, we will generate baits from African wild dog (*Lycaon pictus*) mitogenomes and a list of nuclear genes putatively affecting their pelage [1]. We will then apply basic filtration parameters to remove candidate baits that are too short or contain gaps or Ns. Finally, we will examine the bait filtration parameters output file.  

1. Locate the "lycaon_mito.fa" and "pelage_genes.fa" files in the BaitsTools/example_data/ and move them to the location where you will execute baitstools.rb.  

2. Open both files in a text editor. Note that the two mitogenomes in the lycaon_mito.fa file have "#circ" appended to the end of the header lines. This indicates to BaitsTools to treat these sequences as circular. Also note that the RSPO2 sequences have long gaps at the beginning of the sequences. These gaps will be removed during the bait filtration step.  

3. Concatenate the two files into one fasta file "targets.fa". E.g, using cat:  

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`cat *.fa > targets.fa`.  

4. We will generate 60 bp baits with 30 bp tiling to produce a 2× coverage of the input sequences. The input file is specified using `-i`, the bait length is entered using `-L` and the tiling offset uses `-O`. The `-B` option will output a BED format file giving the locations of the generated baits relative to the input sequences. Enter the following command:  

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`ruby baitstools.rb tilebaits -i targets.fa -L 60 -O 30 -B`  

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;BaitsTools will generate a fasta file "out-baits.fa" that contains the unfiltered candidate bait sequences and a  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;BED file "out-baits.bed" listing the bait sequence coordinates.  

5. Use the `-c` argument to remove candidate baits shorter than the requested bait length. Gaps and Ns can be removed using `-G` and `-N` respectively. Use the `-w` command to output a tsv format table of parameters. Enter the following command:  

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`ruby baitstools.rb tilebaits -i targets.fa -L 60 -O 30 -B -w -c -G -N`  

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;In addition to the output from step 4, baitstools will generate a fasta file out-filtered-baits.fa" and a  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;"out-filtered-baits.bed" containing the filtered bait sequences and coordinates respectively. It will also  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;output a tsv-format file "out-filtered-params.txt" containing the results of bait filtration.  

6. Open the out-filtered-params.txt in a spreadsheet program or text editor. The columns from left to right are:  
* Chromosome (Haplotype): Coordinates: The reference sequence, the haplotype identification (for alternate alleles or alignment data), and the bait coordinates.  
* BaitLength: The length of the generated bait  
* %GC: Bait GC content in percent  
* Tm: Bait melting temperature  
* MeanQuality: Mean Phred-like base quality of the generated bait  
* MinQuality: Minimum Phred-like base quality of the generated bait  
* Kept: Whether the bait was retained in the final filtered set.  

## 2. aln2baits and Sequence Filtration by GC Content and Melting Temperature  
aln2baits generates baits from an alignment file in FASTA or FASTQ format. In this tutorial, we will produce a weighted bait set from an alignment of wild dog mitogenomes [1] and reference data from GenBank (accessions: KT448283.1, NC_008093.1, NC_002008.4 [2-4]) such that more variable regions will have increased bait coverage. We will then filter these baits by GC content and melting temperature.

1. Locate the "canid_mito_aln.fa" file in the BaitsTools/example_data/ and move it to the location where you will execute baitstools.rb.  

2. Open the "canid_mito_aln.fa" file in a text editor. Note that the sequences have NCBI-style headers with a sequence name followed by information about the sample separated by spaces. We will therefore need to specify `-D` to inform BaitsTools that this information exists and is not part of the sequence name.  

3. The specification of bait length (`-L`), tiling offset (`-O`) and input file (`-i`) is the same as tilebaits (see Tutorial 1). This time we will use the defaults (120 bp baits and 20 bp offset), so these two parameters can be omitted. However, we will have to specify how aln2baits will generate baits. It can either produce baits for all unique haplotypes observed within the bait window (`haplotype`) or produce baits corresponding to all permutations of single nucleotide variants observed within that window (`variant`). This is specified using the `-H` argument. *WARNING: Variant permutation is computationally intensive and can be very slow.*  
We will produce baits corresponding to haplotypes. Enter the following command:  

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`ruby baitstools.rb aln2baits -i canid_mito_aln.fa -H haplotype -D -B`  

4. Open the "out-baits.bed" file in a text editor. Notice that all coordinates are given in reference to the source sequence in the alignment since each haplotype derives from a single sequence. Using the `variant` option would list all coordinates relative to the position in the alignment rather than a specific sequence.    

5. We will now filter the baits by target GC content and melting temperature. To use a default value, simply use the argument flag without a specified value. We will use the default minimum (`-n`) and maximum (`-x`) GC contents of 30 and 50% respectively. We will change the minimum (`-q`) and maximum (`-z`) melting temperatures to 80C and 140C and the hybridization type (`-T`) to DNA-DNA. Since we will change these from the defaults, the melting temperature values need to be specified after the corresponding argument flags. Enter the following command:  

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`ruby baitstools.rb aln2baits -i canid_mito_aln.fa -H haplotype -D -B -w -n -x -q 80 -z 140 -T DNA-DNA`  

6. Inspect the parameters file to see which baits were retained after filtration.  

## 3. bed2baits  
bed2baits generates baits from a BED file and the corresponding regions in a FASTA or FASTQ file. In this tutorial, we will use bed2baits to generate baits from a subsection of DNA reference sequences.  

1. Locate the "canid_mito_aln.bed" file in the BaitsTools/example_data/ and move it to the location where you will execute baitstools.rb. This BED file corresponds to the "canid_mito_aln.fa" file from the previous tutorial.  

2. For bed2baits, the input BED file is specified by `-i` while the reference sequence(s) are specified using `-r`. All other parameters are identical to tilebaits. We will generate baits under the default settings for the BED file. Enter the following command:  

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`ruby baitstools.rb bed2baits -i canid_mito_aln.bed -r canid_mito_aln.fa -D -B`  

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;In addition to the standard output, the sequences from the extracted regions will be given in the FASTA file  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;"out-regions.fa".   

## 4. annot2baits
bed2baits generates baits from a GFF/GFT annotation file and the corresponding regions in a FASTA or FASTQ file. In this tutorial, we will use annot2baits to generate baits from extracted genome features from a DNA reference sequences.  

1. Locate the "Ananku.fa" and "Ananku.gff" files in the BaitsTools/example_data/ and move them to the location where you will execute baitstools.rb.  

2. annot2baits works almost identically to bed2baits (see previous tutorial) with the GFF/GTF file specified using `-i`. However, the one significant difference is that the user must specify the features to extract using the `-U` argument. Multiple features can be listed using the a comma-separated list. Features must exactly match the term used in the annotation file. Here we will generate baits from all genes and tRNAs from the *Lycaon pictus* mitogenome sequence "Ananku.fa". Enter the following command:  

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`ruby baitstools.rb annot2baits -i Ananku.gff -r Ananku.fa -U gene,tRNA -B`  

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;As with bed2baits, annot2baits will give the standard output and a FASTA file of the extracted regions (in this  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;case "out-regions.fa").  

## 5. checkbaits  
checkbaits performs quality control and filtration on previously generated bait sequences in FASTA or FASTQ format. Here we will perform quality control on the baits generated in the annot2baits tutorial.  

1. checkbaits requires only a input file in FASTA or FASTQ format (`-i` argument) and a desired probe length (`-L` argument). Note that checkbaits does not output BED files listing bait coordinates as this information is not provided to checkbaits. Bait filtration parameters are the same as with other subcommands. We will filter the baits listed in "out-baits.fa" from the previous tutorial by a bait length of 120 bp (requiring the complete length) and a GC content of 40% to 80%. Enter the following command:  

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`ruby baitstools.rb checkbaits -i out-baits.fa -L 120 -c -n 40 -x 80 -w`  

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Standard output will then be given (except for BED coordinates).  

## 6. vcf2baits  
vcf2baits selects variants from a vcf file and then generates baits from a reference sequence. In this tutorial, we will generate baits from *Lycaon pictus* variants mapped to the CanFam3.1 X chromosome sequence (GenBank accession: NC_006621.3 [5]). 

1. Download the X chromosome sequence (in FASTA format) from GenBank (available [here](https://www.ncbi.nlm.nih.gov/nuccore/NC_006621.3)). Downloads from GenBank are typically named "sequence.fasta". Locate the "WDF20_X.raw.vcf.gz" file in the BaitsTools/example_data/. Move both files to the location where you will execute baitstools.rb.  

2. Decompress the "WDF20_X.raw.vcf.gz" file. E.g., using gunzip:  

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`gunzip WDF20_X.raw.vcf.gz`.  

3. Remove extraneous line breaks from the X chromosome reference sequence. These line breaks will slow the program down tremedously. A simple way is using awk:  

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' \`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`sequence.fasta > canfamX.fa`  

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;This tutorial will assume that the line-break-removed fasta file is named "canfamX.fa".  

4. As with the previous subcommands, specify the input VCF using `-i` and the reference sequence with `-r`. In addition to the common output and quality-control options, vcf2baits can filter VCF variants by their QUAL score (specified using `-V`) to help prevent the inclusion of sequencing errors in the bait set. The total number of desired variants is controlled with `-t` across all sequences in the reference file. The maximum number of variants per contig or scafoold is controlled with `-m`. The minimum distance distance between variants is controlled with `-d`. First, we will generate baits for 40 variants with QUAL score >= 30 and a minimum distance of 20,000 bp apart. Note that both `-t` and `-m` are set to 40 here since there is only one reference sequence in the file. Enter the following command:  

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`ruby baitstools.rb vcf2baits -i WDF20_X.raw.vcf -r canfamX.fa -V 30 -t 40 -m 40 -d 20000 -D`

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;This will produce a set of 120 bp baits with the selected variants at the 61st base position in each bait.  

5. Bait length is controlled with `-L`. To change the position of the variant within the bait, use `-b` and `-a` to specify the number of bases before and after the variant within the bait. For example, to change the previous settings to 80 bp baits with the variants at the 21st base postion, enter the following command:  

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`ruby baitstools.rb vcf2baits -i WDF20_X.raw.vcf -r canfamX.fa -V 30 -t 40 -m 40 -d 20000 -b 20 -a 59 -D`  

6. Baits can be tiled across the selected variants by using the `-u` argument. If tiling, use the `-O` argument to specify the base pair offset between tiled baits and the `-k` argument to specify the number of baits per variant. *WARNING: If tiling, do not use the `-b` and `-a` arguments as these are incompatible.* For example, to generate 80 bp baits from with 3× variant coverage and a 15 bp offset between baits, enter the following command:  

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`ruby baitstools.rb vcf2baits -i WDF20_X.raw.vcf -r canfamX.fa -V 30 -t 40 -m 40 -d 20000 -L 80 -u \` 
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`-O 15 -k 3 -D`  

7. Use `-p` to select variants without generating corresponding baits. No reference sequence is needed if `-p` is specified. Enter the following command:  

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`ruby baitstools.rb vcf2baits -i WDF20_X.raw.vcf -V 30 -t 40 -m 40 -d 20000`

8. Use `-e` to generate baits from every variant within a VCF file. For instance, to generate 120 bp baits from the variants selected in the previous step, enter the following command:  

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`ruby baitstools.rb vcf2baits -i WDF20_X.raw.vcf-selected.vcf -r canfamX.fa -e -D`  

9. Use `-j` to scale the maximum number of selected variants per contig by individual contig length. This argument overrides `-m`. Enter the following command:  

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`ruby baitstools.rb vcf2baits -i WDF20_X.raw.vcf -r canfamX.fa -V 30 -j -D`

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;vcf2baits will select up to 30,000 variants that are spaced minimally 10,000 bp apart within the X chromosome.  

10. Finally, use `-R` to apply alternate alleles to baits to generate a balanced bait set. Enter the following command:  

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`ruby baitstools.rb vcf2baits -i WDF20_X.raw.vcf -r canfamX.fa -V 30 -j -R -D`

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;vcf2baits will select up to 30,000 variants that are spaced minimally 10,000 bp apart within the X chromosome  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;and apply alternate alleles to the bait sequences.  

## 7. stacks2baits  
vcf2baits selects variants from a Stacks [6-7] population summary statistics file and then generates baits from a reference sequence. Here we will use stacks2baits to select and sort variants from a Stacks population summary statistics file.  

1. Locate the "example.sumstats.tsv" file in the BaitsTools/example_data/ and move it to the location where you will execute baitstools.rb.  

2. Specify the population summary statistics file using `-i`. stacks2baits can optionally sort variants into those that vary between populations and those that vary within populations (`-S`) argument. It can then subdivide the within-population variants into those that conform and deviate from Hardy-Weinberg Equilibrium (HWE) (`-H`) according to a chi-squared test. Use `-A` to set the alpha value for this test (supported alphas include 0.1, 0.05, 0.025 and 0.01). Unlike with vcf2baits, the total number of variants (`-t`) is selected *per sorting category*. All other arguments work in the same manner as vcf2baits (see previous tutorial). From our example dataset, we will first sort the variants into those that vary within and between populations and choose a maximum of 5 variants per category. Enter the following command:  

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`ruby baitstools.rb stacks2baits -i example.sumstats.tsv -S -p -t 5`  

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;BaitsTools will produce two tsv files containing up to 5 variants each. "out-betweenpops.tsv"  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;and "out-withinpops.tsv" have the selected between-population and within-population  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;variants respectively.  

3. Sorting by HWE deviation is applied only to the within-population variants. Therefore, using the `-H` option implies `-S`. Enter the following command to sort by HWE deviation with an alpha of 0.025:  

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`ruby baitstools.rb stacks2baits -i example.sumstats.tsv -H -p -t 5 -A 0.025`  

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Three files will be generated this time. "out-inhwe.tsv" contains the selected variants  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;conforming to HWE, while "out-outhwe.txt" has those that deviate from it.  

## 8. References  
1. Campana MG, Parker LD, Hawkins MTR *et al.* (2016) Genome sequence, population history, and pelage genetics of the endangered African wild dog (*Lycaon pictus*). *BMC Genomics*, __17__, 1013.
2. Koepfli K-P, Pollinger J, Godinho R *et al*. (2015) Genome-wide evidence reveals that African and Eurasian golden jackals are distinct species. *Current Biology*, __16__, 2158–2165.
3. Bjornerfeldt S, Webster MT, Vilà C (2006) Relaxation of selective constraint on dog mitochondrial DNA following domestication. *Genome Research*, __16__,990–994.
4. Kim KS, Lee SE, Jeong HW, Ha JH (1998) The complete nucleotide sequence of the domestic dog (*Canis familiaris*) mitochondral genome. *Molecular Phylogenetics and Evolution*, __10__, 210–220.
5. Lindblad-Toh K, Wade CM, Mikkelsen *et al.* (2005) Genome sequence, comparative analysis and haplotype structure of the domestic dog. *Nature*, __438__, 803–819.  
6. Catchen J, Hohenlohe PA, Bassham *et al.* (2013) Stacks: an analysis tool set for population genomics. *Molecular Ecology*, __22__, 3124–3140.  
7. Catchen JM, Amores A, Hohenlohe P *et al.* (2011) Stacks: building and genotyping loci *de novo* from short-read sequences. *G3: Genes, Genomes, Genetics*, __1__, 171–182.  
