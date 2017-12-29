# BaitsTools: software for hybridization capture bait design

Michael G. Campana, 2017  
Smithsonian Conservation Biology Institute  
Contact: campanam@si.edu

## Table of Contents  
1. [Installation](#installation)  
2. [Execution](#execution)  
2a. [Interactive Mode](#interactive-mode)  
2b. [Command-line Mode](#command-line-mode)
3. [Tutorial and Example Data](#tutorial-and-example-data)  
4. [Common Options](#common-options)
5. [Quality Control and Bait Filtration Options](#quality-control-and-bait-filtration-options)
6. [Parameter File](#parameter-file)
7. [Subcommand Arguments](#subcommand-arguments)  
7a. [aln2baits](#aln2baits)  
7b. [annot2baits](#annot2baits)  
7c. [bed2baits](#bed2baits)  
7d. [checkbaits](#checkbaits)  
7e. [stacks2baits](#stacks2baits)  
7f. [tilebaits](#tilebaits)  
7g. [vcf2baits](#vcf2baits)  
8. [Tips and Tricks](#tips-and-tricks)  
9. [Citation](#citation)

## Installation
In a terminal window, execute the following commands:  
`git clone "https://github.com/campanam/BaitsTools/"`  
`cd BaitsTools`  
`chmod +x *.rb`  

Optionally, the Ruby files (.rb files) can be placed within the users $PATH so that they can be executed from any location.  

## Execution  
Executing the command `ruby baitstools.rb` (or `baitstools.rb` if the .rb files are moved to a location within the $PATH) will display the splash screen listing the available subcommands.  

### Interactive Mode  
To enter interactive mode, enter the command `ruby baitstools.rb [subcommand]` without further arguments. This will bring up a series of subcommand-specific interactive prompts. Upon execution, the program will print to the screen the command-line parameters that correspond to the prompt responses for future replication. 

### Command-line Mode  
To use the command-line mode, enter the command `ruby baitstools.rb [subcommand]` with needed arguments. The program will prompt for corrections to parameters that cannot be processed as stated (e.g. missing reference sequences). Upon execution, the program will print to the screen the interpreted command-line parameters for future replication.  

A list of all subcommand-specific arguments (see 'Subcommand Arguments' below) is available with the `-h` or `--help` arguments.  
For instance, to view the help screen for the vcf2baits subcommand, enter:  
`ruby baitstools.rb vcf2baits -h` or `ruby baitstools.rb vcf2baits --help`  

*Entering default values*: For all parameters with defaults, entering the argument flag without specifying a value will cause BaitsTools to assume the default (e.g. entering `-n` will cause BaitsTools to filter by minimum GC content of 30.0%).

## Tutorial and Example Data  
A tutorial and example data are available in the example_data subdirectory of the BaitsTools repository.  

## Common Options  
`-B, --bed`: Output BED file listing the coordinates of baits in comparison to the reference sequences. *Recommended*.  
`-D, --ncbi`: Denotes that FASTA/FASTQ file headers have NCBI-style descriptors after the sample name separated by spaces (e.g. `>{sample_name} {descriptor1} {descriptor2}`).  
`-h, --help`: Print subcommand-specific help to the screen. Use without other arguments (e.g. `ruby baitstools.rb vcf2baits -h`).  
`-v, --version`: Print subcommand version to the screen (which may not correspond with the BaitsTools release version). Use without other arguments (e.g. `ruby baitstools.rb vcf2baits -v`).

## Quality Control and Bait Filtration Options
`-w, --params`: Output table listing bait statistics (GC%, melting temperature, length, etc.). *Highly recommended*.  
`-c, --complete`: Remove candidate baits that are shorter than the requested bait length.  
`-G, --nogaps`: Exclude candidate baits that include gap characters (-).  
`-N, --noNs`: Exclude candidate baits that include unknown bases (Ns).  
`-n, --mingc [VALUE]`: Exclude candidate baits with a GC% below the specified value. Default is 30.0%.  
`-x, --maxgc [VALUE]`: Exclude candidate baits with a GC% above the specified value. Default is 50.0%.  
`-q, --mint [VALUE]`: Exclude candidate baits with melting temperature below the specified value. Default is 0.0°C.  
`-z, --maxt [VALUE]`: Exclude candidate baits with melting temperature above the specified value. Default is 120.0°C.  
`-T, --type [VALUE]`: Hybridization chemistry (one of `DNA-DNA`, `RNA-RNA`, or `RNA-DNA`) for melting temperature estimation. Default is `DNA-RNA`.  
`-s, --na [VALUE]`: Sodium concentration for melting temperature estimation. Default is 0.9M.  
`-f, --formamide [VALUE]`: Formamide concentration for melting temperature estimation. Default is 0.0M.  
`-Q, --meanqual [VALUE]`: Exclude baits with a mean Phred-like base quality below the specified value. Default is 20.0.  
`-M, --minqual [VALUE]`: Exclude baits with any base below the specified Phred-like base quality. Default is 10.  
`-F, --fastascore [VALUE]`: Assume the specified Phred-like base quality for all FASTA sequences (e.g. for mixed FASTQ and FASTA datasets). Default is 0.  

### Parameter File  
The parameter file is a tab-separated file giving bait-specific filtration information. The columns from left to right are:  
* Chromosome (Haplotype): Coordinates: The reference sequence, the haplotype identification (for alternate alleles or alignment data), and the bait coordinates.  
* BaitLength: The length of the generated bait  
* %GC: Bait GC content in percent  
* Tm: Bait melting temperature in Celcius  
* MeanQuality: Mean Phred-like base quality of the generated bait  
* MinQuality: Minimum Phred-like base quality of the generated bait  
* Kept: Whether the bait was retained in the final filtered set. 

## Subcommand Arguments  
### aln2baits  
aln2baits generates baits from a DNA alignment in FASTA or FASTQ format. Bait sets are weighted by variability within each window so that more variable regions have higher coverage.  

`-i, --input [FILE]`: Input alignment file name. Include the path to the file if not in the current directory.  
`-L, --length [VALUE]`: Requested bait length. Default is 120 bp.  
`-O, --offset [VALUE]`: Offset (in bp) between tiled baits. Default is 20 bp.  
`-H, --haplo [VALUE]`: Haplotype definition for bait generation. Entering `haplotype` will cause the program to identify all unique haplotypes within each bait tiling window observed in the data. Alternatively, entering `variant` will cause the program to generate all possible permutations of single nucleotide variants observed within the window.  

### annot2baits  
annot2baits generates baits from an annotation file in GTF or GFF and a corresponding DNA sequence in FASTA or FASTQ format.  

`-i, --input [FILE]`: Input GTF/GFF file name. Include the path to the file if not in the current directory.  
`-r, --refseq [FILE]`: Input reference sequence file name. Include the path to the file if not in the current directory.  
`-L, --length [VALUE]`: Requested bait length. Default is 120 bp.  
`-O, --offset [VALUE]`: Offset (in bp) between tiled baits. Default is 20 bp.  
`-U, --features [FEATURE]`: Comma-separated list of features to extract (e.g. `exon,intron,tRNA`).  

### bed2baits  
bed2baits generates baits from a track file in headerless BED and a a corresponding DNA sequence in FASTA or FASTQ format.  

`-i, --input [FILE]`: Input BED file name. Include the path to the file if not in the current directory.  
`-r, --refseq [FILE]`: Input reference sequence file name. Include the path to the file if not in the current directory.  
`-L, --length [VALUE]`: Requested bait length. Default is 120 bp.  
`-O, --offset [VALUE]`: Offset (in bp) between tiled baits. Default is 20 bp.  

### checkbaits  
checkbaits quality-controls and filters previously generated baits in FASTA or FASTQ format.

`-i, --input [FILE]`: Input sequence file name. Include the path to the file if not in the current directory.  
`-L, --length [VALUE]`: Requested bait length. Default is 120 bp.  

### stacks2baits  
stacks2baits selects variants and generates baits from a Stacks population summary statistics file and a reference sequence.  

`-i, --input [FILE]`: Input Stacks summary tsv file name. Include the path to the file if not in the current directory.  
`-S, --sort`: Sort variants according to variation between or within populations
`-H, --hwe`: Sort variants within populations according to whether they conform to the expectations of Hardy-Weinberg Equilibrium (HWE). This option implies `-S`.  
`-A, --alpha [VALUE]`: Alpha value for the the HWE chi-square test (one of `0.10`,`0.05`, `0.025`, or `0.01`). Default is 0.05.  
`-t, --totalvars [VALUE]`: Total requested variants within each sorted category. For example, without sorting, choosing `-t 30000` will choose a maximum of 30,000 total variants. By choosing to sort variants by within/between populations, it would return a maximum of 30,000 between-population variants and 30,000 within-population variants.  
`-j, --scale`: Scale the maximum number of variants per contig by that contig's length. Overrides the `-m` argument.  
`-m, --maxsnps [VALUE]`: Maximum number of variants per contig. Default is 2.  
`-d, --distance [VALUE]`: Minimum distance (in bp) between variants within a contig. The default is 10,000 bp.  
`-p, --nobaits`: Do not output baits, simply subselect variants. A reference sequence is not required for this analysis.  
`-e, --every`: Output baits for every variant in the input file, skipping subselection.  
`-R, --alt`: Generate baits for alternate alleles. Overrides `-p`.  
`-r, --refseq [FILE]`: Input reference sequence file name. Include the path to the file if not in the current directory.  
`-L, --length [VALUE]`: Requested bait length. Default is 120 bp.  
`-b, --lenbef [VALUE]`: If outputting baits, the number of bases before the variant to include in the bait sequence. Default is 60 bp. Overrides `-L`.  
`-a, --lenaft [VALUE]`: If outputting baits, the number of bases after the variant to include in the bait sequence. Default is 59 bp. Overrides `-L`.  
`-u, --tiling`: If outputting baits, tile bait sequences.  
`-O, --offset [VALUE]`: If tiling baits, offset (in bp) between tiled baits. Default is 20 bp. Overrides `-a` and `-b.`  
`-k, --depth [VALUE]`: If tiling baits, requested baits per variant.  

### tilebaits  
tilebaits generates baits from a list of DNA sequences in FASTA or FASTQ format.

`-i, --input [FILE]`: Input sequence file name. Include the path to the file if not in the current directory.  
`-L, --length [VALUE]`: Requested bait length. Default is 120 bp.  
`-O, --offset [VALUE]`: Offset (in bp) between tiled baits. Default is 20 bp.  

### vcf2baits  
vcf2baits selects variants and generates baits from a VCF file and a reference sequence.  

`-i, --input [FILE]`: Input VCF file name. Include the path to the file if not in the current directory.  
`-V, --varqual [VALUE]`: Minimum variant QUAL score to be included in subselected variants. Default is 30.  
`-t, --totalvars [VALUE]`: Total requested variants.  
`-j, --scale`: Scale the maximum number of variants per contig by that contig's length. Overrides the `-m` argument.  
`-m, --maxsnps [VALUE]`: Maximum number of variants per contig. Default is 2.  
`-d, --distance [VALUE]`: Minimum distance (in bp) between variants within a contig. The default is 10,000 bp.  
`-p, --nobaits`: Do not output baits, simply subselect variants. A reference sequence is not required for this analysis.  
`-e, --every`: Output baits for every variant in the input file, skipping subselection.  
`-R, --alt`: Generate baits for alternate alleles. Overrides `-p`.  
`-r, --refseq [FILE]`: Input reference sequence file name. Include the path to the file if not in the current directory.  
`-L, --length [VALUE]`: Requested bait length. Default is 120 bp.  
`-b, --lenbef [VALUE]`: If outputting baits, the number of bases before the variant to include in the bait sequence. Default is 60 bp. Overrides `-L`.  
`-a, --lenaft [VALUE]`: If outputting baits, the number of bases after the variant to include in the bait sequence. Default is 59 bp. Overrides `-L`.  
`-u, --tiling`: If outputting baits, tile bait sequences.  
`-O, --offset [VALUE]`: If tiling baits, offset (in bp) between tiled baits. Default is 20 bp. Overrides `-a` and `-b.`  
`-k, --depth [VALUE]`: If tiling baits, requested baits per variant.  

## Tips and Tricks
1. While BaitsTools can read wrapped FASTA and FASTQ files, this will slow the program down tremendously. Remove extraneous line breaks in reference sequences before running bait generation. A simple way is using the awk command:  

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }'`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`<input.fasta> > <output.fasta>`  

## Citation  
Please cite:  
Campana, M.G. 2017. *BaitsTools: software for hybridization capture bait design*.
