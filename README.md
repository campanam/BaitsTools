# BaitsTools: software for hybridization capture bait design

Michael G. Campana, 2015-2023  
Smithsonian's National Zoo and Conservation Biology Institute  

BaitsTools is a package to facilitate the design of nucleic acid bait sets for hybridization capture experiments. It can generate RNA and DNA baits from a wide variety of input formats including FASTA/FASTQ sequences and alignments, [Stacks](http://catchenlab.life.illinois.edu/stacks/) population summary statistics files, [PyRAD](http://dereneaton.com/software/pyrad/) and [ipyrad](http://ipyrad.readthedocs.io/) loci files, genome annotations and features (BED/GTF/GFF) and VCF files. BaitsTools provides both a traditional command-line interface with arguments and an interactive interface using text prompts. Please read and cite the accompanying [manuscript](http://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12721/abstract) when using this software.  

## License  
The software is made available under the Smithsonian Institution [terms of use](https://www.si.edu/termsofuse).  

## Table of Contents  
1. [Installation](#installation)  
1a. [Installation using RubyGems and Bundler](#installation-using-rubygems-and-bundler)  
1b. [macOS Installation](#macos-installation)  
1c. [Uninstallation](#uninstallation)  
2. [Execution](#execution)  
2a. [Interactive Mode](#interactive-mode)  
2b. [Command-line Mode](#command-line-mode)  
2c. [GUI Mode](#gui-mode)  
3. [Tutorial and Example Data](#tutorial-and-example-data)  
4. [Common Options](#common-options)
5. [Quality Control and Bait Filtration Options](#quality-control-and-bait-filtration-options)
6. [Parameter File](#parameter-file)
7. [Subcommand Arguments](#subcommand-arguments)  
7a. [aln2baits](#aln2baits)  
7b. [annot2baits](#annot2baits)  
7c. [bed2baits](#bed2baits)  
7d. [blast2baits](#blast2baits)  
7e. [checkbaits](#checkbaits)  
7f. [pyrad2baits](#pyrad2baits)  
7g. [stacks2baits](#stacks2baits)  
7h. [tilebaits](#tilebaits)  
7i. [vcf2baits](#vcf2baits)  
8. [Formula Notes](#formula-notes)  
9. [Tips and Tricks](#tips-and-tricks)  
10. [Bug Reports and Feature Requests](#bug-reports-and-feature-requests)
11. [Citation](#citation)  
12. [References](#References)  

## Installation  
*Ruby >= 2.4.1 is required as of BaitsTools version 1.7.0. BaitsTools versions <= 1.6.8.1 are compatible with Ruby >= 2.0.0.*  

General instructions for installation using RubyGems/Bundler and specific instructions for macOS are provided below. You can test your BaitsTools installation by running the tutorials included in the example_data directory. The archive "tutorial.tgz" includes the expected output of each tutorial. Note that vcf2baits and stacks2baits output will vary slightly due to the random number generator.  

### Installation using RubyGems and Bundler  
The BaitsTools executables can be installed using [RubyGems](https://www.rubygems.org) and [Bundler](https://bundler.io/) (available on most UNIX-like operating systems with [Ruby](https://www.ruby-lang.org) and RubyGems installed). See instructions for macOS below as macOS requires the [Ruby Version Manager](https://rvm.io) to manually install Ruby gems. See the Ruby and RubyGems documentation for installation on other operating systems. Available precompiled gems are listed [here](https://github.com/campanam/BaitsTools/pkgs/rubygems/baitstools).  

First add the GitHub Ruby package repository to your sources, substituting your GitHub username and Personal Access Token for USER and TOKEN respectively:  

`gem sources --add https://USER:TOKEN@rubygems.pkg.github.com/campanam`  

Then download and install the latest gem version by executing the following command in a terminal window:  

`gem install baitstools --version "1.8.1"`  

To manually build and install the gem, execute the following commands in a terminal window:  

`git clone https://github.com/campanam/baitstools`  
`cd baitstools`  
`gem build baitstools.gemspec`  
`gem install baitstools-1.8.1.gem`  

### macOS Installation  
macOS uses a deprecated version of Tcl-Tk as its default Tk framework. For best results, install [ActiveTcl 8.6](https://www.activestate.com/products/activetcl/downloads/) and then reinstall the tk gem (`gem install tk`). Tcl-Tk can also be installed using [Homebrew](https://brew.sh) or [Anaconda](https://anaconda.org/), but the windows are not optimized for these methods.  

_Automatic Installation (macOS):_  
For macOS users, there is a script `osx_install.sh` that will automatically install the [Ruby Version Manager](https://rvm.io/) and the BaitsTools gem. To install using the script: 

`git clone "https://github.com/campanam/BaitsTools/"`  
`cd baitstools`  
`bash osx_install.sh`  

Then restart your terminal.  

_Manual Installation (macOS):_  
Enter the following commands (step annotations are provided after the highlighted text to help debug):

`curl -sSL https://get.rvm.io | bash -s stable`:  Install the Ruby Version Manager.  
`source ~/.rvm/scripts/rvm`: Source the RVM scripts.  
`rvm install 3.2.2`: Install Ruby 3.1.2.  
`rvm --default use 3.2.2`: Set Ruby 3.1.2 as default.  
`git clone https://github.com/campanam/baitstools`: Download the BaitsTools repository.  
`cd baitstools`: Enter the baitstools directory.  
`gem build baitstools.gemspec`: Build the BaitsTools gem.  
`gem install baitstools-1.8.1.gem`: Install the BaitsTools gem.  

_macOS Installation Notes:_
1. The Ruby Version Manager uses [Homebrew](https://brew.sh). During installation you may need to give an administrator password and authorization to install/update Homebrew.  
2. macOS does not include gpg for key verification. Although not necessary, gpg can be installed with Homebrew if you wish to verify your Ruby Version Manager installation using the mpapis public key (see [RVM Installation](https://rvm.io/rvm/install)).  

_GUI Installation:_  
The BaitsTools GUI has only been tested on Intel-based macOS. It is not compatible with Apple Silicon-based machines and has not been extensively tested on other operating systems.  
To install BaitsTools including the GUI, compile the gem manually adding the 'gui' option:  
`git clone https://github.com/campanam/baitstools`  
`cd baitstools`  
`gem build baitstools.gemspec gui`  
`gem install baitstools-1.8.1`  

### Uninstallation  
You can uninstall the BaitsTools gem using:  

`gem uninstall baitstools`  

## Execution  
*As of BaitsTools 1.7.4, the executable is now 'baitstools' rather than 'baitstools.rb'.*  

Executing the command `baitstools` will display the splash screen listing the available subcommands.  

### Interactive Mode  
To enter interactive mode, enter the command `baitstools [subcommand]` without further arguments. This will bring up a series of subcommand-specific interactive prompts. Upon execution, the program will print to the screen the command-line parameters that correspond to the prompt responses for future replication.  

### Command-line Mode  
To use the command-line mode, enter the command `baitstools [subcommand]` with needed arguments. The program will prompt for corrections to parameters that cannot be processed as stated (e.g. missing reference sequences). Upon execution, the program will print to the screen the interpreted command-line parameters for future replication.  

A list of all subcommand-specific arguments (see 'Subcommand Arguments' below) is available with the `-h` or `--help` arguments. For instance, to view the help screen for the vcf2baits subcommand, enter:  

`baitstools vcf2baits -h` or `baitstools vcf2baits --help`  

*Entering default values*: For all parameters with defaults, entering the argument flag without specifying a value will cause BaitsTools to assume the default (e.g. entering `-n` will cause BaitsTools to filter by minimum GC content of 30.0%).  

### GUI Mode
*As of BaitsTools 1.7.4, the executable is now 'baitstoolsgui' rather than 'baitstoolsgui.rb'.*  

Enter the command: `baitstoolsgui`  

Please note that the GUI has only been tested on macOS running Intel processors and may not work well on other operating systems.  

## Tutorial and Example Data  
A tutorial and example data are available in the example_data subdirectory of the BaitsTools repository. The ipyrad.loci file is from the [ipyrad tutorial documentation](https://ipyrad.readthedocs.io/en/latest/tutorial_intro_cli.html) [1].  

## Common Options  
`-o, --outprefix [VALUE]`: Output file prefix. Default is `out`.  
`-Z, --outdir [VALUE]`: Output directory name. Default is current directory.  
`-l, --log`: Output detailed log including subcommand-specific summary information such as retained loci, total selected SNPs, SNP coverage after filtration, etc. *Highly recommended*  
`-B, --bed`: Output BED file listing the absolute coordinates of baits in comparison to the reference sequences. Unless otherwise specified using the #bed annotation (see [Sequence Annotations](#sequence-annotations)), the default starting position of input FASTA/FASTQ sequences is assumed to be position 1 (BED coordinate 0). *Recommended*.   
`-E, --rbed`: Output BED file listing the coordinates of baits relative to input sequences. For example, if a gene is at neighboring positions 570-900 of a reference genome, a bait at positions 570-689 will have BED coordinates 0-120 rather than 569-689.  
`--shuffle`: Shuffle the last bait forward to compensate for reaching the end of a contig  
`-D, --ncbi`: Denotes that FASTA/FASTQ file headers have NCBI-style descriptors after the sample name separated by spaces (e.g. `>{sample_name} {descriptor1} {descriptor2}`).  
`--phred64`: Qualities are encoded in phred64 rather than phred33.  
`-C, --collapse`: Collapse ambiguity codes to a single nucleotide. One nucleotide will be chosen randomly from the possible nucleotides (e.g. either A or G will be selected for a R site). If `-N` is specified, N sites will be treated as missing data and filtered out as specified. Otherwise, N sites will be treated as matching any nucleotide (A, G, C, T(U)).  
`-Y, --rna`: Output bait sequences as RNA rather than DNA.  
`-R, --rc`: Output reverse-complemented baits.  
`-G, --gaps [VALUE]`: Strategy to handle baits that include gap characters (-) (one of `include`, `exclude`,  `extend`). `include` keeps all baits sequences with gaps. `exclude` filters out all baits with gaps. `extend` attempts to extend baits to complete length while removing gap characters.  *WARNING: extended baits will have BED coordinates corresponding to the uncorrected bait sequence.* Default is `include`.  
`--altbaits [VALUES]`: Comma-separated list of additional lengths of baits to generate for same targets. If using checkbaits, this option will truncate previously generated baits by removing bases from 3' end. All alternate bait lengths will be filtered under the same parameters specified for the primary bait length.  
`-5, --5prime [VALUE]`: Sequence to addend to 5' end of baits.  
`-3, --3prime [VALUE]`: Sequence to addend to 3' end of baits.  
`--fillin [VALUE]`: Fill in baits shorter than requested length with specified sequence repeat motif.  
`-X, --threads [VALUE]`: Number of threads. Default is 1.  
`--rng [VALUE]`: Random number seed. Default uses system entropy.  
`--gzip`: Gzip output files.  
`-h, --help`: Print subcommand-specific help to the screen. Use without other arguments (e.g. `ruby baitstools.rb vcf2baits -h`).  
`-v, --version`: Print subcommand version to the screen (which may not correspond with the BaitsTools release version). Use without other arguments (e.g. `ruby baitstools.rb vcf2baits -v`).

## Quality Control and Bait Filtration Options  
`-w, --params`: Output table listing bait statistics (GC%, melting temperature, length, etc.). *Highly recommended*.  
`--disable-lc`: Disable slow linguistic complexity calculations for parameters table  
`-c, --complete`: Remove candidate baits that are shorter than the requested bait length.  
`-N, --noNs`: Exclude candidate baits that include unknown bases (Ns).  
`--noaddenda`: Exclude 5' and 3' addended sequences from bait parameter calculations.  
`-n, --mingc [VALUE]`: Exclude candidate baits with a GC% less than the specified value. Default is 30.0%.  
`-x, --maxgc [VALUE]`: Exclude candidate baits with a GC% greater than the specified value. Default is 50.0%.  
`-q, --mint [VALUE]`: Exclude candidate baits with melting temperature below the specified value. Default is 0.0°C.  
`-z, --maxt [VALUE]`: Exclude candidate baits with melting temperature above the specified value. Default is 120.0°C.  
`-T, --type [VALUE]`: Hybridization chemistry (one of `DNA-DNA`, `RNA-RNA`, or `RNA-DNA`) for melting temperature estimation. Default is `RNA-DNA`.  
`-s, --na [VALUE]`: Sodium concentration for melting temperature estimation. Default is 0.9M.  
`-f, --formamide [VALUE]`: Formamide concentration for melting temperature estimation. Default is 0.0M.  
`-K, --maxmask [VALUE]`: Exclude baits with percentage of masked bases greater than the specified value. *NOTE: BaitsTools considers all bases in lower-case as masked sequence*. Default is 25.0%.  
`-J, --maxhomopoly [VALUE]`: Exclude baits with homopolymers longer than the specified value. Default is 4.  
`-y, --minlc [VALUE]`: Exclude baits with linguistic complexity scores lower than the specified value. Default is 0.9.  
`-Q, --meanqual [VALUE]`: Exclude baits with a mean Phred-like base quality below the specified value. Default is 20.0.  
`-M, --minqual [VALUE]`: Exclude baits with any base below the specified Phred-like base quality. Default is 10.  
`-F, --fastascore [VALUE]`: Assume the specified Phred-like base quality for all FASTA sequences (e.g. for mixed FASTQ and FASTA datasets). Default is 0.  

### Parameter File  
The parameter file is a tab-separated file giving bait-specific filtration information. The columns from left to right are:  
* Chromosome (Haplotype): Coordinates: The reference sequence, the haplotype identification (for alternate alleles or alignment data), and the bait coordinates.  
* BaitLength: The length of the generated bait  
* GC%: Bait GC content in percent  
* Tm: Bait melting temperature in Celcius  
* Masked%: Percent of bait masked  
* MaxHomopolymer: Length of longest homopolymer run  
* SeqComplexity: Bait linguistic complexity  
* MeanQuality: Mean Phred-like base quality of the generated bait  
* MinQuality: Minimum Phred-like base quality of the generated bait  
* Ns: Whether the bait included Ns  
* Gaps: Whether the bait included gaps  
* Kept: Whether the bait was retained in the final filtered set. 

### Sequence Annotations  
BaitsTools reserves the "#" character for program-specific FASTA/FASTQ annotations. Append these annotations to end of sequence headers. The following annotations are currently available:  
`#circ`: Appending this annotation denotes that a sequence is circular. Otherwise the sequence is assumed to be linear.  
`#bed`: Appending this annotation changes the absolute BED starting coordinate from the default of 0 to the specified coordinate (e.g. "#bed80" changes the starting coordinate to 80). This annotation only affects the [tilebaits](#tilebaits) and [aln2baits](#aln2baits) subcommands. If the `variant` definition is used in aln2baits, the absolute BED coordinates will be given relative to the first sequence in the alignment.  
`#loc`: Appending this annotation determines the locus alignment to which a sequence belongs (e.g. "#locSLCA42" assigns the sequence to locus SLCA42). This annotation allows [aln2baits](#aln2baits) to process alignments from multiple loci simultaneously.  

## Subcommand Arguments  
### aln2baits  
aln2baits generates baits from a DNA alignment in FASTA or FASTQ format. Bait sets are weighted by variability within each window so that more variable regions have higher coverage. Sequences can be assigned to locus alignments using the "loc" [sequence annotation](#sequence-annotations).  

`-i, --input [FILE]`: Input alignment file name. Include the path to the file if not in the current directory.  
`-L, --length [VALUE]`: Requested bait length. Default is 120 bp.  
`-O, --offset [VALUE]`: Offset (in bp) between tiled baits. Default is 60 bp.  
`-H, --haplo [VALUE]`: Alignment window haplotype definition (`haplotype` or `variant`). `haplotype` will cause the program to identify all unique haplotypes within each bait tiling window observed in the data. `variant` will cause the program to generate random permutations of single nucleotide variants observed within the window. Default is `haplotype`.  
`--maxvars [VALUE]`: Maximum number of variant permutations to retain within each alignment window when using the `variant` haplotype definition. Default is 24.  

### annot2baits  
annot2baits generates baits from an annotation file in GTF or GFF and a corresponding DNA sequence in FASTA or FASTQ format.  

`-i, --input [FILE]`: Input GTF/GFF file name. Include the path to the file if not in the current directory.  
`-r, --refseq [FILE]`: Input reference sequence file name. Include the path to the file if not in the current directory.  
`-P, --pad [VALUE]`: Length to pad beginning and ending of extracted regions. Padding the annotation coordinates helps ensure that baits cover the entire region evenly. Default is 0 bp.  
`-L, --length [VALUE]`: Requested bait length. Default is 120 bp.  
`-O, --offset [VALUE]`: Offset (in bp) between tiled baits. Default is 60 bp.  
`-U, --features [FEATURE]`: Comma-separated list of features to extract (e.g. `exon,intron,tRNA`).  

### bed2baits  
bed2baits generates baits from a track file in BED or GATK/Picard interval list format and a a corresponding DNA sequence in FASTA or FASTQ format.  

`-i, --input [FILE]`: Input BED or interval list file name. Include the path to the file if not in the current directory.  
`--list [VALUE]`: BED or interval list file format (bed, GATK, or Picard). Default is BED.  
`-r, --refseq [FILE]`: Input reference sequence file name. Include the path to the file if not in the current directory.  
`-P, --pad [VALUE]`: Length to pad beginning and ending of extracted regions. Padding the BED coordinates helps ensure that baits cover the entire region evenly. Default is 0 bp.  
`-L, --length [VALUE]`: Requested bait length. Default is 120 bp.  
`-O, --offset [VALUE]`: Offset (in bp) between tiled baits. Default is 60 bp.  

### blast2baits  
blast2baits generates baits from a BLAST hit tabular file and a corresponding DNA sequence in FASTA or FASTQ format.  

`-i, --input [FILE]`: Input BLAST hit table. Include the path to the file if not in the current directory.  
`-r, --refseq [FILE]`: Input reference sequence file name. Include the path to the file if not in the current directory.  
 `-P, --pad [VALUE]`: Length to pad beginning and ending of extracted regions. Padding the BED coordinates helps ensure that baits cover the entire region evenly. Default is 0 bp.  
`-L, --length [VALUE]`: Requested bait length. Default is 120 bp.  
`-O, --offset [VALUE]`: Offset (in bp) between tiled baits. Default is 60 bp.  
`--percid [VALUE]`: Minimum percent identity to include BLAST hit. Default is 0.0%.  
`--blastlen [VALUE]`: Minimum length of BLAST hit in nucleotides. Default is 1 bp.  
`--evalue [VALUE]`: Maximum E-value to include BLAST hit. Default is 0.1.  

### checkbaits  
checkbaits quality-controls and filters previously generated baits in FASTA or FASTQ format.

`-i, --input [FILE]`: Input sequence file name. Include the path to the file if not in the current directory.  
`--inbed [FILE]`: Optional BED file corresponding to baits in the sequence file. BED file will be filtered alongside baits sequence file.  
`-L, --length [VALUE]`: Requested bait length. Default is 120 bp.  

### pyrad2baits  
pyrad2baits selects variants and generates baits from a PyRAD/ipyrad loci file. It can either treat loci as FASTA alignments using the `alignment` strategy or as variant calls using the `SNPs` or `informative` strategies.  

`-i, --input [FILE]`: Input LOCI file  
`-L, --length [VALUE]`: Requested bait length. Default is 120 bp.  
`-O, --offset [VALUE]`: Base pair offset between tiled baits. Default is 60 bp.  
`-I, --minind [VALUE]`: Minimum number of individuals to include locus. Default is 1.  
`-W, --strategy [VALUE]`: Strategy to generate baits from loci (`alignment`, `SNPs`, or `informative`). `alignment` treats the individual loci as FASTA alignments and passes the alignments to [aln2baits](#aln2baits) to generate weighted alignments. `SNPs` and `informative` select and generate baits for identified variable sites. `SNPs` includes all identified sites, whereas `informative` includes only phylogenetically informative sites. Default is `alignment`.  
`-H, --haplo [VALUE]`: If using `alignment` strategy, alignment window haplotype definition (`haplotype` or `variant`). `haplotype` will cause the program to identify all unique haplotypes within each bait tiling window observed in the data. `variant` will cause the program to generate random permutations of single nucleotide variants observed within the window. Default is `haplotype`.  
`--maxvars [VALUE]`: Maximum number of variant permutations to retain within each alignment window when using the `variant` haplotype definition. Default is 24.  
`--uncollapsedref`: If using `SNPs` or `informative` strategies, choose a random reference sequence and keep ambiguities for each locus.  
`-a, --alt`: If using `SNPs` or `informative` strategies, generate baits for alternate alleles.  
`-t, --totalvars [VALUE]`: If using `SNPs` or `informative` strategies, total requested variants. Default is 30,000.  
`-m, --maxsnps [VALUE]`: If using `SNPs` or `informative` strategies, maximum number of SNPs per locus. Default is 1.  
`-d, --distance [VALUE]`: If using `SNPs` or `informative` strategies, minimum distance between variants within a locus. Default is 100 bp.  
`-k, --depth [VALUE]`: If using `SNPs` or `informative` strategies, requested tiled baits per variant. Default is 1.  

### stacks2baits  
stacks2baits selects variants and generates baits from a Stacks population summary statistics file and a reference sequence.  

*NOTE: Stacks uses 0-based indexing, while BaitsTools uses 1-based indexing. Variants coordinates will thus differ by 1 between the input/output Stacks TSV files and the remaining BaitsTools output.*  

`-i, --input [FILE]`: Input Stacks summary TSV file name. Include the path to the file if not in the current directory.  
`-S, --sort`: Sort variants according to variation between or within populations
`-H, --hwe`: Sort variants within populations according to whether they conform to the expectations of Hardy-Weinberg Equilibrium (HWE). This option implies `-S`.  
`-A, --alpha [VALUE]`: Alpha value for the the HWE chi-square test. Default is 0.05.  
`-t, --totalvars [VALUE]`: Total requested variants within each sorted category. For example, without sorting, choosing `-t 20000` will choose a maximum of 20,000 total variants. By choosing to sort variants by within/between populations, it would return a maximum of 20,000 between-population variants and 20,000 within-population variants. Default is 30,000.  
`-j, --scale`: Scale the maximum number of variants per contig by that contig's length. Overrides the `-m` argument.  
`-m, --maxsnps [VALUE]`: Maximum number of variants per contig. Default is 2.  
`-d, --distance [VALUE]`: Minimum distance (in bp) between variants within a contig. The default is 10,000 bp.  
`-p, --nobaits`: Do not output baits, simply subselect variants. A reference sequence is not required for this analysis.  
`-e, --every`: Output baits for every variant in the input file, skipping subselection.  
`-a, --alt`: Generate baits for alternate alleles. Overrides `-p`.  
`-r, --refseq [FILE]`: Input reference sequence file name. Include the path to the file if not in the current directory.  
`-L, --length [VALUE]`: Requested bait length. Default is 120 bp.  
`-b, --lenbef [VALUE]`: If outputting baits, the number of bases before the variant to include in the bait sequence. Default is 60 bp.   
`-O, --offset [VALUE]`: Base pair offset between tiled baits. Default is 60 bp.  
`-k, --depth [VALUE]`: Requested tiled baits per variant. Default is 1.  

### tilebaits  
tilebaits generates baits from a list of DNA sequences in FASTA or FASTQ format.

`-i, --input [FILE]`: Input sequence file name. Include the path to the file if not in the current directory.  
`-L, --length [VALUE]`: Requested bait length. Default is 120 bp.  
`-O, --offset [VALUE]`: Offset (in bp) between tiled baits. Default is 60 bp.  

### vcf2baits  
vcf2baits selects variants and generates baits from a VCF file and a reference sequence.  

`-i, --input [FILE]`: Input VCF file name. Include the path to the file if not in the current directory.  
`--taxafile [FILE]`: Balance variants by taxa specified in the specified optional TSV (tab-separated values) file. The TSV file is a two column file with sample name as it appears in the VCF in the first column and population assignment in the second, e.g.:  

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;sample1&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;pop1  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;sample2&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;pop1  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;sample3&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;pop2  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;sample4&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;pop2  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;sample5&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;pop3  

`--taxacount [VALUES]`: Comma-separated list of values for balancing variants by variation category (Order: AllPopulations,BetweenPopulations,WithinPopulations). AllPopulations are those variants for that are variable across all taxa. BetweenPopulations are variants that are homozygous within taxa, but variable across taxa. WithinPopulations are variable within a subset of the taxa (but not across all taxa).  
`--popcategories [VALUES]`: Comma-separated list of maximum number of population-specific variants in order of appearance in taxa TSV file.  
`--previousbaits [FILE]:` Generate baits that complement previously generated baits specified in BED file.  
`-V, --varqual [VALUE]`: Minimum variant QUAL score to be included in subselected variants. Default is 30.  
`-t, --totalvars [VALUE]`: Total requested variants. Default is 30,000.  
`-j, --scale`: Scale the maximum number of variants per contig by that contig's length. Overrides the `-m` argument.  
`-m, --maxsnps [VALUE]`: Maximum number of variants per contig. Default is 2.  
`-d, --distance [VALUE]`: Minimum distance (in bp) between variants within a contig. The default is 10,000 bp.  
`-p, --nobaits`: Do not output baits, simply subselect variants. A reference sequence is not required for this analysis.  
`-e, --every`: Output baits for every variant in the input file, skipping subselection. Overrides -t, -j, -d, -m , -p, --taxafile, --popcategories, --previousbaits.  
`-a, --alt`: Generate baits for alternate alleles. Overrides `-p`.  
`-r, --refseq [FILE]`: Input reference sequence file name. Include the path to the file if not in the current directory.  
`-L, --length [VALUE]`: Requested bait length. Default is 120 bp.  
`-b, --lenbef [VALUE]`: If outputting baits, the number of bases before the variant to include in the bait sequence. Default is 60 bp.  
`-O, --offset [VALUE]`: Base pair offset between tiled baits. Default is 60 bp.  
`-k, --depth [VALUE]`: Requested tiled baits per variant. Default is 1.  

## Formula Notes  
Bait melting temperatures are calculated according to salt-adjusted formulas for molecules longer than 50 nucleotides [2-3] as given in [4-5].  

## Tips and Tricks
1. The "#" character is reserved for BaitsTools annotations of sequence headers. Do not include this character in sequence identifiers.  

2. Sequence characters in lowercase are considered masked by BaitsTools. This only impacts the optional `--maxmask` filter.  

3. While BaitsTools can read wrapped FASTA and FASTQ files, this will slow the program down tremendously. Remove extraneous line breaks in reference sequences before running bait generation. A simple way using awk is available [here](https://stackoverflow.com/questions/15857088/remove-line-breaks-in-a-fasta-file) [6].  

4. If you need to enter a negative value into the command-line interface (e.g. a negative melting temperature), omit the space between the parameter and the value (e.g. `-q-30.0` rather than `-q -30.0`). Including the space will cause a Ruby parsing error.  

5. It is possible to generate baits from paired-end sequences that do not overlap in a single BaitsTools run. Concatenate the two reads (reverse-complementing as appropriate) with a pad of Ns between the two reads. Using the `-N` option will then exclude candidate baits that overlap the N pad. This could be helpful for generating baits from FASTA/FASTQ alignments (e.g. using `aln2baits`) in which some sequences overlap and can be merged, but others have an unsequenced gap between the paired reads. This trick can also be used for unmerged paired-end PyRAD/ipyrad loci in `pyrad2baits`.  

6. As of version 1.5.0, BaitsTools considers files with the terminal extension '.gz' to be gzip compressed. All other file names are assumed to be uncompressed.  

## Bug Reports and Feature Requests  
BaitsTools is a complex program under active development. Bugs and technical issues are inevitable. Please report any issues and associated error reports using the issues template. Feature requests can also be filed using the feature request issue template. Please see the [CONTRIBUTING](CONTRIBUTING.Md) guidelines.  

## Citation  
Please cite:  
Campana, M.G. (2018) BaitsTools: software for hybridization capture bait design. *Molecular Ecology Resources*, __18__, 356-361. doi: [10.1111/1755-0998.12721](http://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12721/abstract).  

## References  
1. Eaton, D., Overcast, I. (2019) Introductory tutorial - CLI. (https://ipyrad.readthedocs.io/en/latest/tutorial_intro_cli.html).
2. Howley, P.M., Israel, M.A., Law, M.F., Martin, M.A. (1979) A rapid method for detecting and mapping homology between heterologous DNAs. Evaluation of polyomavirus genomes. *The Journal of Biological Chemistry*, __254__, 4876-4883.  
3. Sambrook, J.F., Russell, D.W. (eds). (2001) Molecular Cloning: A Laboratory Manual. Cold Spring Harbor Laboratory Press: Cold Spring Harbor, NY.  
4. Kibbe, W.A. (2007) OligoCalc: an online oligonucleotide properties calculator. *Nucleic Acids Res*, __35__, W43-W46.  
5. Kibbe, W.A. (2015) Oligo Calc: Oligonucleotide Properties Calculator. Version 3.27. (http://biotools.nubic.northwestern.edu/OligoCalc.html.)   
6. User 'Johnsyweb' (6 April 2013) Stack Overflow. Remove line breaks in a FASTA file. (https://stackoverflow.com/questions/15857088/remove-line-breaks-in-a-fasta-file.)  

  
