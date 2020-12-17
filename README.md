# BaitsTools: software for hybridization capture bait design

Michael G. Campana, 2017-2020  
Smithsonian Conservation Biology Institute  

BaitsTools is an open-source package to facilitate the design of nucleic acid bait sets for hybridization capture experiments. It can generate RNA and DNA baits from a wide variety of input formats including FASTA/FASTQ sequences and alignments, [Stacks](http://catchenlab.life.illinois.edu/stacks/) population summary statistics files, [PyRAD](http://dereneaton.com/software/pyrad/) and [ipyrad](http://ipyrad.readthedocs.io/) loci files, genome annotations and features (BED/GTF/GFF) and VCF files. BaitsTools provides both a traditional command-line interface with arguments and an interactive interface using text prompts. Please read and cite the accompanying [manuscript](http://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12721/abstract) when using this software.  

## License  
The software is made available under the Smithsonian Institution [terms of use](https://www.si.edu/termsofuse).  

## Table of Contents  
1. [Installation](#installation)  
1a. [Basic Installation](#basic-installation)  
1b. [GUI Installation](#gui-installation)  
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
### Basic Installation  
In a terminal window, execute the following commands:  

`git clone "https://github.com/campanam/BaitsTools/"`  
`cd BaitsTools`  
`chmod +x *.rb`  

Optionally, the Ruby files (.rb files) can be placed within the user’s $PATH so that they can be executed from any location. Depending on your operating system, you may need to change the shebang lines in the scripts (first lines starting with #!) to specify the path of your Ruby executable.  

You can test your BaitsTools installation by running the tutorials included in the example data. The archive "tutorial.tgz" includes the expected output of each tutorial. Note that vcf2baits and stacks2baits output will vary slightly due to the random number generator.  

### GUI Installation  
The GUI requires the Ruby gem 'tk' (typically installed using `gem install tk` on most UNIX-like operating systems with the appropriate RubyGems package installed). macOS requires the [Ruby Version Manager](https://get.rvm.io) to manually install Ruby gems. The GUI requires BaitsTools Ruby files to be in your $PATH. Please note that the GUI has only been tested on macOS and may not work well on other operating systems. Due to its requirement of external dependencies, the GUI can be difficult to install. I have included instructions for automatic and manual installation on macOS. Please report any encountered bugs using the Bug Report issues template.  

_Installation on macOS Mojave:_  
The tk gem will not compile correctly on macOS Mojave using the default Tk framework. Install [ActiveTcl 8.5](https://www.activestate.com/products/activetcl/downloads/) and then follow the steps described below. Please note that the tk gem is not currently compatible with Tcl-Tk 8.6, so you will need to install the 8.5 version.  

_Automatic GUI Installation (macOS):_  
For macOS users, there is a script `osx_install.sh` that will automatically install the scripts and modify your $PATH variable as needed. It will also install the [Ruby Version Manager](https://rvm.io/) and the tk gem. To install using the script: 

`git clone "https://github.com/campanam/BaitsTools/"`  
`cd BaitsTools`  
`bash osx_install.sh`  

Afterwards, relaunch your terminal window.  

_Manual GUI Installation (macOS):_  
Enter the following commands (step annotations are provided after the highlighted text to help debug):

`curl -sSL https://get.rvm.io | bash -s stable`:  Install the Ruby Version Manager.  
`source ~/.rvm/scripts/rvm`: Source the RVM scripts.  
`rvm install 2.4.1`: Install RVM Ruby 2.4.1.  
`rvm --default use 2.4.1`: Set Ruby 2.4.1 as default.  
`gem install tk`: Install the latest 'tk' version.  
`git clone "https://github.com/campanam/BaitsTools/"`: Download BaitsTools.  
`cd BaitsTools`  
`chmod +x *.rb`: Make BaitsTools' scripts executable.  
`sudo mv *.rb /usr/local/bin/`: Move scripts to your $PATH.  

_GUI Installation Notes:_
1. The Ruby Version Manager uses Homebrew. During installation you may need to give an administrator password and authorization to install/update Homebrew.  
2. macOS does not include gpg for key verification. Although not necessary, gpg can be installed with Homebrew if you wish to verify your Ruby Version Manager installation using the mpapis public key (see [RVM Installation](https://rvm.io/rvm/install)).  

## Execution  
Executing the command `ruby baitstools.rb` (or `baitstools.rb` if the .rb files are moved to a location within the $PATH) will display the splash screen listing the available subcommands.  

### Interactive Mode  
To enter interactive mode, enter the command `ruby baitstools.rb [subcommand]` without further arguments. This will bring up a series of subcommand-specific interactive prompts. Upon execution, the program will print to the screen the command-line parameters that correspond to the prompt responses for future replication.  

### Command-line Mode  
To use the command-line mode, enter the command `ruby baitstools.rb [subcommand]` with needed arguments. The program will prompt for corrections to parameters that cannot be processed as stated (e.g. missing reference sequences). Upon execution, the program will print to the screen the interpreted command-line parameters for future replication.  

A list of all subcommand-specific arguments (see 'Subcommand Arguments' below) is available with the `-h` or `--help` arguments. For instance, to view the help screen for the vcf2baits subcommand, enter:  

`ruby baitstools.rb vcf2baits -h` or `ruby baitstools.rb vcf2baits --help`  

*Entering default values*: For all parameters with defaults, entering the argument flag without specifying a value will cause BaitsTools to assume the default (e.g. entering `-n` will cause BaitsTools to filter by minimum GC content of 30.0%).  

### GUI Mode
Enter the command:  

`baitstoolsgui.rb`  

## Tutorial and Example Data  
A tutorial and example data are available in the example_data subdirectory of the BaitsTools repository.  

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
`-H, --haplo [VALUE]`: Alignment window haplotype definition (`haplotype` or `variant`). `haplotype` will cause the program to identify all unique haplotypes within each bait tiling window observed in the data. `variant` will cause the program to generate all possible permutations of single nucleotide variants observed within the window. Default is `haplotype`.  

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
`-L, --length [VALUE]`: Requested bait length. Default is 120 bp.  

### pyrad2baits  
pyrad2baits selects variants and generates baits from a PyRAD/ipyrad loci file. It can either treat loci as FASTA alignments using the `alignment` strategy or as variant calls using the `SNPs` or `informative` strategies.  

`-i, --input [FILE]`: Input LOCI file  
`-L, --length [VALUE]`: Requested bait length. Default is 120 bp.  
`-O, --offset [VALUE]`: Base pair offset between tiled baits. Default is 60 bp.  
`-I, --minind [VALUE]`: Minimum number of individuals to include locus. Default is 1.  
`-W, --strategy [VALUE]`: Strategy to generate baits from loci (`alignment`, `SNPs`, or `informative`). `alignment` treats the individual loci as FASTA alignments and passes the alignments to [aln2baits](#aln2baits) to generate weighted alignments. `SNPs` and `informative` select and generate baits for identified variable sites. `SNPs` includes all identified sites, whereas `informative` includes only phylogenetically informative sites. Default is `alignment`.  
`-H, --haplo [VALUE]`: If using `alignment` strategy, alignment window haplotype definition (`haplotype` or `variant`). `haplotype` will cause the program to identify all unique haplotypes within each bait tiling window observed in the data. `variant` will cause the program to generate all possible permutations of single nucleotide variants observed within the window. Default is `haplotype`.  
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
`--taxafile [FILE]`: Balance variants by taxa specified in optional TSV file  
`--taxacount [VALUES]`: Comma-separated list of values for balancing variants by variation category (Order: AllPopulations,BetweenPopulations,WithinPopulations). AllPopulations are those variants for that are variable across all taxa. BetweenPopulations are variants that are homozygous within taxa, but variable across taxa. WithinPopulations are variable within a subset of the taxa (but not across all taxa).  
`--popcategories [VALUES]`: Comma-separated list of maximum number of population-specific variants in order of appearance in taxa TSV file.  
`-V, --varqual [VALUE]`: Minimum variant QUAL score to be included in subselected variants. Default is 30.  
`-t, --totalvars [VALUE]`: Total requested variants. Default is 30,000.  
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

## Formula Notes  
Bait melting temperatures are calculated according to salt-adjusted formulas for molecules longer than 50 nucleotides [1-2] as given in [3-4].  

## Tips and Tricks
1. The "#" character is reserved for BaitsTools annotations of sequence headers. Do not include this character in sequence identifiers.  

2. Sequence characters in lowercase are considered masked by BaitsTools. This only impacts the optional `--maxmask` filter.  

3. While BaitsTools can read wrapped FASTA and FASTQ files, this will slow the program down tremendously. Remove extraneous line breaks in reference sequences before running bait generation. A simple way is using the awk command [5]:  

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' \`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`<input.fasta> > <output.fasta>`  

4. If BaitsTools is running slowly, try using the Rubinius Ruby compiler (https://rubinius.com/). Rubinius can multithread across multiple processors, while the standard Ruby interpreter only uses one processor (despite the number of requested threads).  

5. If you need to enter a negative value into the command-line interface (e.g. a negative melting temperature), omit the space between the parameter and the value (e.g. `-q-30.0` rather than `-q -30.0`). Including the space will cause a Ruby parsing error.  

6. It is possible to generate baits from paired-end sequences that do not overlap in a single BaitsTools run. Concatenate the two reads (reverse-complementing as appropriate) with a pad of Ns between the two reads. Using the `-N` option will then exclude candidate baits that overlap the N pad. This could be helpful for generating baits from FASTA/FASTQ alignments (e.g. using `aln2baits`) in which some sequences overlap and can be merged, but others have an unsequenced gap between the paired reads. This trick can also be used for unmerged paired-end PyRAD/ipyrad loci in `pyrad2baits`.  

7. As of version 1.5.0, BaitsTools considers files with the terminal extension '.gz' to be gzip compressed. All other file names are assumed to be uncompressed.  

## Bug Reports and Feature Requests  
BaitsTools is a complex program under active development. Bugs and technical issues are inevitable. Please report any issues and associated error reports using the issues template. Feature requests can also be filed using the feature request issue template. Please see the CONTRIBUTING guidelines.  

## Citation  
Please cite:  
Campana, M.G. (2018) BaitsTools: software for hybridization capture bait design. *Molecular Ecology Resources*, __18__, 356-361. doi: [10.1111/1755-0998.12721](http://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12721/abstract).  

## References  
1. Howley, P.M., Israel, M.A., Law, M.F., Martin, M.A. (1979) A rapid method for detecting and mapping homology between heterologous DNAs. Evaluation of polyomavirus genomes. *The Journal of Biological Chemistry*, __254__, 4876-4883.  
2. Sambrook, J.F., Russell, D.W. (eds). (2001) Molecular Cloning: A Laboratory Manual. Cold Spring Harbor Laboratory Press: Cold Spring Harbor, NY.  
3. Kibbe, W.A. (2007) OligoCalc: an online oligonucleotide properties calculator. *Nucleic Acids Res*, __35__, W43-W46.  
4. Kibbe, W.A. (2015) Oligo Calc: Oligonucleotide Properties Calculator. Version 3.27. (http://biotools.nubic.northwestern.edu/OligoCalc.html.)   
5. User 'Johnsyweb' (6 April 2013) Stack Overflow. Remove line breaks in a FASTA file. (https://stackoverflow.com/questions/15857088/remove-line-breaks-in-a-fasta-file.)  

  
