# BaitsTools Change Log  
Michael G. Campana, 2017-2023  
Smithsonian's National Zoo and Conservation Biology Institute  
Contact: campanam@si.edu  

## Table of Contents  
[aln2baits](#aln2baits)  
[annot2baits](#annot2baits)  
[baitslib](#baitslib)  
[baitstools](#baitstools)  
[baitstoolsgui](#baitstoolsgui)  
[bed2baits](#bed2baits)  
[blast2baits](#blast2baits)  
[checkbaits](#checkbaits)  
[osx_install](#osx_install)  
[pyrad2baits](#pyrad2baits)  
[stacks2baits](#stacks2baits)  
[tilebaits](#tilebaits)  
[vcf2baits](#vcf2baits)  
[Deprecated](#deprecated)  

## aln2baits  
### Version 1.8.0  
Added --maxvars option to make variant haplotype definition more efficient/function properly  

### Version 1.7.8  
Fixed bug in --shuffle option that caused infinite loop  
Fixed hanging bug when running variant options  
Discovered -- but did not fix -- that variant permutation was running incorrectly  

### Version 1.7.7  
Fixed Thread.current[:window].var_permutations(aln) glitch using variant haplotype definition  

### Version 1.7.0  
altbaits option  

### Version 1.6.3  
Log outputs total numbers of generated/filtered baits in addition to coverage  

### Version 1.6.0  
5' and 3' addenda  
Sequence motif fill-in  
Correction of bugs/improvement of reversecomp/RNA qualities and ordering  
Used thread bug fix  
Collapsed baits output  

### Version 1.4.0  
Circular sequence handling  

### Version 1.2.1  
Fixed bugs in gap extension  
Used appending for sequence variant generation  

### Version 1.2.0  
Fixed possible issues due to uniq!  

### Version 1.1.0  
Reworking to use new temp file methods  

### Version 1.0.4  
Added ability to shuffle final bait if contig end reached  

### Version 1.0  
Updated parameters table  
Can reverse complement baits  

### Version 0.7  
Absolute/relative BED coordinates  
Gap handling options  
Param file notes Ns and gaps, mask %  
Mask handling  
Detailed log  
Ambiguity handling transferred to get_ambiguity in baitstools  
Ability to handle multiple loci  

### Version 0.6  
Rubinius threading compatibility  
BED start index bug fixed  

### Version 0.5  
Upcasing removed since now read as upcased  
RNA output handling added  
Added multi-threading support  
Improved algorithm of generating baits  

### Version 0.4  
Version constant added to header  
Outputs bed files rather than coordinates tables  

### Version 0.3  
Parameter header updated  
Program now updates to screen when major steps started  
Completely revised Hap_Window.var_permutations to get all possible permutations and allow wobble bases  
Haplotype definition 'variants' is now 'variant' for consistency  
Updated the 'variant' options in aln2baits to use all permutations rather than random selection of original sequences. Output headers now say 'Alignment' since they do not correspond to original sequences.  
Option 'no_indels' moved to general filtration option (now -I rather than -N)  

### Version 0.2  
FASTA headers of generated baits now match the original source sequence  
Corrected double-tabbed output in parameters output  

### Version 0.1  
Preliminary script to generate weighted bait set from a fasta alignment  

## annot2baits  
### Version 1.7.3  
Fixing missing reference sequence bug  

### Version 1.6.8  
improved gz_file_open takes block  

### Version 1.6.6  
Removed eval statement from gz_file_open call  

### Version 1.5.0  
Read/write gzip files  

### Version 1.4.1  
Reduced code redundancy using baitslib functions get_looped_sequences, get_padded faseq  

### Version 1.4.0  
Circular sequence handling for padding that extends beyond sequence break  

### Version 1.2.3  
Region tiling code transferred to baitslib  

### Version 1.1.0  
Reworking to use new temp file methods  
Log gives extracted region chromosome names  
Absolute bed now gives correct chromosome name  

### Version 0.4  
Absolute BED coordinates recorded  
Padding parameter  
Detailed log  

### Version 0.3  
Upcasing removed  
Output directory and prefix can be named  

### Version 0.2  
Version constant added to header  

### Version 0.1  
Preliminary script to generate baits from an annotation file and a reference sequence  

## baitslib  
### Version 1.8.0  
Handling for --maxvars option  

### Version 1.7.4  
mean method uses .sum rather than .reduce(:+)  

### Version 1.7.2  
Fixed a glitch in selectsnps which caused a crash if previous baits specify a chr not found in current VCF  

### Version 1.7.1  
Fixed a glitch in selectsnps which caused a crash if no suitable SNPs were found after taxonomic sorting  

### Version 1.7.0  
Improved selectsnps distance-based SNP deletion  
Previously generated bait complementation in selectsnps  
taxafile/popcategories no longer erroneously overrides every option  
altbaits option required edits of write_files, cat_files  
method checkpop removed  
Fixed cat_files bug in snp_to_baits   
Fixed file output bugs with stacks2baits  

### Version 1.6.8  
improved gz_file_open takes block  
Fa_Seq#make_DNA converts '?' bases to 'N' to prevent downstream bugs  

### Version 1.6.6  
Removed eval statement from gz_file_open call  
Revised gz_file_open to return reader class rather than a string  

### Version 1.6.5  
Added pipe character (|) to reserved characters list in resolve_unix_path  

### Version 1.6.3  
snp_to_baits logs the total numbers of generated/filtered baits in addition to coverage  
mean uses inject for arraying summing (Array#sum not compatible with Ruby 2.0)  

### Version 1.6.2  
Minor code improvements  

### Version 1.6.1  
max_homopolymer bug fix  

### Version 1.6.0  
5' and 3' addenda  
Sequence motif fill-in  
Correction of bugs/improvement of reversecomp/RNA qualities and ordering  
Numerous new methods for improved handling  
Used thread bug fix  
read_fasta bug fix  
Revision of ambiguity handling code  
Collapsed baits output  

### Version 1.5.0  
Read/write gzip files  

### Version 1.4.1  
Reduced code redundancy using baitslib functions get_looped_sequences, get_padded faseq  

### Version 1.4.0  
get_sequence_tags method cleans up redundant code in read_fasta  
Bug fixes to prevent crashes with unnamed sequences, sequences named 'circ' and sequence names starting with 'loc'  
extend_baits can wrap around circular sequences  

### Version 1.3.3  
Multithreading concatenation bug fix introduced by unix path resolution  

### Version 1.3.2  
Bug fix for reserved characters in file paths  

### Version 1.3.0  
RNG seed option  
Interval list option  

### Version 1.2.3  
Fixed dup error for $options.popcategories == nil  
Transferred redundant region tiling code from annot2baits, bed2baits, blast2baits  
New tile_regions code prints error and exits when no suitable region found  
Fixed bug in taxacount handling that allowed negative values  

### Version 1.2.2  
Fixed bug in max_homopolymer for Y bases  

### Version 1.2.0  
Fixed bugs in snp.line causing extra line breaks with write_file  
Taxon-based SNP sorting  
Fixed possible issues due to uniq!  
FASTQ value hash converted to library method  
phred64 encoding  

### Version 1.1.0  
get_command_line and read_fasta use << to improve performance  
File start called with setup_output  
setup_temp method sets up temporary files  
concat_file concatenates temporary files  
selectsnps compatible with temporary files  
snp_to_baits compatible with temporary files  
removal of write_baits  
filter_baits speed improvement using count rather than loops  
Fa_Seq has bedheader for region extraction chromosome names  
filter_baits can now exclude linguistic complexity calculations  
get_command_line updated for blast2baits  

### Version 1.0.4  
pyrad2baits can produce an uncollapsed reference sequence  
Added ability to shuffle final bait if contig end reached  

### Version 1.0.3  
Fixed bug in snp_to_baits that added 2 bp to bait lengths.  
baitstools and baitslib version will now track package version (even if no updates to these scripts) for package clarity. Individual version numbers left only to subcommands/GUI.  

### Version 1.0  
Library restored for general release  
Includes all general methods and classes from baitstools  
filter_baits can filter by linguistic complexity and maximum homopolymer length  
Filter statistics calculated only when needed since new filters much slower than earlier ones  

### Version 0.4  
Added an exception to write_baits so that checkbaits does not output unfiltered baits  
Added @qual to class Fa_Seq to allow it to handle FASTQ files  
Modified read_fasta so that it can read standard format FASTQ  
read_fasta can now separate NCBI-style descriptions from sequence names  
Fixed filter_baits bug that always filtered by bait length completeness  
GC content output in parameters file in % rather than decimal  
filter_baits can filter by quality scores  
filter_baits can remove sequences with Ns and gap characters (-)  
Improved melting temperature calculations to use formamide concentrations and calculate for DNA hybridization  
Added method 'mean' to calculate mean value of an array  
Chromo_SNP updated to allow vcf metadata and store original vcf/stacks lines  
selectsnps now uses Chromo_SNPs for the SNP data to allow metadata  
snps_to_baits uses Chromo_SNPs and improved vcf2baits algorithm  
snps_to_baits now updates to screen when major steps started  
snps_to_baits can use quality scores as necessary. Parameter header updated  
snps_to_baits can apply alternate alleles (including indels) from stacks2baits and vcf2baits  
Fixed bug in snps_to_baits that reread reference sequence every time run in stacks2baits  
New method get_command_line gets complete command-line for user reference  

### Version 0.3  
Added class Chromo_SNP as a container for SNP data  
selectsnps method added (moved from vcf2baits)  
snps_to_baits method added (moved from vcf2baits)  
Modified write_baits to allow for different filestems (but defaulting to $options.infile)  
Fixed bug in params file output that reported raw number of GCes rather than GC%  

### Version 0.2
filter_probes renamed baitslib  
The word 'probe' changed to 'baits' in all instances for clarity  
New method read_fasta handles basic fasta input  
read_fasta removes #circ from appropriate sequence headers  
New class Fa_Seq is container for fasta format  
Fa_Seq is compatible with multi-line sequence fasta files  
New method write_probes handles basic output  

### Version 0.1  
filter_probes definition removed into separate script for access by other scripts  

## baitstools  
### Version 1.8.1  
Made script executable in repository  

### Version 1.8.0  
Handling for --maxvars option  

### Version 1.7.4  
Conversion to a RubyGem  

### Version 1.7.0  
Complement previous baits option  
Interactive mode no longer prompts you to enter linguistic complexity if LC calculation disabled  
vcf2baits/pyrad2baits/stacks2baits tiling offset can now be longer than bait length  
vcf2baits/pyrad2baits/stacks2baits tiles per SNP can now be greater than length of bait/offset  
altbaits option  
gzip of log moved to end of baitstools script  
Use of sum and any? methods for array requires Ruby >= 2.4.1  
blast2baits and annot2baits no longer have erroneous --list option  
method ynq to clean up yes/no questions in interactive prompt  
checkbaits previously made BED filtering (--inbed option)  
fixed bug that caused infinite loop in checkbaits interactive mode  

### Version 1.6.7  
Capitalization of options that select between choices no longer matters  
bed2baits description updated to include interval lists  
bed2baits list option capitalized to BED for consistency  

### Version 1.6.3  
Program version numbers printed to log and upon final command execution  

### Version 1.6.0  
5' and 3' addenda  
Sequence motif fill-in  
--collapse moved to general options  
Ambiguities and reverse complement hashes generated at start-up  

### Version 1.5.0  
Read/write gzip files  

### Version 1.3.0  
RNG seed option  
Interval list option  

### Version 1.2.3  
Fixed popcategories bug when neither pop-specific variants or taxa file specified  

### Version 1.2.1  
Fixed bug requesting population-specific variants when no taxa file specified  
Replaced word SNP to variant in interactive mode  

### Version 1.2.0  
Taxon-based SNP-sorting  
Fixed bugs requesting filters when no baits are going to be output  
phred64 quality encoding  

### Version 1.1.0  
File start called with setup_output  
Reworking to use new temp file methods  
HWE alpha no longer hard-coded  
Updated for blast2baits  
Linguistic complexity calculations can be disabled for fast parameters  

### Version 1.0.4  
pyrad2baits can produce an uncollapsed reference sequence  
Added ability to shuffle final bait if contig end reached  

### Version 1.0.3  
baitstools and baitslib version will now track package version (even if no updates to these scripts) for package clarity. Individual version numbers left only to subcommands/GUI.  

### Version 1.0.1  
Fixed inaccurate help text for --maxmask  

### Version 1.0  
GC% calculated for uncollapsed ambiguities  
General methods and classes moved to restored baitslib  
Linguistic complexity and maximum homopolymer length filters  
Reverse complement baits option  
snp_to_baits can reverse complement  

### Version 0.9  
Tiling offset default changed to 60 bp  
Added pyrad2baits method  
Fixed bug in get_command_line for negative parameter values  
Quality scores limited to possible range of 0 to 93  
TSV capitalized in help text  
write_baits can output relative BED positions  
Detailed log  
Gap handling options  
Param file notes Ns and gaps, mask %  
Mask handling  
Revision of variant tiling parameters  
Alt alleles now -a for intuitiveness  
selectsnps speed improvements  
Ambiguity handling transferred to get_ambiguity in baitstools  
snp_to_baits bait headers more self-explanatory  
read_fasta can handle multiple loci with loc tag  

### Version 0.8  
Rubinius threading compatibility  
Output directory and prefix can be named  
Fixed bug in interactive prompts to ensure formamide entered and no negative Na/Formamide concentrations  
BED start index bug fixed  

### Version 0.7  
-H defaults to haplotype  
Ambiguity collapsing added to bait filtration  
Fa_Seq method 'make_dna' converts uracil residues to thymines for internal consistency  
make_dna upcases all sequences for later analyses  
read_fasta calls make_dna after each sequence read in  
frontend handling for RNA bait output  
RNA output handling for snp_to_baits and collapse_ambiguity  
Basic multithreading added for read_fasta, snp_to_baits  
Frontend control for threading (-X)  

### Version 0.6  
Incorporated baitslib into baitstools main script  
require statements for subcommands now require_relative  
Version constant added to header  
Subcommand script versions now reported by constants to minimize having to manually update  
Updated control program to output BED files rather than proprietary coordinates tables  
snp_to_baits now returns BED format information rather than coordinates tableswrite_baits outputs BED files rather than 'coords.txt' files  
Chromo_SNP now includes a reference allele variable  
Fixed option -A to accept values as needed  
Fixed options -t, -m, -L, -a, -b to avoid crashing when called without specifying a value  
Option -o, --coords is now -B, --bed to reflect updated file type  
Filter option -I, --no_Indels now -G, --nogaps to be more accurate and less confusing  
Option --no_Ns now --noNs for naming consistency  
Option --no_baits now --nobaits for naming consistency  
Option -G, --features now -U,--features  
--totalsnps now --totalvars to be more accurate  
get_command_line updated to reflect options revisions  
Updated snp_to_baits to properly apply indels in VCF files  
Help text and prompts updated for accuracy  
Fixed interactive checking of alpha values from command-line interface  
Fixed unnecessary request for NCBI headers when no_baits specified  

### Version 0.5  
Added checkbaits method  
Added annot2baits method  
Added code for new command-line options  
Removed non-functional code in baitstools main help output  
Help/interactive mode text cleaned up  
baitstools builds a FASTQ translation hash at the start-up to speed up FASTQ numeric translation in subprograms  
Changed the vcf2baits/stacks2baits default to outputting baits  
Options -e and -R (new) override -p (no baits) to fix conflicts  
Choosing parameters forces the required filtration options  
Commented code more thoroughly to describe the variables  
Program outputs the complete command line to screen on start-up  
Program outputs 'Program complete' on completion  

### Version 0.4  
Added stacks2baits method  
Fixed bug in selectsnps algorithm call that asked to filter baits even when bait output not selected  
Set tilebaits, coords2baits, and aln2baits to have $options.baits = true since these algorithms only output baits

### Version 0.3  
Added aln2baits to general help output for baitstools  
Fixed no_tiling bug for command-line -N option  

### Version 0.2  
First version but to make consistent with package numbering  
Allowed tiling offsets to be greater than bait length for everything except selectsnps (since it works differently)  
Changed bait length to -L to make consistent across scripts  
Allowed selectsnps to directly determine bait length (but lenaft and lenbef override if specified)  
Changed tiling offset to -O from -v to allow -v to be version  
New front-end for baitstools package (hence inconsistent version number)  
The word 'probe' changed to 'baits' in all instances for clarity  
Set default for tiling offset as 20 bp (from 60 for select_snps and 25 for tile_probes)  

### baitstools.gemspec  
### Version 1.8.1  
GUI must be specifically requested during build  

## baitstoolsgui  
### Version 1.8.1  
Made script executable in repository  
Updated P&H easter egg  

### Version 1.8.0  
Handling for --maxvars option  

### Version 1.7.5  
Fixed bug calling baitstools.rb rather than updated baitstools executable  

### Version 1.7.4  
Conversion to a RubyGem  

### Version 1.7.0  
Handling for previousbaits option  
Every option for vcf2baits now bypasses taxafile/previousbaits screen  
vcf2baits/pyrad2baits/stacks2baits tiling offset can now be longer than bait length  
vcf2baits/pyrad2baits/stacks2baits tiles per SNP can now be greater than length of bait/offset  
--inbed option handling for checkbaits  
--altbaits option handling  

### Version 1.6.7  
bed2baits list option capitalized to BED for consistency  

### Version 1.6.0  
Sequence addenda/fill-in handling added  

### Version 1.5.0  
Read/write gzip files  

### Version 1.3.2  
Bug fix for reserved characters in file paths  

### Version 1.3.0  
RNG seed option  
Interval list option  

### Version 1.2.0  
Fixed RC option bug  
phred64 option  
Taxon-based SNP sorting  

### Version 1.1.1  
Fixed evalue window bug  
Fixed start_baitstools reference sequence/padding parameter for blast2baits bug  

### Version 1.1.0  
start_baitstools uses << to improve performance  
blast2baits handling added  
HWE alpha no longer hard-coded  
disable-lc option added  

### Version 1.0.4  
pyrad2baits can produce an uncollapsed reference sequence  
Added ability to shuffle final bait if contig end reached  

### Version 1.0  
Linguistic complexity and maximum homopolymer length filters  
Reverse complement option  

### Version 0.2  
Tiling offset default changed to 60 bp  
Fixed bug in start_baitstools for negative parameter values  
Quality scores limited to possible range of 0 to 93  
Added relative BED option  
Detailed log  
Gap handling options  
Mask handling  
Revision of variant tiling parameters  
pyrad2baits added  

### Version 0.1  
Basic GUI compatible with baitstools 0.8  

## bed2baits  
### Version 1.6.8  
improved gz_file_open takes block  

### Version 1.6.7  
Headed BEDs are now readable  
Fixed bug that did not process whitespace separated values correctly  
Revision of list options to account for capitalization  

### Version 1.6.6  
Removed eval statement from gz_file_open call  

### Version 1.5.0  
Read/write gzip files  

### Version 1.4.1  
Reduced code redundancy using baitslib functions get_looped_sequences, get_padded faseq  

### Version 1.4.0  
Circular sequence handling for padding that extends beyond sequence break  

### Version 1.3.0  
Interval list option  

### Version 1.2.3  
Region tiling code transferred to baitslib  

### Version 1.2.2  
chromo bug fix for seq.bedheader  

### Version 1.1.0  
Reworking to use new temp file methods  
Log gives extracted region chromosome names  
Absolute bed now gives correct chromosome name  

### Version 0.5  
Absolute BED coordinates recorded  
Padding parameter  
Detailed log  

### Version 0.4  
Output directory and prefix can be named  
BED start index bug fixed  

### Version 0.3  
BaitsTools-specific coordinates table was more-or-less BED format with excess symbols. Converted the coordinates table to BED format.  
Version constant added to header  

### Version 0.2  
Updated algorithm to use quality scores as necessary  
Program now updates to screen when major steps started  
Popvar @alleles replaces @pnuc and @qnuc, which were otherwise unused  

### Version 0.1  
Preliminary script to turn a coordinates table and a reference sequence into baits  

## blast2baits  
### Version 1.6.8  
improved gz_file_open takes block  

### Version 1.6.6  
Removed eval statement from gz_file_open call  

### Version 1.6.0  
Correction of bugs/improvement of reversecomp/RNA qualities and ordering  

### Version 1.5.0  
Read/write gzip files  

### Version 1.4.2  
Reverse-complement quality bug fix  

### Version 1.4.1  
Reduced code redundancy using baitslib functions get_looped_sequences, get_padded faseq  

### Version 1.4.0  
Circular sequence handling for padding that extends beyond sequence break  

### Version 1.3.1  
Prevent crash for missing sequence in reference file  
Notifications of coordinate adjustments for out-of-bounds coordinates  

### Version 1.2.3  
Region tiling code transferred to baitslib  

### Version 1.1.0  
Preliminary script to generate baits from BLAST hit tables and a reference sequence  

## checkbaits  
### Version 1.7.6  
Fixed omission of word bp from ** Filtering XX baits **  

### Version 1.7.0  
checkbaits previously made BED filtering (--inbed option)  
--altbaits option  

### Version 1.6.0  
5' and 3' addenda  
Sequence motif fill-in  
Code-update and bug fixes  
Used thread bug fix  
Collapsed baits output  

### Version 1.1.0  
Reworking to use new temp file methods  

### Version 1.0  
Updated parameters table  
Can reverse complement baits  

### Version 0.5  
Param file notes Ns and gaps, mask %  
Mask handling  
Detailed log  

### Version 0.4  
Rubinius threading compatibility  

### Version 0.3  
RNA output handling added  
Added multi-threading support  

### Version 0.2  
Version constant added to header  

### Version 0.1  
Preliminary script to filter predefined baits through quality filters  

## osx_install  
### Version 1.8.1  
Updated Ruby version to latest (3.2.2)  
Installs latest baitstools gem (1.8.1)  

### Version 1.8.0  
Installs latest baitstools gem (1.8.0)  

### Version 1.7.8  
Installs latest baitstools gem (1.7.8)  

### Version 1.7.7  
Installs latest baitstools gem (1.7.7)  

### Version 1.7.5  
Updated Ruby version to latest (3.1.2)  
Installs latest baitstools gem (1.7.5)  

### Version 1.7.4  
Updated Ruby version to latest (3.0.0)  
Installs latest baitstools gem (1.7.4)  

### Version 1.6.8  
Updated Ruby version to latest (2.7.2)  
Installs 'shell' gem due to removal from bundle in Ruby 2.7  
Update for zsh as macOS default shell  

### Version 1.6.5  
Updated Ruby version to latest (2.6.5)  

### Version 1.5.0  
Updated Ruby version to latest (2.6.4)    

### Version 1.0.3  
Changes now recorded in change log.  
Versioning added to header of script  

### Version 1.0.2  
RVM default bug fixed  

## pyrad2baits  
### Version 1.7.0  
improved gz_file_open takes block  
--altbaits option  

### Version 1.6.6  
Removed eval statement from gz_file_open call  

### Version 1.6.0  
Diplotype sequence bug fix  
Removal of unused get_variants method from Pyrad_Locus  

### Version 1.6.0  
Updated to use $ambig_hash  

### Version 1.5.0  
Read/write gzip files  

### Version 1.1.0  
Reworking to use new temp file methods  

### Version 1.0.4  
pyrad2baits can produce an uncollapsed reference sequence  

### Version 0.1  
Preliminary script to generate baits from PyRAD/ipyrad loci files  

## stacks2baits  
### Version 1.7.0  
--altbaits option  
Fixed bug in which filtered TSV printed even if not filtering  
Files labelled 'all' rather than 'betweenpops' when not sorting  
Fixed file output bugs  

### Version 1.6.8  
improved gz_file_open takes block  

### Version 1.6.6  
Removed eval statement from gz_file_open call  

### Version 1.5.0  
Read/write gzip files  

### Version 1.1.0  
Reworking to use new temp file methods  
HWE test uses floats rather than integers for precision  
HWE alpha no longer hard-coded (CDF calculated on fly using chi_cum_prob)  

### Version 0.4  
Output directory and prefix can be named  
Detailed log  

### Version 0.3  
Version constant added to header  

### Version 0.2  
Bait output now default  
Fixed possible bug in Popvar.in_hwe? method that used local variable rather pfreq rather than @pfreq  
Chromo_SNPs passed directly to selectsnps and snps_to_baits  
stacks2baits outputs filtered stacks summary tsvs (new method write_stacks)  
Improved method to read tsvs using array splitting rather than a cumbersome loop  
Program now updates to screen when major steps started  
Fixed bug in that reread reference sequence every time snps_to_baits run in stacks2baits  

### Version 0.1  
Preliminary script to turn a Stacks summary tsv file and a reference sequence into baits  

## tilebaits  
### Version 1.7.8  
Fixed bug in --shuffle option that caused infinite loop  

### Version 1.7.0  
altbaits option  

### Version 1.6.3  
Log now notes total numbers of baits generated/filtered and the mean bait coverages across loci (rough mean: does not account for differences in source sequence lengths)  

### Version 1.6.0  
5' and 3' addenda  
Sequence motif fill-in  
Correction of bugs/improvement of reversecomp/RNA qualities and ordering  
Used thread bug fix  
Collapsed baits output  

### Version 1.4.0  
Baits can extend multiple times around ultra-short reference sequences  
Fixed bug in quality scores for FASTQ files  
Fixed bug that permitted going off sequence end of short linear sequences when shuffle option turned on.  

### Version 1.1.0  
Reworking to use new temp file methods  

### Version 1.0.4  
Added ability to shuffle final bait if contig end reached  

### Version 1.0  
Updated parameters table  
Can reverse complement baits  

### Version 0.7  
Absolute/relative BED coordinates handled  
Fixed :seqst-1 bug  
Gap handling options  
Param file notes Ns and gaps, mask %  
Mask handling  
Detailed log  

### Version 0.7  
Rubinius threading compatibility  
BED start index bug fixed  

### Version 0.6  
RNA output handling added  
Added multi-threading support  

### Version 0.5  
Returns BED format information rather than coordinates tables  
Version constant added to header  

### Version 0.4  
Parameter header updated  
Updated algorithm to use quality scores as necessary  
Program now updates to screen when major steps started  

### Version 0.3  
Corrected double-tabbed output in parameters output  

### Version 0.2  
tilebaits can now take internally passed input rather than just from command-line  
Input now controlled by baitstools frontend. Converted to method tilebaits.  
renamed tilebaits for consistency  
The word 'probe' changed to 'baits' in all instances for clarity  
Input now controlled by baitslib read_fasta function  
Output now controlled by baitslib write_probes function  

### Version 0.1  
Preliminary script to divide genic/genomic sequences into tiled baits

## vcf2baits  
### Version 1.7.0  
--altbaits option  
Fixed glitches in Gzipping VCF output  

### Version 1.6.8  
improved gz_file_open takes block  

### Version 1.6.6  
Removed eval statement from gz_file_open call  

### Version 1.6.4  
Revised output VCF header to include baitslib version  

### Version 1.5.0  
Read/write gzip files  

### Version 1.3.1  
fixed bug in header column causing extra line break with write_file  

### Version 1.2.0  
fixed bug in snp.line causing extra line breaks with write_file  
Taxon-based SNP sorting  

### Version 1.1.0  
Reworking to use new temp file methods  

### Version 0.10  
Output directory and prefix can be named  
Detailed log  

### Version 0.9  
Now records reference allele to allow proper application of indel variants  
Version constant added to header  
VCF header modification now uses version constants  

### Version 0.8  
Bait output now default  
Now uses Chromo_SNPs  to allow metadata  
Program now updates to screen when major steps started  
Revised algorithm to store original lines when reading vcf. No longer needs to reread vcfs when writing output  
Output vcfs now have BaitsTools command lines written into them for reference  
Added option to filter vcf variants by QUAL scores  

### Version 0.7  
Method involved in actual snp selection moved to baitslib as method selectsnps  
Method involved in bait generation moved to baitslib as method snp_to_baits  
Renamed to vcf2baits to be more consistent with nomenclature of other subprograms  

### Version 0.6  
Input now controlled by baitstools frontend. Converted to method selectsnps  
Renamed selectsnps for consistency  
The word 'probe' changed to 'baits' in all instances for clarity  
Most output now controlled by baitslib write_probes function  
Reference sequence now read in through read_fasta function in baitslib  
Can now tile/get complete baits across circular genomes  
Corrected probe locations for 0-based indexing of Ruby  

### Version 0.5  
Incorporated into BaitsTools package  
filter_probes definition removed into separate script for access by other scripts  

### Version 0.4  
Added the option to produce tiled baits for SNPs  
Added the option to output coordinates tables for the probes  
Added the option to scale number of SNPs per contig by contig length  
Added the option to filter by bait melting temperature  
Added -e option to output probe sequences without selecting SNPs  
Added the ability to output a table of bait filtering statistics  
Made filtering option output filtered VCF as well as selected VCF  
Changed GC filtering defaults and allowed this filter to be skipped if desired  
Fixed 'ouput' typo  
Removed unnecessary ::Version variable  

### Version 0.3  
GC content is now floating point  
Altered the code to use OptionParser/ostruct.  
Set defaults for all major variables  
Fixed a bug that terminated reference sequence chromosome names only at spaces rather than line breaks  

### Version 0.2  
Filters probes based on GC content, length.  
Improved handling of Y/N prompts.  

### Version 0.1  
Selects SNPs from VCFs  
Outputs probes based on reference sequence

# Deprecated  

## coords2baits  
Now bed2baits  

## filter_probes  
Incorporated into baitslib  

## select_snps  
Now vcf2baits  

## selectsnps  
Now vcf2baits  

## tile_probes
Now tilebaits  
