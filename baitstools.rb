#!/usr/bin/env ruby
#-----------------------------------------------------------------------------------------------
# baitstools
BAITSTOOLSVER = "1.7.1"
# Michael G. Campana, 2017-2021
# Smithsonian Conservation Biology Institute
#-----------------------------------------------------------------------------------------------

require 'optparse'
require 'ostruct'
require_relative 'vcf2baits'
require_relative 'tilebaits'
require_relative 'bed2baits'
require_relative 'annot2baits'
require_relative 'aln2baits'
require_relative 'stacks2baits'
require_relative 'checkbaits'
require_relative 'pyrad2baits'
require_relative 'blast2baits'
require_relative 'baitslib'
#-----------------------------------------------------------------------------------------------
class Parser
	def self.parse(options)
		# Set defaults
		args = OpenStruct.new
		algorithms = ["aln2baits", "annot2baits", "bed2baits", "blast2baits", "checkbaits", "pyrad2baits", "stacks2baits", "tilebaits", "vcf2baits"]
		args.algorithm = ARGV[0] # Define subcommands
		if algorithms.include?(args.algorithm) #Set interactive mode or whether to kill
			args.interact = true if ARGV.size == 1
			args.kill = false
		else
			args.kill = true # kill the program if neither interactive nor command line
		end
		ARGV.delete_at(0) # Remove subcommand from input otherwise crashes
		args.threads = 1 # Number of threads requested
		args.used_threads = 1 # Number of threads used
		args.infile = "" # Primary input file
		args.every = false # Flag that baits will be generated from every SNP
		args.totalsnps = 30000 # Maximum requested SNPs
		args.scale = false # Flag to scale SNPs per contig by contig length
		args.scalehash = {} # Hash to get length from contig headers
		args.algorithm == "pyrad2baits" ? args.maxsnps = 1 : args.maxsnps = 2 # Maximum SNPs per contig
		args.algorithm == "pyrad2baits" ? args.distance = 100 : args.distance = 10000 # Minimum distance between SNPs in a contig
		args.no_baits = false # Flag to omit generating baits
		args.refseq = "" # Reference sequence file
		args.lenbef = 60 # Length before SNP in bait
		args.tileoffset = 60 # Offset between tiled baits
		args.tiledepth = 1 # Tiling depth
		args.minind = 1 # Minimum number of individuals to include pyrad locus
		args.strategy = "alignment" # Strategy for pyrad loci
		args.filter = false # Flag to determine whether to run bait filtration
		args.completebait = false # Flag to filter by complete bait length
		args.baitlength = 120 # Bait length
		args.mingc_filter = false # Flag to filter by minimum gc content
		args.mingc = 30.0 # Minimum GC content
		args.maxgc_filter = false # Flag to filter by maximum gc content
		args.maxgc = 50.0 # Maximum GC content
		args.na = 0.9 # Sodium concentration
		args.formamide = 0.0 # Formamide concentration
		args.bait_type = "RNA-DNA" # Hybridization type
		args.mint_filter = false # Flag to filter by minimum melting temperature
		args.mint = 0.0 # Minimum melting temperature
		args.maxt_filter = false # Flag to filter by maximum melting temperature
		args.maxt = 120.0 # Maximum melting temperature
		args.params = false # Flag to output filtration parameters
		args.list_format = "bed" # File format for bed2baits intervals
		args.coords = false # Flag to output BED file of baits (absolute coordinates)
		args.rbed = false # Flag to output BED file of baits (relative coordinates)
		args.percid = 0.0 # Minimum percent identity to include blast hit
		args.blastlen = 1 # Minimum sequence length in nucleotides to include blast hit
		args.evalue_filter = false # Filter by maximum evalue
		args.evalue = 0.1 # Maximum evalue to include blast hit
		args.gaps = "include" # Gapped bait handling
		args.no_Ns = false # Flag to omit bait sequences with Ns
		args.collapse_ambiguities = false # Flag to collapse ambiguities to a single nucleotide
		args.haplodef = "haplotype" # Haplotype definition for aln2baits
		args.uncollapsed_ref = false # Flag to keep ambiguities in pyrad reference sequence
		args.sort = false # Flag to sort stack2baits SNPs by between/within population variation
		args.hwe = false # Flag to sort stacks2baits SNPs by Hardy-Weinberg Equilibrium
		args.alpha = 0.05 # HWE test alpha value
		args.meanqual = 20.0 # Minimum mean base quality
		args.meanqual_filter = false # Flag to filter by minimum mean base quality
		args.minqual = 10 # Minimum base quality
		args.minqual_filter = false # Flag to filter by minimum base quality
		args.fasta_score = 0 # Asssumed base quality score for FASTA sequences
		args.maxmask_filter = false # Flag to filter by maximum masked sequence
		args.maxmask = 25.0 # Maximum percent masked sequence
		args.maxhomopoly_filter = false # Flag to filter by maximum homopolymer length
		args.maxhomopoly = 4 # Maximum homopolymer length
		args.lc_filter = false # Flag to filter by minimum sequence complexity
		args.lc = 0.9 # Minimum sequence complexity
		args.no_lc = false # Disable linguistic complexity calculation in parameters-only calculation
		args.features = [] # Array holding desired features
		args.pad = 0 # BP to pad ends of extracted regions
		args.shuffle = false # Flag whether baits are shuffled to account for end of contigs
		args.ncbi = false # Flag whether FASTA/FASTQ headers include NCBI-style descriptors
		args.rna = false # Flag whether baits are output as RNA
		args.alt_alleles = false # Flag to apply alternate alleles
		args.rc = false # Flag to generate reverse-complement bait
		args.phred64 = false # Flag to use phred64 quality encoding
		args.varqual_filter = false # Flag to determine whether to filter vcf variants by QUAL scores
		args.varqual = 30 # Minimum vcf variant QUAL score
		args.taxafile = nil # File holding taxa IDs
		args.taxa = {} # Taxa hash
		args.taxacount = [0,0,0] # Array of values for taxa balancing
		args.popcategories = nil # Maximum numbers of population-specific variants
		args.previousbaits = nil # Previously generated baits BED file
		args.outdir = File.expand_path("./") # Output directory
		args.outprefix = "out" # Output prefix
		args.log = false # Flag to output detailed log
		args.altbaits = nil # Array of alternative bait lengths. nil turns off flag
		args.infix = '' # Infix for filestems with multiple bait lengths.
		args.filestem = args.outdir + args.outprefix # File stem for output
		args.default_files = [] # Default files to be written
		args.rng = srand # Random number seed
		args.gzip = false # Gzip output files
		args.fiveprime = "" # Sequence to addend to 5' end
		args.threeprime = "" # Sequence to addend to 3' end
		args.fillin = "" # Motif to fill-in short baits
		args.fillin_switch = false # Switch to turn off if filling in
		args.noaddenda = false # Exclude 5' and 3' addended sequences from QC calculations
		args.inbed = nil # Input BED for checkbaits for revised filtering
		opt_parser = OptionParser.new do |opts|
			if algorithms.include?(args.algorithm) # Process correct commands
				opts.banner = "Command-line usage: ruby baitstools.rb "+args.algorithm+" [options]"
				opts.separator ""
				opts.separator args.algorithm+"-specific options:"
				opts.separator ""
				case args.algorithm
				when "vcf2baits" # vcf2baits options
					opts.on("-i","--input [FILE]", String, "Input VCF file") do |vcf|
						args.infile = vcf
					end
					opts.on("--taxafile [FILE]", String, "Balance variants by taxa specified in TSV file") do |taxafile|
						args.taxafile = taxafile
					end
					opts.on("--taxacount [VALUES]", String, "Comma-separated list of values for taxa balancing (Order: AllPopulations,BetweenPopulations,WithinPopulations)") do |taxacount|
						args.taxacount = taxacount.split(",").map! { |value| value.to_i } if taxacount != nil
					end
					opts.on("--popcategories [VALUES]", String, "Comma-separated list of maximum number of population-specific variants in order of appearance in taxa TSV file") do |popcat|
						args.popcategories = popcat.split(",").map! { |value| value.to_i } if popcat != nil
					end
					opts.on("--previousbaits [FILE]", String, "Complement previously generated baits in BED file") do |previousbaits|
						args.previousbaits = previousbaits
					end
					opts.on("-V", "--varqual [VALUE]", Integer, "Minimum variant QUAL score (Default = 30)") do |varf|
						args.varqual = varf if varf != nil
						args.varqual_filter = true
					end
					opts.on("-t", "--totalvars [VALUE]", Integer, "Total requested variants (Default = 30,000)") do |tsnps|
						args.totalsnps = tsnps if tsnps != nil
					end
				when "stacks2baits" # stacks2baits only options
					opts.on("-i","--input [FILE]", String, "Input Stacks summary TSV file") do |fa|
						args.infile = fa
					end
					opts.on("-S", "--sort", "Sort variants according to variation between/within populations") do
						args.sort = true
					end
					opts.on("-H", "--hwe", "Sort variants within populations according to Hardy-Weinberg Equilibrium (Implies -S)") do
						args.sort = true
						args.hwe = true
					end
					opts.on("-A", "--alpha [VALUE]", Float, "Alpha value for HWE test (Default = 0.05)") do |alph|
						args.alpha = alph if alph != nil
					end
					opts.on("-t", "--totalvars [VALUE]", Integer, "Total requested variants per category (Default = 30,000)") do |tsnps|
						args.totalsnps = tsnps if tsnps != nil
					end
				when "pyrad2baits" # pyrad2baits only options
					opts.on("-i","--input [FILE]", String, "Input LOCI file") do |fa|
						args.infile = fa
					end
					opts.on("-L", "--length [VALUE]", Integer, "Requested bait length (Default = 120)") do |prblength|
						args.baitlength = prblength if prblength != nil
					end
					opts.on("-O", "--offset [VALUE]", Integer, "Base pair offset between tiled baits (Default = 60)") do |toff|
						args.tileoffset = toff if toff != nil
					end
					opts.on("-I", "--minind [VALUE]", Integer, "Minimum number of individuals to include locus (Default = 1)") do |mind|
						args.minind = mind if mind != nil
					end
					opts.on("-W", "--strategy [VALUE]", String, "Strategy to generate baits from loci (alignment, SNPs, or informative) (Default = alignment)") do |strat|
						args.strategy = strat.downcase if strat != nil
					end
					opts.on("-H","--haplo [VALUE]", String, "If using alignment strategy, window haplotype definition (haplotype or variant) (Default = haplotype)") do |fa|
						args.haplodef = fa.downcase if fa != nil
					end
					opts.on("--uncollapsedref","Keep ambiguities in pyrad2baits reference sequence") do
						args.uncollapsed_ref = true
					end
					opts.on("-a", "--alt", "If using SNPs or informative strategies, generate baits for alternate alleles") do
 						args.alt_alleles = true
 					end
 					opts.on("-t", "--totalvars [VALUE]", Integer, "If using SNPs or informative strategies, total requested variants (Default = 30,000)") do |tsnps|
						args.totalsnps = tsnps if tsnps != nil
					end
					opts.on("-m","--maxsnps [VALUE]", Integer, "If using SNPs or informative strategies, maximum number of SNPs per locus (Default = 1)") do |msnps|
						args.maxsnps = msnps if msnps != nil
					end
					opts.on("-d","--distance [VALUE]", Integer, "If using SNPs or informative strategies, minimum distance between variants within a locus (Default = 100)") do |dist|
						args.distance = dist if dist != nil
					end
 					opts.on("-k", "--depth [VALUE]", Integer, "If using SNPs or informative strategies, requested baits per variant (Default = 1)") do |tdep|
						args.tiledepth = tdep if tdep != nil
					end
				else
					case args.algorithm
					when "bed2baits"
						inputstr = "Input BED or interval list file"
					when "blast2baits"
						inputstr = "Input BLAST hit table"
					when "annot2baits"
						inputstr = "Input GFF/GTF table"
					else
						inputstr = "Input FASTA/FASTQ file"
					end
					opts.on("-i","--input [FILE]", String, inputstr) do |fa|
						args.infile = fa
					end
					if args.algorithm == "bed2baits" or args.algorithm == "annot2baits" or args.algorithm == "blast2baits"
						if args.algorithm == 'bed2baits'
							opts.on("--list [VALUE]", String, "Interval list file format (BED, GATK, Picard) (Default = BED)") do |strat|
								args.list_format = strat.downcase if strat != nil
							end
						end
						opts.on("-r","--refseq [FILE]", String, "Reference FASTA/FASTQ sequence file") do |ref|
							args.refseq = ref
						end
						opts.on("-P", "--pad [VALUE]", Integer, "Length in bp to pad beginning and ending of extracted regions (Default = 0)") do |pad|
							args.pad = pad if pad != nil
						end
					end
					if args.algorithm == 'checkbaits'
						opts.on("--inbed [FILE]", String, "Optional BED file to be filtered with baits") do |inbed|
							args.inbed = inbed
						end
					end
					opts.on("-L", "--length [VALUE]", Integer, "Requested bait length (Default = 120)") do |prblength|
						args.baitlength = prblength if prblength != nil
					end
					unless args.algorithm == "checkbaits"
						opts.on("-O", "--offset [VALUE]", Integer, "Base pair offset between tiled baits (Default = 60)") do |toff|
							args.tileoffset = toff if toff != nil
						end
					end
					if args.algorithm == "annot2baits"
						opts.on("-U", "--features [FEATURE]", String, "Comma-separated list of feature types") do |feat|
							args.features = feat.upcase.split(",")
						end
					elsif args.algorithm == "aln2baits"
						opts.on("-H","--haplo [VALUE]", String, "Window haplotype definition (haplotype or variant) (Default = haplotype)") do |fa|
							args.haplodef = fa if fa.downcase != nil
						end
					end
					if args.algorithm == "blast2baits"
						opts.on("--percid [VALUE]", Float, "Minimum percent identity to include BLAST hit (Default = 0.0)") do |percid|
							args.percid = percid if percid != nil
						end
						opts.on("--blastlen [VALUE]", Integer, "Minimum length of BLAST hit in nucleotides (Default = 1)") do |blastlen|
							args.blastlen = blastlen if blastlen != nil
						end
						opts.on("--evalue [VALUE]", Float, "Maximum E-value to include BLAST hit (Default = 0.1)") do |evalue|
							args.evalue = evalue if evalue != nil
							args.evalue_filter = true
						end
					end
				end
				if args.algorithm == "vcf2baits" or args.algorithm == "stacks2baits" #vcf2baits, stacks2baits shared options
 					opts.on("-j", "--scale", "Scale maximum number of variants per contig by contig length (Overrides -m)") do
						args.scale = true
					end
					opts.on("-m","--maxsnps [VALUE]", Integer, "Maximum number of variants per contig (Default = 2)") do |msnps|
						args.maxsnps = msnps if msnps != nil
					end
					opts.on("-d","--distance [VALUE]", Integer, "Minimum distance between variants within a contig (Default = 10,000)") do |dist|
						args.distance = dist if dist != nil
					end
					opts.on("-p","--nobaits", "Do not output bait sequences") do
						args.no_baits = true
					end
					 opts.on("-e","--every", "Output bait sequences for every variant in the input file") do
						args.every = true
					end
					 opts.on("-a", "--alt", "Generate baits for alternate alleles (Overrides -p)") do
 						args.alt_alleles = true
 					end
					opts.on("-r","--refseq [FILE]", String, "If outputting baits, the reference FASTA/FASTQ sequence file") do |ref|
						args.refseq= ref
					end
					opts.on("-L", "--length [VALUE]", Integer, "Requested total bait length (Default = 120)") do |prblength|
						args.baitlength = prblength if prblength != nil
					end
					opts.on("-b", "--lenbef [VALUE]", Integer, "Base pairs before the variant to start generating baits (Default = 60)") do |b4|
						args.lenbef = b4 if b4 != nil
					end
					opts.on("-O", "--offset [VALUE]", Integer, "Base pair offset between tiled baits (Default = 60)") do |toff|
						args.tileoffset = toff if toff != nil
					end
					opts.on("-k", "--depth [VALUE]", Integer, "Requested tiles per variant (Default = 1)") do |tdep|
						args.tiledepth = tdep if tdep != nil
					end
				end
				opts.separator "" #bait filtering options
				opts.separator "Bait filtration options:"
				opts.separator ""
				opts.on("-w", "--params", "Output bait statistics table") do
					args.params = true
				end
				opts.on("--disable-lc", "Disable slow linguistic complexity calculations for parameters table") do
					args.no_lc = true
				end
				opts.on("-c","--complete", "Require baits be full length") do
					args.completebait = true
				end
				opts.on("-N", "--noNs", "Exclude bait sequences with Ns") do
					args.no_Ns = true
				end
				opts.on("--noaddenda", "Exclude 5' and 3' addended sequences from bait parameter calculations") do
					args.noaddenda = true
				end
				opts.on("-n","--mingc [VALUE]", Float, "Minimum GC content in % (Default = 30.0)") do |mgc|
					args.mingc = mgc if mgc != nil
					args.mingc_filter = true
				end
				opts.on("-x","--maxgc [VALUE]", Float, "Maximum GC content in % (Default = 50.0)") do |xgc|
					args.maxgc = xgc if xgc != nil
					args.maxgc_filter = true
				end
				opts.on("-q","--mint [VALUE]", Float, "Minimum melting temperature in °C (Default = 0.0)") do |mt|
					args.mint = mt if mt != nil
					args.mint_filter = true # Force bait filtration if called
				end
				opts.on("-z","--maxt [VALUE]", Float, "Maximum melting temperature in °C (Default = 120.0)") do |xt|
					args.maxt = xt if xt != nil
					args.maxt_filter = true # Force bait filtration if called
				end
				opts.on("-T", "--type [VALUE]", String, "Melting temperature for DNA-DNA, RNA-RNA, or RNA-DNA hybridization (Default = RNA-DNA)") do |hyb|
					args.bait_type = hyb.upcase if hyb != nil
				end
				opts.on("-s","--na [VALUE]", Float, "Melting temperature sodium concentration in M (Default = 0.9)") do |na|
					args.na = na if na != nil
				end
				opts.on("-f", "--formamide [VALUE]", Float, "Melting temperature formamide concentration in M (Default = 0.0)") do |form|
					args.formamide = form if form != nil
				end
				opts.on("-K", "--maxmask [VALUE]", Float, "Maximum % sequence consisting of masked elements (Default = 25.0)") do |maxrep|
					args.maxmask = maxrep if maxrep != nil
					args.maxmask_filter = true # Force bait filtration if called
				end
				opts.on("-J", "--maxhomopoly [VALUE]", Integer, "Maximum homopolymer length (Default = 4)") do |maxhomopoly|
					args.maxhomopoly = maxhomopoly if maxhomopoly != nil
					args.maxhomopoly_filter = true
				end
				opts.on("-y", "--minlc [VALUE]", Float, "Minimum sequence linguistic complexity  (Default = 0.9)") do |minlc|
					args.lc = minlc if minlc != nil
					args.lc_filter = true
				end
				if args.algorithm != "aln2baits"
					opts.on("-Q", "--meanqual [VALUE]", Float, "Minimum mean base quality (Default = 20.0)") do |meanqual|
						args.meanqual = meanqual if meanqual != nil
						args.meanqual_filter = true
					end
					opts.on("-M", "--minqual [VALUE]", Integer, "Minimum base quality (Default = 10)") do |minqual|
						args.minqual = minqual if minqual != nil
						args.minqual_filter = true
					end
					opts.on("-F", "--fastascore [VALUE]", Integer, "Assumed FASTA-base quality (Default = 0)") do |fasqual|
						args.fasta_score = fasqual if fasqual != nil
					end
				end
				opts.separator ""
				opts.separator "Common options:"
				opts.separator ""
				opts.on("-o", "--outprefix [VALUE]", String, "Output file prefix (Default = out)") do |out|
					args.outprefix = out if out != nil
				end
				opts.on("-Z", "--outdir [VALUE]", String, "Output directory path (Default is ./)") do |outdir|
					args.outdir = File.expand_path(outdir)  if outdir != nil
				end
				opts.on("-l", "--log", "Output detailed log") do |log|
					args.log = true
				end
				unless args.algorithm == "checkbaits" # checkbaits does not get coordinate information to output
					opts.on("-B", "--bed", "Output BED file for the baits") do
						args.coords = true
					end
					if args.algorithm == "aln2baits" or args.algorithm == "annot2baits" or args.algorithm == "bed2baits" or args.algorithm == "tilebaits" or args.algorithm == "blast2baits"
						opts.on("-E", "--rbed", "Output BED file for the baits relative to extracted sequences") do
							args.rbed = true
						end
					end
					opts.on("--shuffle", "Shuffle baits to compensate for extending beyond contig ends") do
						args.shuffle = true
					end
				end
				unless args.algorithm == "pyrad2baits"
					opts.on("-D", "--ncbi", "FASTA/FASTQ file headers have NCBI-style descriptions") do
						args.ncbi = true
					end
				end
				opts.on("--phred64", "Quality scores are in Phred64") do
					args.phred64 = true
				end
				opts.on("-C", "--collapse", "Collapse ambiguities to a single nucleotide") do
					args.collapse_ambiguities = true
				end
				opts.on("-Y","--rna", "Output baits as RNA rather than DNA") do
					args.rna = true
				end
				opts.on("-R", "--rc", "Output reverse-complemented baits") do
					args.rc = true
				end
				opts.on("-G", "--gaps [VALUE]", String, "Strategy to handle sequence gaps (-) (include, exclude, or extend) (Default = include)") do |gap|
					args.gaps = gap.downcase if gap != nil
				end
				opts.on("--altbaits [VALUES]", String, "Comma-separated list of alternative bait lengths") do |altbaits|
					args.altbaits = altbaits.split(",").map! { |value| value.to_i } if altbaits != nil
				end
				opts.on("-5", "--5prime [VALUE]", String, "Sequence to addend to 5' end of baits") do |fivepr|
					args.fiveprime = fivepr if fivepr != nil
				end
				opts.on("-3", "--3prime [VALUE]", String, "Sequence to addend to 3' end of baits") do |threepr|
					args.threeprime = threepr if threepr != nil
				end
				opts.on("--fillin [VALUE]", String, "Fill in baits shorter than requested length with specified sequence repeat motif") do |fillin|
					if fillin != nil
						args.fillin = fillin
						args.fillin_switch = true
					end
				end
				opts.on("-X", "--threads [VALUE]", Integer, "Number of threads (Default = 1)") do |thr|
					args.threads = thr if thr != nil
				end
				opts.on("--rng [VALUE]", Integer, "Random number seed (Default uses system entropy)") do |rng|
					args.rng = rng if rng != nil
				end
				opts.on("--gzip", "Gzip output files") do
					args.gzip = true
				end
				opts.on_tail("-h","--help", "Show help") do
					print "Welcome to baitstools " + args.algorithm + '.' +"\n\n"
					print "To use the interactive interface, enter <ruby baitstools.rb " + args.algorithm + "> without command-line options.\n"
					puts opts
					exit
				end
				opts.on_tail("-v","--version","Show subcommand version") do
					case args.algorithm
					when "vcf2baits"
						print "vcf2baits "+ VCF2BAITSVER + "\n"
					when "tilebaits"
						print "tilebaits "+ TILEBAITSVER + "\n"
					when "bed2baits"
						print "bed2baits "+ BED2BAITSVER + "\n"
					when "annot2baits"
						print "annot2baits "+ ANNOT2BAITSVER + "\n"
					when "aln2baits"
						print "aln2baits " + ALN2BAITSVER + "\n"
					when "stacks2baits"
						print "stacks2baits " + STACKS2BAITSVER + "\n"
					when "checkbaits"
						print "checkbaits " + CHECKBAITSVER + "\n"
					when "pyrad2baits"
						print "pyrad2baits " + PYRAD2BAITSVER + "\n"
					when "blast2baits"
						print "blast2baits " + BLAST2BAITSVER + "\n"
					end
					exit
				end
			else
				opts.banner = "Common options:"
				print "Welcome to baitstools " + BAITSTOOLSVER + "\n\nTo use the interactive interface, enter <ruby baitstools.rb [subcommand]> without command-line options.\nCommand-line usage: ruby baitstools.rb [subcommand] [options]"
				print "\nAdd '-h' or '--help' to subcommands (without other options) to see their relevant options.\n\nAvailable subcommands:\n\n"
				print "    aln2baits\t\t\t\tGenerate weighted baits from a FASTA/FASTQ alignment\n"
				print "    annot2baits\t\t\t\tGenerate tiled baits from a GFF/GTF file and a reference sequence\n"
				print "    bed2baits\t\t\t\tGenerate tiled baits from BED/interval list file and a reference sequence\n"
				print "    blast2baits\t\t\t\tGenerate tiled baits from a BLAST hit table and a reference sequence\n"
				print "    checkbaits\t\t\t\tFilter a FASTA/FASTQ of candidate baits by quality\n"
				print "    pyrad2baits\t\t\t\tSelect variants and generate baits from a PyRAD/ipyrad LOCI file\n"
				print "    stacks2baits\t\t\tSelect variants and generate baits from a Stacks summary TSV file\n"
				print "    tilebaits\t\t\t\tGenerate tiled baits from FASTA/FASTQ sequences\n"
				print "    vcf2baits\t\t\t\tSelect variants and generate baits from a VCF\n\n"

			end
		end
		opt_parser.parse!(options)
		return args
	end
end
#-----------------------------------------------------------------------------------------------
def ynq # Yes-No question handling
	t = gets.chomp.upcase
	(t == "Y" || t == "YES" ) ? val = true : val = false
	return val
end
#-----------------------------------------------------------------------------------------------
begin
	$options = Parser.parse(ARGV)
	exit if $options.kill
	# Interactive mode block
	if $options.interact
		print "Enter input file.\n"
		$options.infile = gets.chomp
	end
	while !FileTest.exist?($options.infile)
		print "Input file not found. Please re-enter.\n"
		$options.infile = gets.chomp
	end
	if $options.algorithm == 'checkbaits'
		if $options.interact
			print "Filter corresponding BED file? (y/n)\n"
			if ynq
				print "Enter BED file.\n"
				$options.inbed = gets.chomp
			end
		end
		unless $options.inbed.nil?
			while !FileTest.exist?($options.inbed)
				print "BED file not found. Please re-enter.\n"
				$options.inbed = gets.chomp
			end
		end
	end
	if $options.algorithm == "bed2baits"
			if $options.interact
				print "File format (BED, GATK, or Picard)?\n"
				$options.list_format = gets.chomp.downcase
			end
			while $options.list_format != 'bed' and $options.list_format != 'gatk' and $options.list_format != 'picard'
				print "Please choose a recognized file format (BED, GATK, or Picard)?\n"
				$options.list_format = gets.chomp.downcase
			end
		end
	if $options.interact
		print "Enter output file prefix.\n"
		$options.outprefix = gets.chomp
		print "Enter output file directory.\n"
		$options.outdir = File.expand_path(gets.chomp)
		print "Output detailed log? (y/n)\n"
		$options.log = ynq
	end
	Dir.mkdir ($options.outdir) if !FileTest.directory?($options.outdir)
	if $options.interact
		print "Number of threads?\n"
		$options.threads = gets.chomp.to_i
	end
	while $options.threads < 1
		print "Number of threads must be at least 1. Re-enter\n"
		$options.threads = gets.chomp.to_i
	end
	if $options.algorithm ==  "vcf2baits" or $options.algorithm == "stacks2baits" #algorithm-specific options
		if $options.interact
			print "Output baits for all variants in input file? (y/n)\n"
			$options.every = ynq
			if $options.algorithm == "stacks2baits"
				print "Sort variants according to variation between/within populations? (y/n)\n"
				if ynq
					$options.sort = true
					print "Sort variants within populations according to Hardy-Weinberg Equilibrium? (y/n)\n"
					if ynq
						$options.hwe = true
						print "Alpha value for HWE?\n"
						$options.alpha = gets.chomp.to_f
					end
				end
			end
		end
		if $options.algorithm == "stacks2baits"
			while $options.alpha < 0.0 or $options.alpha > 1.0
				print "Alpha value must be between 0.0 and 1.0. Re-enter.\n"
				$options.alpha = gets.chomp.to_f
			end
		end
		if $options.algorithm == "vcf2baits" and !$options.every
			if $options.interact
				print "Complement previously generated baits?\n"
				if ynq
					print "Enter baits BED file name.\n"
					$options.previousbaits = gets.chomp
					while !FileTest.exist?($options.previousbaits)
						print "BED file not found. Please re-enter.\n"
						$options.previousbaits = gets.chomp
					end
				end
				print "Sort variants by variation within, between and among taxa?\n"
				if ynq
					print "Enter taxa file name.\n"
					$options.taxafile = gets.chomp
				end
			end
			unless $options.taxafile.nil?
				while !FileTest.exist?($options.taxafile)
					print "Taxa file not found. Please re-enter.\n"
					$options.taxafile = gets.chomp
				end
				$options.totalsnps = $options.taxacount.sum
				while $options.totalsnps < 1 or $options.taxacount.any? { |x| x < 0 }
					if $options.interact
						print "Enter requested number of variants variable across populations.\n"
						$options.taxacount[0] = gets.chomp.to_i
					end
					while $options.taxacount[0] < 0
						print "The number of variants variable across populations must be at least 0. Re-enter.\n"
						$options.taxacount[0] = gets.chomp.to_i
					end
					if $options.interact
						print "Enter requested number of variants between populations.\n"
						$options.taxacount[1] = gets.chomp.to_i
					end
					while $options.taxacount[1] < 0
						print "The number of variants variable between populations must be at least 0. Re-enter.\n"
						$options.taxacount[1] = gets.chomp.to_i
					end
					if $options.interact
						print "Enter requested number of variants variable within populations.\n"
						$options.taxacount[2] = gets.chomp.to_i
					end
					while $options.taxacount[2] < 0
						print "The number of variants variable within populations must be at least 0. Re-enter.\n"
						$options.taxacount[2] = gets.chomp.to_i
					end
					$options.totalsnps = $options.taxacount[0] + $options.taxacount[1] + $options.taxacount[2] # Ruby 2.0 does not have Array#sum method
					if $options.totalsnps < 1
						print "The total number of variants must be greater than 0. Please re-enter.\n"
						$options.taxacount = [-1,-1,-1] # Force taxacount redo
					end
				end
			end
			if $options.interact and !$options.taxafile.nil?
				print "Control population-specific variants by quantity?\n"
				if ynq
					print "Enter comma-separated list of population-specific variants in order of appearance in taxa TSV file\n"
					$options.popcategories = gets.chomp.split(",").map! { |val| val.to_i }
				end
			end
			unless $options.popcategories.nil? # Check that population values are ok
				if !$options.taxafile.nil? # Only check the values if a taxafile exists
					while $options.popcategories.any? { |x| x < 0 or x > $options.taxacount[2] }
						print "Population-specific variant numbers must be between 0 and the total number of within-population variants. Re-enter.\n"
						$options.popcategories = gets.chomp.split(",").map! { |val| val.to_i }
					end
				else
					$options.popcategories = nil # Reset bad popcategories to nil to prevent crashes
				end
			end
		end
		if $options.interact and !$options.every and $options.taxafile.nil?
			if $options.algorithm == "vcf2baits"
				print "Enter total number of requested variants.\n"
			else $options.algorithm == "stacks2baits"
				print "Enter total number of requested variants per category.\n"
			end
			$options.totalsnps = gets.chomp.to_i
		end
		while $options.totalsnps < 1
			print "The total number of variants must be greater than 0. Please re-enter.\n"
			$options.totalsnps = gets.chomp.to_i
		end
		if $options.interact and !$options.every
			print "Scale maximum number of variants per contig by contig length? (y/n)\n"
			$options.scale = ynq
		end
		if $options.interact and !$options.every and !$options.scale
			print "Enter maximum number of variants per contig.\n"
			$options.maxsnps = gets.chomp.to_i
		end
		while $options.maxsnps < 1 and !$options.scale
			print "The maximum number of variants per contig must be greater than 0. Please re-enter.\n"
			$options.maxsnps = gets.chomp.to_i
		end
		if $options.interact and !$options.every
			print "Enter minimum distance between variants within a contig.\n"
			$options.distance = gets.chomp.to_i
		end
		while $options.distance < 1
			print "The minimum distance between variants must be greater than 0. Please re-enter.\n"
			$options.distance = gets.chomp.to_i
		end
		if $options.interact and !$options.every
			print "Output baits? (y/n)\n"
			ynq ? $options.no_baits = false : $options.no_baits = true
		end
		unless $options.no_baits
			if $options.interact
				print "Enter reference sequence.\n"
				$options.refseq = gets.chomp
			end
			while !FileTest.exist?($options.refseq)
				print "Reference sequence not found. Please re-enter.\n"
				$options.refseq = gets.chomp
			end
			if $options.interact
				print "Generate baits for alternate alleles? (y/n)\n"
				$options.alt_alleles = ynq
				print "Enter bait length.\n"
				$options.baitlength = gets.chomp.to_i
			end
			while $options.baitlength < 1
				print "Baits must be at least 1 bp long. Re-enter.\n"
				$options.baitlength = gets.chomp.to_i
			end
			if $options.interact
				print "Enter number of bases before variant in bait.\n"
				$options.lenbef = gets.chomp.to_i
			end
			while $options.lenbef < 0
				print "The number of bait bases before the variant must be at least 0. Please re-enter.\n"
				$options.lenbef = gets.chomp.to_i
			end
			if $options.interact
				print "Enter tiling bp offset.\n"
				$options.tileoffset = gets.chomp.to_i
			end
			while  $options.tileoffset < 1
				print "Tiling offset cannot be less than 1. Please re-enter.\n"
				$options.tileoffset = gets.chomp.to_i
			end
			if $options.interact
				print "Enter number of baits per SNP.\n"
				$options.tiledepth = gets.chomp.to_i
			end
			while $options.tiledepth < 1
				print "Tiles per SNP cannot be less than 1. Re-enter.\n"
				$options.tiledepth = gets.chomp.to_i
			end
		end
	else
		if $options.algorithm == "bed2baits" or $options.algorithm == "annot2baits" or $options.algorithm == "blast2baits"
			if $options.interact
				print "Enter reference sequence.\n"
				$options.refseq = gets.chomp
			end
			while !FileTest.exist?($options.refseq)
				print "Reference sequence not found. Please re-enter.\n"
				$options.refseq = gets.chomp
			end
			if $options.interact
				print "Length to pad extracted regions.\n"
				$options.pad = gets.chomp.to_i
			end
			while $options.pad < 0
				print "Pad length cannot be less than 0. Please re-enter.\n"
				$options.pad = gets.chomp.to_i
			end
		end
		if $options.interact
			print "Enter bait length.\n"
			$options.baitlength = gets.chomp.to_i
		end
		while $options.baitlength < 1
			print "Baits must be at least 1 bp long. Please re-enter.\n"
			$options.baitlength = gets.chomp.to_i
		end
		if $options.interact && !$options.algorithm == "checkbaits"
			print "Enter tiling bp offset.\n" 
			$options.tileoffset = gets.chomp.to_i
		end
		while $options.tileoffset < 1
			print "Tiling offset cannot be less than 1. Please re-enter.\n"
			$options.tileoffset = gets.chomp.to_i
		end
		if $options.algorithm == "pyrad2baits"
			while $options.tileoffset > $options.baitlength 
				print "Tiling offset cannot be greater than bait length. Please re-enter.\n"
				$options.tileoffset = gets.chomp.to_i
			end
			if $options.interact
				print "Enter minimum number of individuals to retain locus.\n"
				$options.minind = gets.chomp.to_i
			end
			while $options.minind < 1
				print "Minimum number of individuals must be at least 1. Re-enter.\n"
				$options.minind = gets.chomp.to_i
			end
			if $options.interact
				print "Strategy? (alignment, SNPs, or informative)\n"
				$options.strategy = gets.chomp.downcase
			end
			while $options.strategy != 'alignment' and $options.strategy != 'snps' and $options.strategy != 'informative'
				print "Please choose a strategy (alignment, SNPs, or informative)\n"
				$options.strategy = gets.chomp.downcase
			end
			if $options.strategy != "alignment"
				if $options.interact
					print "Keep ambiguities in pyrad2baits reference sequence? (y/n)\n"
					$options.uncollapsed_ref = ynq
					print "Enter total number of requested variants.\n"
					$options.totalsnps = gets.chomp.to_i
				end
				while $options.totalsnps < 1
					print "The total number of variants must be greater than 0. Please re-enter.\n"
					$options.totalsnps = gets.chomp.to_i
				end
				if $options.interact
					print "Enter maximum number of variants per locus.\n"
					$options.maxsnps = gets.chomp.to_i
				end
				while $options.maxsnps < 1
					print "The maximum number of variants per locus must be greater than 0. Please re-enter.\n"
					$options.maxsnps = gets.chomp.to_i
				end
				if $options.interact
					print "Enter minimum distance between variants within a contig.\n"
					$options.distance = gets.chomp.to_i
				end
				while $options.distance < 1
					print "The minimum distance between variants must be greater than 0. Please re-enter.\n"
					$options.distance = gets.chomp.to_i
				end
				if $options.interact
					print "Enter number of bases before variant in bait.\n"
					$options.lenbef = gets.chomp.to_i
				end
				while $options.lenbef < 0
					print "The number of bait bases before the variant must be at least 0. Please re-enter.\n"
					$options.lenbef = gets.chomp.to_i
				end
				if $options.interact
					print "Enter number of baits per SNP.\n"
					$options.tiledepth = gets.chomp.to_i
				end
				while $options.tiledepth > $options.baitlength/$options.tileoffset or $options.tiledepth < 1
					print "Tiling depth cannot be less than 1 or greater than bait length/tiling offset ratio. Re-enter.\n"
					$options.tiledepth = gets.chomp.to_i
				end
			end
		elsif $options.algorithm == "annot2baits" and $options.interact
			print "Enter comma-separated list of features.\n"
			$options.features = gets.chomp.upcase.split(",")
		end
		if $options.algorithm == "aln2baits" or ($options.algorithm == "pyrad2baits" && $options.strategy == "alignment")
			if $options.interact
				print "Haplotype definition? (haplotype or variant)\n"
				$options.haplodef = gets.chomp.downcase
			end
			while $options.haplodef != 'haplotype' and $options.haplodef != 'variant'
				print "Please choose a haplotype definition (haplotype or variant)\n"
				$options.haplodef = gets.chomp.downcase
			end
		end
		if $options.algorithm == "blast2baits"
			if $options.interact
				print "Enter minimum percent identity to include BLAST hit.\n"
				$options.percid = gets.chomp.to_f
			end
			while $options.percid < 0.0 or $options.percid > 100.0
				print "Percent identity must be between 0.0 and 100.0. Re-enter.\n"
				$options.percid = gets.chomp.to_f
			end
			if $options.interact
				print "Enter minimum length of BLAST hit in nucleotides.\n"
				$options.blastlen = gets.chomp.to_i
			end
			while $options.blastlen < 0
				print "Minimum length must be greater than 0. Re-enter.\n"
				$options.blastlen = gets.chomp.to_i
			end
			if $options.interact
				print "Filter BLAST hits by E-value?\n"
				$options.evalue_filter = ynq
				if $options.evalue_filter
					print "Enter maximum E-value to include BLAST hit.\n"
					$options.evalue = gets.chomp.to_f
				end
			end
			while $options.evalue_filter && $options.evalue < 0
				print "Maximum E-value must be at least 0.0 Re-enter.\n"
				$options.evalue = gets.chomp.to_f
			end
		end
	end
	if $options.interact and !$options.no_baits
		if $options.algorithm != "pyrad2baits"
			print "Do FASTA/FASTQ sequence headers include NCBI-style descriptions? (y/n)\n"
			$options.ncbi = ynq
		end
		if $options.algorithm != "checkbaits"
			print "Output BED file(s) for candidate baits? (y/n)\n"
			$options.coords = ynq
			if $options.algorithm == "aln2baits" or $options.algorithm == "annot2baits" or $options.algorithm == "bed2baits" or $options.algorithm == "tilebaits" or $options.algorithm == "blast2baits"
				print "Output BED file(s) for the baits relative to extracted sequences? (y/n)\n"
				$options.rbed = ynq
			end
			print "Shuffle baits to compensate for extending beyond contig ends? (y/n)\n"
			$options.shuffle = ynq
		end
		print "Output bait statistics table? (y/n)\n"
		$options.params = ynq
		if $options.params
			print "Disable linguistic complexity calculations to improve program speed? (y/n)\n"
			$options.no_lc = ynq
		end
		print "Require complete length bait? (y/n)\n"
		$options.completebait = ynq
		$options.gaps = "" # Turn off default to allow value entry entry
	end
	while $options.gaps != 'include' and $options.gaps != 'exclude' and $options.gaps != 'extend'
		print "Strategy to handle sequence gaps (-) ('include', 'exclude', or 'extend')\n"
		$options.gaps = gets.chomp.downcase
	end
	if $options.interact and !$options.no_baits
		print "Exclude baits with Ns? (y/n)\n"
		$options.no_Ns = ynq
		print "Collapse ambiguity codes to a single nucleotide? (y/n)\n"
		$options.collapse_ambiguities = ynq
		print "Output baits as RNA rather than DNA? (y/n)\n"
		$options.rna = ynq
		print "Generate reverse complement baits? (y/n)\n"
		$options.rc = ynq
		print "Generate alternate length baits? (y/n)\n"
		if ynq
			print "Enter alternate lengths requested separated by commas.\n"
			$options.altbaits = gets.chomp.split(",").map! { |val| val.to_i }
		end
	end
	unless $options.altbaits.nil?
		while $options.altbaits.any? { |x| x < 1 or ($options.algorithm == 'checkbaits' && x > $options.baitlength)}
			if $options.altbaits.any? { |x| x < 1 }
				print "Baits must at least 1 bp. Re-enter.\n"
			else
				print "Alternate baits cannot be longer than previously generated baits. Re-enter.\n"
			end
			$options.altbaits = gets.chomp.split(",").map! { |val| val.to_i }
		end
	end
	if $options.interact and !$options.no_baits
		print "Addend a sequence to 5' end of baits? (y/n)\n"
		if ynq
			print "Enter 5' sequence.\n"
			$options.fiveprime = gets.chomp
		end
		print "Addend a sequence to 3' end of baits? (y/n)\n"
		if ynq
			print "Enter 3' sequence.\n"
			$options.threeprime = gets.chomp
		end
		unless $options.fiveprime == "" && $options.threeprime == ""
			print "Exclude sequence addenda from bait filtration parameter calculations? (y/n)\n"
			$options.noaddenda = ynq
		end
		print "Fill in incomplete baits with a specified sequence repeat motif? (y/n)\n"
		if ynq
			print "Enter sequence motif.\n"
			$options.fillin = gets.chomp
			$options.fillin_switch = true unless $options.fillin == ""
		end
		print "Filter by minimum GC content? (y/n)\n"
		if ynq
			$options.mingc_filter = true
			print "Enter minimum GC content.\n"
			$options.mingc = gets.chomp.to_f
			while $options.mingc < 0.0 or $options.mingc > 100.0
				print "Minimum GC content must be between 0 and 100. Please re-enter.\n"
				$options.mingc = gets.chomp.to_f
			end
		end
		print "Filter by maximum GC content? (y/n)\n"
		if ynq
			$options.maxgc_filter = true
			print "Enter maximum GC content.\n"
			$options.maxgc = gets.chomp.to_f
			while $options.maxgc < $options.mingc or $options.maxgc > 100.0
				print "Maximum GC content must be between minimum GC content and 100. Re-enter.\n"
				$options.maxgc = gets.chomp.to_f
			end
		end
		print "Filter by minimum melting temperature? (y/n)\n"
		if ynq
			$options.mint_filter = true
			print "Enter minimum melting temperature.\n"
			$options.mint = gets.chomp.to_f
		end
		print "Filter by maximum melting temperature? (y/n)\n"
		if ynq
			$options.maxt_filter = true
			print "Enter maximum melting temperature.\n"
			$options.maxt = gets.chomp.to_f
		end
		if ($options.mint_filter && $options.maxt_filter)
			while $options.maxt < $options.mint
				print "Maximum melting temperature cannot be less than minimum melting temperature. Re-enter maximum temperature.\n"
				$options.maxt = gets.chomp.to_f
			end
		end
	end
	if ($options.mint_filter or $options.maxt_filter)
		if $options.interact
			$options.bait_type = "" # Turn off bait_type default so that interactive mode can choose
			print "Enter sodium concentration.\n"
			$options.na = gets.chomp.to_f
		end
		while $options.na < 0.0
			print "Sodium concentrations cannot be negative. Re-enter.\n"
			$options.na = gets.chomp.to_f
		end
		if $options.interact
			print "Enter formamide concentration.\n"
			$options.formamide = gets.chomp.to_f
		end
		while $options.formamide < 0.0
			print "Formamide concentrations cannot be negative. Re-enter.\n"
			$options.formamide = gets.chomp.to_f
		end
		while $options.bait_type != 'RNA-RNA' and $options.bait_type != 'DNA-DNA' and $options.bait_type != 'RNA-DNA'
			print "Please choose a hybridization type ('DNA-DNA', 'RNA-RNA', or 'RNA-DNA')\n"
			$options.bait_type = gets.chomp.upcase
		end
	end
	if $options.interact and !$options.no_baits
		print "Filter by maximum percent masked sequence? (y/n)\n"
		if ynq
			$options.maxmask_filter = true
			print "Enter maximum percent masked sequence\n"
			$options.maxmask = gets.chomp.to_f
		end
	end
	while $options.maxmask < 0.0 or $options.maxmask > 100.0
		print "Percent masked sequence must be between 0.0 and 100.0. Re-enter.\n"
		$options.maxmask = gets.chomp.to_f
	end
	if $options.interact and !$options.no_baits
		print "Filter by maximum homopolymer length? (y/n)\n"
		if ynq
			$options.maxhomopolymer_filter = true
			print "Enter maximum homopolymer length\n"
			$options.maxhomopoly = gets.chomp.to_i
		end
	end
	while $options.maxhomopoly < 1
		print "Maximum homopolymer length must be greater than 0. Re-enter.\n"
		$options.maxhomopoly = gets.chomp.to_i
	end
	if $options.interact and !$options.no_baits and !$options.no_lc
		print "Filter by minimum sequence linguistic complexity? (y/n)\n"
		if ynq
			$options.lc_filter = true
			print "Enter minimum linguistic complexity\n"
			$options.lc = gets.chomp.to_f
		end
	end
	while $options.lc < 0.0 or $options.lc > 1.0
		print "Linguistic complexity must be between 0.0 and 1.0. Re-enter.\n"
		$options.lc = gets.chomp.to_f
	end			
	if $options.interact and $options.algorithm != "aln2baits" and !$options.no_baits
		print "Filter by minimum mean base quality? (y/n)\n"
		if ynq
			$options.meanqual_filter = true
			print "Enter minimum mean base quality.\n"
			$options.meanqual = gets.chomp.to_f
		end
	end
	while $options.meanqual < 0.0 or $options.meanqual > 93.0
		print "Mean quality must be between 0.0 and 93.0. Re-enter.\n"
		$options.meanqual = gets.chomp.to_f
	end
	if $options.interact and $options.algorithm != "aln2baits" and !$options.no_baits
		print "Filter by minimum base quality? (y/n)\n"
		if ynq
			$options.minqual_filter = true
			print "Enter minimum mean base quality.\n"
			$options.minqual = gets.chomp.to_i
		end
	end
	while $options.minqual < 0 or $options.minqual > 93
		print "Minimum quality must be between 0 and 93. Re-enter.\n"
		$options.minqual = gets.chomp.to_i
	end
	if $options.interact and $options.algorithm != "aln2baits" and !$options.no_baits
		if ($options.meanqual_filter or $options.minqual_filter)
			print "Enter assumed FASTA quality score.\n"
			$options.fasta_score = gets.chomp.to_i
		end
	end
	while $options.fasta_score < 0 or $options.fasta_score > 93
		print "Assumed quality must be between 0 and 93. Re-enter.\n"
		$options.fasta_score = gets.chomp.to_i
	end	
	if $options.interact
		print "Are FASTQ qualities in Phred64? (y/n)\n"
		$options.phred64 = ynq
		print "Set random number seed (otherwise uses system entropy)? (y/n)\n"
		if ynq
			print "Enter random number seed.\n"
			$options.rng = gets.chomp.to_i
		end
		print "Gzip output files? (y/n)\n"
		$options.gzip = ynq
	end	
	$options.no_baits = false if ($options.every or $options.alt_alleles) # Override -p when needed
	$options.filter = true if ($options.completebait or $options.params or $options.algorithm == "checkbaits" or $options.lc_filter or $options.mingc_filter or $options.maxgc_filter or $options.mint_filter or $options.maxt_filter or $options.maxmask_filter or $options.maxhomopoly_filter or $options.meanqual_filter or $options.minqual_filter or $options.gaps == "exclude" or $options.no_Ns or $options.collapse_ambiguities) # Force filtration as necessary
	cmdline = get_command_line
	srand($options.rng) # Set random number seed
	print "** Starting program with the following options: **\n"
	print "** Basic command: " + cmdline[0] + " **\n"
	print "** Filtration options:" + cmdline[1] + " **\n" # filtered line always starts with a space if present
	versline = "Using BaitsTools v. " + BAITSTOOLSVER + ", " + $options.algorithm + " v. " + eval($options.algorithm.upcase + "VER") + ", and baitslib v. #{BAITSLIBVER}"
	print "** " + versline + " **\n"
	setup_output
	build_fq_hash # Build FASTQ quality translation hash
	build_rc_hash # Build reverse complementation hash
	build_ambig_hash # Build ambiguity hash
	write_file(".log.txt", "Parsed commands: " + cmdline[0] + cmdline[1] + "\n" + versline + "\n") if $options.log
	case $options.algorithm
	when "aln2baits"
		aln2baits($options.infile)
	when "annot2baits"
		annot2baits
	when "bed2baits"
		bed2baits
	when "blast2baits"
		blast2baits
	when "checkbaits"
		checkbaits
	when "pyrad2baits"
		pyrad2baits
	when "stacks2baits"
		stacks2baits
	when "tilebaits"
		tilebaits($options.infile)
	when "vcf2baits"
		vcf2baits
	end
	if $options.gzip && $options.log # Gzip log after program complete
		system("gzip #{resolve_unix_path($options.filestem + '.log.txt')}")
	end
	print "** Program complete **\n"
end
