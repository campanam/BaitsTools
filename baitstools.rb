#!/usr/bin/ruby
#-----------------------------------------------------------------------------------------------
# baitstools 0.5
# Michael G. Campana, 2017
# Smithsonian Conservation Biology Institute
#-----------------------------------------------------------------------------------------------

require 'optparse'
require 'ostruct'
require './baitslib.rb'
require './vcf2baits.rb'
require './tilebaits.rb'
require './coords2baits.rb'
require './annot2baits.rb'
require './aln2baits.rb'
require './stacks2baits.rb'
require './checkbaits.rb'
#-----------------------------------------------------------------------------------------------
class Parser
	def self.parse(options)
		# Set defaults
		args = OpenStruct.new
		algorithms = ["aln2baits", "annot2baits", "checkbaits", "coords2baits", "stacks2baits", "tilebaits", "vcf2baits"]
		args.algorithm = ARGV[0] # Define subcommands
		if algorithms.include?(args.algorithm) #Set interactive mode or whether to kill
			args.interact = true if ARGV.size == 1
			args.kill = false
		else
			args.kill = true # kill the program if neither interactive nor command line
		end
		ARGV.delete_at(0) # Remove subcommand from input otherwise crashes
		args.infile = "" # Primary input file
		args.every = false # Flag that baits will be generated from every SNP
		args.totalsnps = 30000 # Maximum requested SNPs
		args.scale = false # Flag to scale SNPs per contig by contig length
		args.scalehash = {} # Hash to get length from contig headers
		args.maxsnps = 2 # Maximum SNPs per contig
		args.distance = 10000 # Minimum distance between SNPs in a contig
		args.no_baits = false # Flag to omit generating baits
		args.refseq = "" # Reference sequence file
		args.lenbef = 60 # Length before SNP in bait
		args.lenaft = 59 # Length after SNP in bait
		args.tiling = false # Flag to determine whether to tile bait generation
		args.tileoffset = 20 # Offset between tiled baits
		args.tiledepth = 2 # Tiling depth
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
		args.coords = false # Flag to output coordinates table of baits
		args.no_indels = false # Flag to omit bait sequences with indels
		args.no_Ns = false # Flag to omit bait sequences with Ns
		args.haplodef = "" # Haplotype definition for aln2baits
		args.sort = false # Flag to sort stack2baits SNPs by between/within population variation
		args.hwe = false # Flag to sort stacks2baits SNPs by Hardy-Weinberg Equilibrium
		args.alpha = 0.05 # HWE test alpha value
		args.meanqual = 20.0 # Minimum mean base quality
		args.meanqual_filter = false # Flag to filter by minimum mean base quality
		args.minqual = 10 # Minimum base quality
		args.minqual_filter = false # Flag to filter by minimum base quality
		args.fasta_score = 0 # Asssumed base quality score for FASTA sequences
		args.features = [] # Array holding desired features
		args.ncbi = false # Flag whether FASTA/FASTQ headers include NCBI-style descriptors
		args.alt_alleles = false # Flag to apply alternate alleles
		args.varqual_filter = false # Flag to determine whether to filter vcf variants by QUAL scores
		args.varqual = 30 # Minimum vcf variant QUAL score
		opt_parser = OptionParser.new do |opts|
			if algorithms.include?(args.algorithm) # Process correct commands
				opts.banner = "Command-line usage: ruby baitstools.rb "+args.algorithm+" [options]"
				opts.separator ""
				opts.separator args.algorithm+"-specific options:"
				opts.separator ""
				case args.algorithm
				when "vcf2baits" #vcf2baits options
					opts.on("-i","--input [FILE]", String, "Input VCF file") do |vcf|
						args.infile = vcf
					end
					opts.on("-V", "--varqual [VALUE]", Integer, "Minimum variant QUAL score (Default = 30)") do |varf|
						args.varqual = varf if varf != nil
						args.varqual_filter = true
					end
					opts.on("-t", "--totalsnps [VALUE]", Integer, "Total requested SNPs (Default = 30,000)") do |tsnps|
						args.totalsnps = tsnps
					end
				when "stacks2baits" #stacks2baits only options
					opts.on("-i","--input [FILE]", String, "Input Stacks summary tsv") do |fa|
						args.infile = fa
					end
					opts.on("-S", "--sort", "Sort SNPs according to variation between/within populations") do
						args.sort = true
					end
					opts.on("-H", "--hwe", "Sort SNPs within populations according to Hardy-Weinberg Equilibrium (Implies -S)") do
						args.sort = true
						args.hwe = true
					end
					opts.on("-A", "--alpha", Float, "Set alpha value for HWE test (Either 0.10, 0.05, 0.025, 0.01) (Default = 0.05)") do |alph|
						args.alpha = alph
					end
					opts.on("-t", "--totalsnps [VALUE]", Integer, "Total requested SNPs per category (Default = 30,000)") do |tsnps|
						args.totalsnps = tsnps
					end
				else
					case args.algorithm
					when "coords2baits"
						inputstr = "Input coordinates table"
					when "annot2baits"
						inputstr = "Input GFF/GTF table"
					else 
						inputstr = "Input FASTA/FASTQ file"
					end
					opts.on("-i","--input [FILE]", String, inputstr) do |fa|
						args.infile = fa
					end
					if args.algorithm == "coords2baits" or args.algorithm == "annot2baits"
						opts.on("-r","--refseq [FILE]", String, "Reference FASTA/FASTQ sequence file") do |ref|
							args.refseq = ref
						end
					end
					opts.on("-L", "--length [VALUE]", Integer, "Requested bait length (Default = 120)") do |prblength|
						args.baitlength = prblength
					end
					unless args.algorithm == "checkbaits"
						opts.on("-O", "--offset [VALUE]", Integer, "Base pair offset between tiled baits (Default = 20)") do |toff|
							args.tileoffset = toff
						end
					end
					if args.algorithm == "annot2baits"
						opts.on("-G", "--features [FEATURE]", String, "Comma-separated list of feature types") do |feat|
							args.features = feat.upcase.split(",")
						end
					elsif args.algorithm == "aln2baits"
						opts.on("-H","--haplo [VALUE]", String, "Window haplotype definition (either 'haplotypes' or 'variant')") do |fa|
							args.haplodef = fa
						end
					end
				end	
				if args.algorithm == "vcf2baits" or args.algorithm == "stacks2baits" #vcf2baits, stacks2baits shared options
 					opts.on("-j", "--scale", "Scale maximum number of SNPs per contig by contig length (Overrides -m)") do
						args.scale = true
					end
					opts.on("-m","--maxsnps [VALUE]", Integer, "Maximum number of SNPs per contig (Default = 2)") do |msnps|
						args.maxsnps = msnps
					end
					opts.on("-d","--distance [VALUE]", Integer, "Minimum distance between SNPs within a contig (Default = 10,000)") do |dist|
						args.distance = dist
					end
					opts.on("-p","--no_baits", "Do not output bait sequences") do
						args.no_baits = true
					end
					 opts.on("-e","--every", "Output bait sequences for every SNP in the input file (Overrides -t, -j, -d, -m , -p)") do
						args.every = true
					end
					 opts.on("-R", "--alt", "Generate baits for alternate alleles (Overrides -p)") do
 						args.alt_alleles = true
 					end
					opts.on("-r","--refseq [FILE]", String, "If outputting baits, the reference FASTA/FASTQ sequence file") do |ref|
						args.refseq= ref
					end
					opts.on("-L", "--length [VALUE]", Integer, "Requested total bait length (Default = 120)") do |prblength|
						args.baitlength = prblength
						args.lenbef =	prblength - 1 #Sets these values to the appropriate position
						args.lenaft =  0 #Sets these values to the appropriate position
					end
					opts.on("-b", "--lenbef [VALUE]", Integer, "If outputting baits, the number of bases before the SNP (Default = 60) (Overrides -L)") do |b4|
						args.lenbef = b4
					end
					opts.on("-a", "--lenaft [VALUE]", Integer, "If outputting baits, the number of bases after the SNP (Default = 59) (Overrides -L)") do |af|
						args.lenaft = af
					end
					opts.on("-u", "--tiling", "If outputting baits, tile bait sequences") do
						args.tiling = true
					end
					opts.on("-O", "--offset [VALUE]", Integer, "If tiling baits, base pair offset between tiled baits (Default = 20) (Overrides -a and -b)") do |toff|
						args.tileoffset = toff
					end
					opts.on("-k", "--depth [VALUE]", Integer, "If tiling baits, requested baits per SNP") do |tdep|
						args.tiledepth = tdep
					end
				end
				opts.separator "" #bait filtering options
				opts.separator "Bait filtration options:"
				opts.separator ""
				opts.on("-w", "--params", "Output bait statistics table") do
					args.params = true
				end
				opts.on("-c","--complete", "Require baits be full length") do
					args.completebait = true
				end
				opts.on("-I", "--no_indels", "Exclude bait sequences with indels (-)") do
					args.no_indels = true
				end	
				opts.on("-N", "--no_Ns", "Exclude bait sequences with Ns") do
					args.no_Ns = true
				end	
				opts.on("-n","--mingc [VALUE]", Float, "Minimum GC content (Default = 30.0)") do |mgc|
					args.mingc = mgc if mgc != nil
					args.mingc_filter = true
				end
				opts.on("-x","--maxgc [VALUE]", Float, "Maximum GC content (Default = 50.0)") do |xgc|
					args.maxgc = xgc if xgc != nil
					args.maxgc_filter = true
				end			
				opts.on("-q","--mint [VALUE]", Float, "Minimum melting temperature (Default = 0.0)") do |mt|
					args.mint = mt if mt != nil
					args.mint_filter = true # Force bait filtration if called
				end
				opts.on("-z","--maxt [VALUE]", Float, "Maximum melting temperature (Default = 120.0)") do |xt|
					args.maxt = xt if xt != nil
					args.maxt_filter = true # Force bait filtration if called
				end
				opts.on("-T", "--type [VALUE]", String, "Melting temperature for DNA-DNA, RNA-RNA or RNA-DNA hybridization (Default = RNA-DNA)") do |hyb|
					args.bait_type = hyb if hyb != nil
				end
				opts.on("-s","--na [VALUE]", Float, "Melting temperature sodium concentration (Default = 0.9)") do |na|
					args.na = na if na != nil
				end
				opts.on("-f", "--formamide [VALUE]", Float, "Melting temperature formamide concentration (Default = 0.0)") do |form|
					args.formamide = form if form != nil
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
				unless args.algorithm == "checkbaits" # checkbaits does not get coordinate information to output
					opts.on("-o", "--coords", "Output coordinate tables for the candidate baits") do
						args.coords = true
					end
				end
				opts.on("-D", "--ncbi", "FASTA/FASTQ file headers have NCBI-style descriptions") do
					args.ncbi = true
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
						print "vcf2baits 0.7\n"
					when "tilebaits"
						print "tilebaits 0.3\n"
					when "coords2baits"
						print "coords2baits 0.1\n"
					when "annot2baits"
						print "annot2baits 0.1\n"
					when "aln2baits"
						print "aln2baits 0.2\n"
					when "stacks2baits"
						print "stacks2baits 0.1\n"
					when "checkbaits"
						print "checkbaits 0.1\n"
					end
					exit
				end
			else
				opts.banner = "Common options:"
				print "Welcome to baitstools 0.5.\n\nTo use the interactive interface, enter <ruby baitstools.rb [subcommand]> without command-line options.\nCommand-line usage: ruby baitstools.rb [subcommand] [options]"
				print "\nAdd '-h' or '--help' to subcommands (without other options) to see their relevant options.\n\nAvailable subcommands:\n\n"
				print "    aln2baits\t\t\t\tGenerate weighted baits from a FASTA/FASTQ alignment\n"
				print "    annot2baits\t\t\t\tGenerate tiled baits from a GFF/GTF file and a reference sequence\n"				
				print "    checkbaits\t\t\t\tFilter a FASTA/FASTQ of candidate baits by quality\n"
				print "    coords2baits\t\t\tGenerate tiled baits from a coordinates table and a reference sequence\n"
				print "    stacks2baits\t\t\tSelect variants and generate baits from a Stacks summary tsv file\n"
				print "    tilebaits\t\t\t\tGenerate tiled baits from FASTA/FASTQ sequences\n"
				print "    vcf2baits\t\t\t\tSelect variants and generate baits from a VCF\n\n"

			end
		end
		opt_parser.parse!(options)
		return args
	end
end
#-----------------------------------------------------------------------------------------------
begin
	# Build FASTQ quality translation hash
	fq_val = ["!","\"","#","$","%","&","\'","(",")","*","+",",","-",".","/","0","1","2","3","4","5","6","7","8","9",
		":",";","<","=",">","?","@","A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T",
		"U","V","W","X","Y","Z","[","\\","]","^","_","`","a","b","c","d","e","f","g","h","i","j","k","l","m","n","o",
		"p","q","r","s","t","u","v","w","x","y","z","{","|","}","~"]
	$fq_hash = {}
	for i in 0..93
		$fq_hash[fq_val[i]] = i
	end
	# Interactive mode block
	$options = Parser.parse(ARGV)
	exit if $options.kill
	if $options.interact
		print "Enter input file.\n"
		$options.infile = gets.chomp
	end
	while !FileTest.exist?($options.infile)
		print "Input file not found. Please re-enter.\n"
		$options.infile = gets.chomp
	end
	case $options.algorithm #algorithm-specific options
	when "vcf2baits","stacks2baits"
		if $options.interact
			print "Output baits for all variants in input file? (y/n)\n"
			t = gets.chomp.upcase
			if t == "Y" or t == "YES"
				$options.every = true
			end
			if $options.algorithm == "stacks2baits"
				print "Sort SNPs according to variation between/within populations? (y/n)\n"
				t = gets.chomp.upcase
				if t == "Y" or t == "YES"
					$options.sort = true
					print "Sort SNPs within populations according to Hardy-Weinberg Equilibrium? (y/n)\n"
					t = gets.chomp.upcase
					if t == "Y" or t == "YES"
						$options.hwe = true
						print "Alpha value for HWE (alpha = 0.10, 0.05, 0.025, 0.01)?\n"
						$options.alpha = gets.chomp.to_f
						alphas = [0.10, 0.05, 0.025, 0.01]
						while !alphas.include?($options.alpha)
							print "Alpha value must be 0.10, 0.05, 0.025, 0.01. Re-enter.\n"
							$options.alpha = gets.chomp.to_f
						end
					end
				end
			end
		end
		if $options.interact and !$options.every
			if $options.algorithm == "vcf2baits"
				print "Enter total number of requested variants.\n"
			else
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
			t = gets.chomp.upcase
			if t == "Y" or t == "YES"
				$options.scale = true
			end
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
			t = gets.chomp.upcase
			if t == "N" or t == "NO"
				$options.no_baits = true
			end
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
				t = gets.chomp.upcase
				if t == "Y" or t == "YES"
					$options.alt_alleles = true
				end
				print "Tile baits? (y/n)\n"
				t = gets.chomp.upcase
				if t == "Y" or t == "YES"
					$options.tiling = true
					print "Enter bait length.\n"
					$options.baitlength = gets.chomp.to_i
					while $options.baitlength < 1
						print "Baits must be at least 1 bp long. Re-enter.\n"
						$options.baitlength = gets.chomp.to_i
					end
					$options.lenbef = $options.baitlength - 1
					$options.lenaft = 0
					print "Enter tiling bp offset.\n"
					$options.tileoffset = gets.chomp.to_i
					print "Enter number of baits per SNP.\n"
					$options.tiledepth = gets.chomp.to_i
				end
			end
			if $options.interact and !$options.tiling
				print "Enter number of bases before variant in bait.\n"
				$options.lenbef = gets.chomp.to_i
			end
			while $options.lenbef < 0
				print "The number of bait bases before the variant must be at least 0. Please re-enter.\n"
				$options.lenbef = gets.chomp.to_i
			end
			if $options.interact and !$options.tiling
				print "Enter number of bases after variant in bait.\n"
				$options.lenaft = gets.chomp.to_i
			end
			while $options.lenaft < 0
				print "The number of bait bases after the variant must be at least 0. Please re-enter.\n"
				$options.lenaft = gets.chomp.to_i
			end
			$options.baitlength = $options.lenbef + $options.lenaft + 1
			if $options.tiling
				$options.lenbef = $options.baitlength - 1
				$options.lenaft = 0
				while $options.tileoffset > $options.baitlength or $options.tileoffset < 1
					print "Tiling offset cannot be less than 1 or greater than bait length. Please re-enter.\n"
					$options.tileoffset = gets.chomp.to_i
				end
				while $options.tiledepth > $options.baitlength/$options.tileoffset or $options.tiledepth < 1
					print "Tiling depth cannot be less than 1 or greater than bait length/tiling offset ratio. Re-enter.\n"
					$options.tiledepth = gets.chomp.to_i
				end
			end
		end
	else
		if $options.algorithm == "coords2baits" or $options.algorithm == "annot2baits"
			if $options.interact
				print "Enter reference sequence.\n"
				$options.refseq = gets.chomp
			end
			while !FileTest.exist?($options.refseq)
				print "Reference sequence not found. Please re-enter.\n"
				$options.refseq = gets.chomp
			end
		end
		if $options.interact
			print "Enter bait length.\n"
			$options.baitlength = gets.chomp.to_i
			while $options.baitlength < 1
				print "Baits must be at least 1 bp long. Re-enter.\n"
				$options.baitlength = gets.chomp.to_i
			end
			unless $options.algorithm == "checkbaits"
				print "Enter tiling bp offset.\n"
				$options.tileoffset = gets.chomp.to_i
				while $options.tileoffset < 1
					print "Tiling offset cannot be less than 1. Please re-enter.\n"
					$options.tileoffset = gets.chomp.to_i
				end
			end
			if $options.algorithm == "annot2baits"
				print "Enter comma-separated list of features.\n"
				$options.features = gets.chomp.upcase.split(",")
			end
		end
		if $options.algorithm == "aln2baits"
			if $options.interact
				print "Haplotype definition? ('haplotype' or 'variant')\n"
				$options.haplodef = gets.chomp
			end
			while $options.haplodef != 'haplotype' and $options.haplodef != 'variant'
				print "Please choose a haplotype definition ('haplotype' or 'variant')\n"
				$options.haplodef = gets.chomp
			end
		end	
	end
	if $options.interact
		print "Do FASTA/FASTQ sequence headers include NCBI-style descriptions? (y/n)\n"
		t = gets.chomp.upcase
		if t == "Y" or t == "YES"
			$options.ncbi = true
		end
		if !$options.no_baits
			if $options.algorithm != "checkbaits"
				print "Output coordinates table(s) for candidate baits? (y/n)\n"
				t = gets.chomp.upcase
				if t == "Y" or t == "YES"
					$options.coords = true
				end
			end
			print "Output bait statistics table? (y/n)\n"
			t = gets.chomp.upcase
			if t == "Y" or t == "YES"
				$options.params = true
			end
			print "Require complete length bait? (y/n)\n"
			t = gets.chomp.upcase
			if t == "Y" or t == "YES"
				$options.completebait = true
			end
			print "Exclude baits with indels (-)? (y/n)\n"
			t = gets.chomp.upcase
			if t == "Y" or t == "YES"
				$options.no_indels = true
			end
			print "Exclude baits with Ns? (y/n)\n"
			t = gets.chomp.upcase
			if t == "Y" or t == "YES"
				$options.no_Ns = true
			end
			print "Filter by minimum GC content? (y/n)\n"
			t = gets.chomp.upcase
			if t == "Y" or t == "YES"
				$options.mingc_filter = true
				print "Enter minimum GC content.\n"
				$options.mingc = gets.chomp.to_f
				while $options.mingc < 0.0 or $options.mingc > 100.0
					print "Minimum GC content must be between 0 and 100. Please re-enter.\n"
					$options.mingc = gets.chomp.to_f
				end
			end
			print "Filter by maximum GC content? (y/n)\n"
			t = gets.chomp.upcase
			if t == "Y" or t == "YES"
				$options.maxgc_filter = true
				print "Enter maximum GC content.\n"
				$options.maxgc = gets.chomp.to_f
				while $options.maxgc < $options.mingc or $options.maxgc > 100.0
					print "Maximum GC content must be between minimum GC content and 100. Re-enter.\n"
					$options.maxgc = gets.chomp.to_f
				end
			end
			print "Filter by minimum melting temperature? (y/n)\n"
			t = gets.chomp.upcase
			if t == "Y" or t == "YES"
				$options.mint_filter = true
				print "Enter minimum melting temperature.\n"
				$options.mint = gets.chomp.to_f
			end
			print "Filter by maximum melting temperature? (y/n)\n"
			t = gets.chomp.upcase
			if t == "Y" or t == "YES"
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
	end
	if ($options.mint_filter or $options.maxt_filter)
		if $options.interact
			$options.bait_type = "" # Turn off bait_type default so that interactive mode can choose
			print "Enter sodium concentration.\n"
			$options.na = gets.chomp.to_f
			while $options.na < 0.0
				print "Sodium concentrations cannot be negative. Re-enter.\n"
				$options.na = gets.chomp.to_f
			end
		end
		while $options.bait_type != 'RNA-RNA' and $options.bait_type != 'DNA-DNA' and $options.bait_type != 'RNA-DNA'
			print "Please choose a hybridization type ('DNA-DNA', 'RNA-RNA' or 'RNA-DNA')\n"
			$options.bait_type = gets.chomp
		end
	end
	if $options.interact and $options.algorithm != "aln2baits" and !$options.no_baits
		print "Filter by minimum mean base quality? (y/n)\n"
		t = gets.chomp.upcase
		if t == "Y" or t == "YES"
			$options.meanqual_filter = true
			print "Enter minimum mean base quality.\n"
			$options.meanqual = gets.chomp.to_f
		end
		print "Filter by minimum base quality? (y/n)\n"
		t = gets.chomp.upcase
		if t == "Y" or t == "YES"
			$options.minqual_filter = true
			print "Enter minimum mean base quality.\n"
			$options.min = gets.chomp.to_i
		end
		if ($options.meanqual_filter or $options.minqual_filter)
			print "Enter assumed FASTA quality score.\n"
			$options.fasta_score = gets.chomp.to_i
		end
	end
	$options.fasta_score = 0 if $options.fasta_score < 0 # Correct fasta_score parameter
	$options.fasta_score = 93 if $options.fasta_score > 93 # Limit to currently possible scores
	$options.no_baits = false if ($options.every or $options.alt_alleles) # Override -p when needed
	$options.filter = true if ($options.completebait or $options.params or $options.algorithm == "checkbaits" or $options.mingc_filter or $options.maxgc_filter or $options.mint_filter or $options.maxt_filter or $options.meanqual_filter or $options.minqual_filter or $options.no_indels or $options.no_Ns) # Force filtration as necessary
	cmdline = get_command_line
	print "** Starting program with the following options: **\n"
	print "** Basic command: " + cmdline[0] + " **\n"
	print "** Filtration options:" + cmdline[1] + " **\n" # filtered line always starts with a space if present
	case $options.algorithm
	when "aln2baits"
		aln2baits
	when "annot2baits"
		annot2baits
	when "checkbaits"
		checkbaits
	when "coords2baits"
		coords2baits
	when "stacks2baits"
		stacks2baits
	when "tilebaits"
		tilebaits($options.infile)
	when "vcf2baits"
		vcf2baits
	end
	print "** Program complete **\n"
end
