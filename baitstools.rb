#!/usr/bin/ruby
#-----------------------------------------------------------------------------------------------
# baitstools 0.3
# Michael G. Campana, 2016
# Smithsonian Conservation Biology Institute
#-----------------------------------------------------------------------------------------------

# To be implemented
# Output non-reference alleles
# Include/Exclude multi-allelic sites

require 'optparse'
require 'ostruct'
require './baitslib.rb'
require './selectsnps.rb'
require './tilebaits.rb'
require './coords2baits.rb'
require './aln2baits.rb'
#-----------------------------------------------------------------------------------------------
class Parser
	def self.parse(options)
		# Set defaults
		args = OpenStruct.new
		algorithms = ["selectsnps", "tilebaits","coords2baits", "aln2baits"]
		args.algorithm = ARGV[0] # Define subcommands
		if algorithms.include?(args.algorithm) #Set interactive mode or whether to kill
			args.interact = true if ARGV.size == 1
			args.kill = false
		else
			args.kill = true #kill the program if neither interactive nor command line
		end
		ARGV.delete_at(0) #Remove subcommand from input otherwise crashes
		args.infile = ""
		args.every = false
		args.totalsnps = 30000
		args.scale = false
		args.scalehash = {}
		args.maxsnps = 2
		args.distance = 10000
		args.baits = false
		args.refseq = ""
		args.lenbef = 60
		args.lenaft = 59
		args.tiling = false
		args.tileoffset = 20
		args.tiledepth = 2
		args.filter = false
		args.completebait = false
		args.baitlength = 120
		args.gc = false
		args.mingc = 30.0
		args.maxgc = 50.0
		args.melt = false
		args.na = 0.9
		args.mint = 0.0
		args.maxt = 120.0
		args.params = false
		args.coords = false
		args.no_indels = false
		args.haplodef = "haplotype"
		opt_parser = OptionParser.new do |opts|
			if algorithms.include?(args.algorithm) #Process correct commands
				opts.banner = "Command-line usage: ruby baitstools.rb "+args.algorithm+" [options]"
				opts.separator ""
				opts.separator args.algorithm+"-specific options:"
				opts.separator ""
				case args.algorithm
				when "selectsnps" #selectsnps options
					opts.on("-i","--input [FILE]", String, "Input VCF file") do |vcf|
						args.infile = vcf
					end
					opts.on("-e","--every", "Output bait sequences for every SNP in the input VCF file") do
						args.every = true
						args.baits = true
					end
						opts.on("-t", "--totalsnps [VALUE]", Integer, "Total requested SNPs (Default = 30,000)") do |tsnps|
						args.totalsnps = tsnps
					end
					opts.on("-j", "--scale", "Scale maximum number of SNPs per contig by contig length (overrides -m)") do
						args.scale = true
					end
					opts.on("-m","--maxsnps [VALUE]", Integer, "Maximum number of SNPs per contig (Default = 2)") do |msnps|
						args.maxsnps = msnps
					end
					opts.on("-d","--distance [VALUE]", Integer, "Minimum distance between SNPs within a contig (Default = 10,000)") do |dist|
						args.distance = dist
					end
					opts.on("-p","--baits", "Output bait sequences") do
						args.baits = true
					end
					opts.on("-r","--refseq [FILE]", String, "If outputting baits, the reference FASTA sequence file") do |ref|
						args.refseq= ref
					end
					opts.on("-L", "--length [VALUE]", Integer, "Requested total bait length (Default = 120).") do |prblength|
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
					opts.on("-u", "--tiling", "If outputting baits, tile bait sequences.") do
						args.tiling = true
					end
					opts.on("-O", "--offset [VALUE]", Integer, "If tiling baits, base pair offset between tiled baits (Default = 20). (Overrides -a and -b)") do |toff|
						args.tileoffset = toff
					end
					opts.on("-k", "--depth [VALUE]", Integer, "If tiling baits, requested baits per SNP.") do |tdep|
						args.tiledepth = tdep
					end
				when "tilebaits" #tilebaits only options
					opts.on("-i","--input [FILE]", String, "Input FASTA file") do |fa|
						args.infile = fa
					end
					opts.on("-L", "--length [VALUE]", Integer, "Requested bait length (Default = 120).") do |prblength|
						args.baitlength = prblength
					end
					opts.on("-O", "--offset [VALUE]", Integer, "Base pair offset between tiled baits (Default = 20).") do |toff|
						args.tileoffset = toff
					end
				when "coords2baits" #coords2baits only options
					opts.on("-i","--input [FILE]", String, "Input coordinates table") do |fa|
						args.infile = fa
					end
					opts.on("-r","--refseq [FILE]", String, "Reference FASTA sequence file") do |ref|
						args.refseq = ref
					end
					opts.on("-L", "--length [VALUE]", Integer, "Requested bait length (Default = 120).") do |prblength|
						args.baitlength = prblength
					end
					opts.on("-O", "--offset [VALUE]", Integer, "Base pair offset between tiled baits (Default = 20).") do |toff|
						args.tileoffset = toff
					end
				when "aln2baits" #aln2baits only options
					opts.on("-i","--input [FILE]", String, "Input FASTA alignment") do |fa|
						args.infile = fa
					end
					opts.on("-L", "--length [VALUE]", Integer, "Requested bait length (Default = 120).") do |prblength|
						args.baitlength = prblength
					end
					opts.on("-O", "--offset [VALUE]", Integer, "Base pair offset between haplotype windows (Default = 20).") do |toff|
						args.tileoffset = toff
					end
					opts.on("-N", "--no_indels", "Exclude sequences with indels.") do
						args.no_indels = true
					end
					opts.on("-H","--haplo [VALUE]", String, "Window haplotype definition (either 'haplotypes' or 'variants')") do |fa|
						args.haplodef = fa
					end
				end	
				opts.separator "" #bait filtering options
				opts.separator "Bait filtration options:"
				opts.separator ""
				opts.on("-f","--filter", "If outputting baits, filter resulting bait sequences") do
					args.filter = true
				end
				opts.on("-w", "--params", "If filtering, output bait statistics table") do
					args.params = true
				end
				opts.on("-c","--complete", "If filtering, require baits be full length") do
					args.completebait = true
				end
				opts.on("-g", "--gc", "If filtering, select_snpsfilter by GC content") do
					args.gc = true
				end				
				opts.on("-n","--mingc [VALUE]", Float, "If filtering, minimum GC content (Default = 30.0)") do |mgc|
					args.mingc = mgc
				end
				opts.on("-x","--maxgc [VALUE]", Float, "If filtering, maximum GC content (Default = 50.0)") do |xgc|
					args.maxgc = xgc
				end
				opts.on("-l", "--melt", "If filtering, filter by melting temperature") do
					args.melt = true
				end				
				opts.on("-q","--mint [VALUE]", Float, "If filtering, minimum melting temperature (Default = 0.0)") do |mt|
					args.mint = mt
				end
				opts.on("-z","--maxt [VALUE]", Float, "If filtering, maximum melting temperature (Default = 120.0)") do |xt|
					args.maxt = xt
				end
				opts.on("-s","--na [VALUE]", Float, "If filtering by melting temperature, sodium concentration (Default = 0.9)") do |na|
					args.na = na
				end
				opts.on("-o", "--coords", "Output coordinate tables for the candidate baits.") do
					args.coords = true
				end
				opts.separator ""
				opts.separator "Common options:"
				opts.on_tail("-h","--help", "Show help") do
					print "Welcome to baitstools " + args.algorithm + '.' +"\n\n"
					print "To use the interactive interface, enter <ruby baitstools.rb " + args.algorithm + "> without command-line options.\n"
					puts opts
					exit
				end
				opts.on_tail("-v","--version","Show subcommand version") do
					case args.algorithm
					when "selectsnps"
						print "selectsnps 0.6\n"
					when "tilebaits"
						print "tilebaits 0.3\n"
					when "coords2baits"
						print "coords2baits 0.1\n"
					when "aln2baits"
						print "aln2baits 0.2\n"
					end
					exit
				end
			else
				if !ARGV.include?("-v") and !ARGV.include?("--version")
					opts.banner = "Common options:"
					print "Welcome to baitstools.\n\nTo use the interactive interface, enter <ruby baitstools.rb [subcommand]> without command-line options.\nCommand-line usage: ruby baitstools.rb [subcommand] [options]"
					print "\nAdd '-h' or '--help' to subcommands (without other options) to see their relevant options.\n\nAvailable subcommands:\n\n"
					print "    selectsnps\t\t\t\tSelect SNPs from a VCF\n"
					print "    tilebaits\t\t\t\tGenerate tiled baits from FASTA sequences\n"
					print "    coords2baits\t\t\tGenerate tiled baits from a coordinates table and a reference sequence\n"
					print "    aln2baits\t\t\t\tGenerate weighted baits from a FASTA alignment\n\n"
				end
				opts.on_tail("-h","--help", "Show help") do
					puts opts
					exit
				end				
				opts.on_tail("-v","--version","Show baitstools version") do
					print "baitstools 0.3\n"
					exit
				end
			end
		end
		opt_parser.parse!(options)
		return args
	end
end
#-----------------------------------------------------------------------------------------------
begin
	# Interactive mode block
	$options = Parser.parse(ARGV)
	exit if $options.kill
	if $options.interact
		print "Enter input file.\n"
		$options.infile = gets.chomp
	end
	while !FileTest.exist?($options.infile)
		print "File not found. Please re-enter.\n"
		$options.infile = gets.chomp
	end
	case $options.algorithm #algorithm-specific options
	when "selectsnps" 
		if $options.interact
			print "Output baits for all SNPs in VCF? (y/n)\n"
			t = gets.chomp.upcase
			if t == "Y" or t == "YES"
				$options.every = true
				$options.baits = true
			end
		end
		if $options.interact and !$options.every
			print "Enter total number of requested SNPs.\n"
			$options.totalsnps = gets.chomp.to_i
		end
		while $options.totalsnps < 1
			print "The total number of SNPs must be greater than 0. Please re-enter.\n"
			$options.totalsnps = gets.chomp.to_i
		end
		if $options.interact and !$options.every
			print "Scale maximum number of SNPs per contig by contig length? (y/n)\n"
			t = gets.chomp.upcase
			if t == "Y" or t == "YES"
				$options.scale = true
			end
		end
		if $options.interact and !$options.every and !$options.scale
			print "Enter maximum number of SNPs per contig.\n"
			$options.maxsnps = gets.chomp.to_i
		end
		while $options.maxsnps < 1 and !$options.scale
			print "The maximum number of SNPs per contig must be greater than 0. Please re-enter.\n"
			$options.maxsnps = gets.chomp.to_i
		end
		if $options.interact and !$options.every
			print "Enter minimum distance between SNPs within a contig.\n"
			$options.distance = gets.chomp.to_i
		end
		while $options.distance < 1
			print "The minimum distance between SNPs must be greater than 0. Please re-enter.\n"
			$options.distance = gets.chomp.to_i
		end
		if $options.interact and !$options.every
			print "Output baits? (y/n)\n"
			t = gets.chomp.upcase
			if t == "Y" or t == "YES"
				$options.baits = true
			end
		end
		if $options.baits
			if $options.interact
				print "Enter reference sequence.\n"
				$options.refseq = gets.chomp
			end
			while !FileTest.exist?($options.refseq)
				print "Reference sequence not found. Please re-enter.\n"
				$options.refseq = gets.chomp
			end
			if $options.interact
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
				print "Enter number of bases before SNP in bait.\n"
				$options.lenbef = gets.chomp.to_i
			end
			while $options.lenbef < 0
				print "The number of bait bases before the SNP must be at least 0. Please re-enter.\n"
				$options.lenbef = gets.chomp.to_i
			end
			if $options.interact and !$options.tiling
				print "Enter number of bases after SNP in bait.\n"
				$options.lenaft = gets.chomp.to_i
			end
			while $options.lenaft < 0
				print "The number of bait bases after the SNP must be at least 0. Please re-enter.\n"
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
	when "tilebaits"
		if $options.interact
			print "Enter bait length.\n"
			$options.baitlength = gets.chomp.to_i
			while $options.baitlength < 1
				print "Baits must be at least 1 bp long. Re-enter.\n"
				$options.baitlength = gets.chomp.to_i
			end
			print "Enter tiling bp offset.\n"
			$options.tileoffset = gets.chomp.to_i
			while $options.tileoffset < 1
				print "Tiling offset cannot be less than 1. Please re-enter.\n"
				$options.tileoffset = gets.chomp.to_i
			end
		end
	when "coords2baits"
		if $options.interact
			print "Enter reference sequence.\n"
			$options.refseq = gets.chomp
		end
		while !FileTest.exist?($options.refseq)
			print "Reference sequence not found. Please re-enter.\n"
			$options.refseq = gets.chomp
		end
		if $options.interact
			print "Enter bait length.\n"
			$options.baitlength = gets.chomp.to_i
			while $options.baitlength < 1
				print "Baits must be at least 1 bp long. Re-enter.\n"
				$options.baitlength = gets.chomp.to_i
			end
			print "Enter tiling bp offset.\n"
			$options.tileoffset = gets.chomp.to_i
			while $options.tileoffset < 1
				print "Tiling offset cannot be less than 1. Please re-enter.\n"
				$options.tileoffset = gets.chomp.to_i
			end
		end
	when "aln2baits"
		if $options.interact
			print "Enter bait length.\n"
			$options.baitlength = gets.chomp.to_i
			while $options.baitlength < 1
				print "baits must be at least 1 bp long. Re-enter.\n"
				$options.baitlength = gets.chomp.to_i
			end
			print "Enter window bp offset.\n"
			$options.tileoffset = gets.chomp.to_i
			while $options.tileoffset < 1
				print "Window offset cannot be less than 1. Please re-enter.\n"
				$options.tileoffset = gets.chomp.to_i
			end
			print "Exclude indel sites? (y/n)\n"
			t = gets.chomp.upcase
			if t == "Y" or t == "YES"
				$options.no_indels = true
			end
			print "Haplotype definition? ('haplotype' or 'variants')\n"
			$options.haplodef = gets.chomp
		end
		while $options.haplodef != 'haplotype' and $options.haplodef != 'variants'
			print "Please choose a haplotype definition ('haplotype' or 'variants')\n"
			$options.haplodef = gets.chomp
		end
	end
	if $options.interact
		print "Filter baits? (y/n)\n"
		t = gets.chomp.upcase
		if t == "Y" or t == "YES"
			$options.filter = true
		end
		print "Output coordinates table(s) for candidate baits? (y/n)\n"
		t = gets.chomp.upcase
		if t == "Y" or t == "YES"
			$options.coords = true
		end
	end
	if $options.filter
		if $options.interact
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
			print "Filter by GC content? (y/n)\n"
			t = gets.chomp.upcase
			if t == "Y" or t == "YES"
				$options.gc = true
			end
			print "Filter by melting temperature? (y/n)\n"
			t = gets.chomp.upcase
			if t == "Y" or t == "YES"
				$options.melt = true
			end
		end
		if $options.gc
			if $options.interact
				print "Enter minimum GC content.\n"
				$options.mingc = gets.chomp.to_f
			end
			while $options.mingc < 0.0 or $options.mingc > 100.0
				print "Minimum GC content must be between 0 and 100. Please re-enter.\n"
				$options.mingc = gets.chomp.to_f
			end
			if $options.interact
				print "Enter maximum GC content.\n"
				$options.maxgc = gets.chomp.to_f
			end
			while $options.maxgc < $options.mingc or $options.maxgc > 100.0
				print "Maximum GC content must be between minimum GC content and 100. Re-enter.\n"
				$options.maxgc = gets.chomp.to_f
			end
		end
		if $options.melt
			if $options.interact
				print "Enter minimum melting temperature.\n"
				$options.mint = gets.chomp.to_f
				print "Enter maximum melting temperature.\n"
				$options.maxt = gets.chomp.to_f
			end
			while $options.maxt < $options.mint
				print "Maximum melting temperature cannot be less than minimum melting temperature. Re-enter.\n"
				$options.maxt = gets.chomp.to_f
			end
			if $options.interact
				print "Enter sodium concentration.\n"
				$options.na = gets.chomp.to_f
			end
			while $options.na < 0.0
				print "Sodium concentrations cannot be negative. Re-enter.\n"
				$options.na = gets.chomp.to_f
			end
		end
	end
	case $options.algorithm
	when "selectsnps"
		selectsnps
	when "tilebaits"
		tilebaits($options.infile)
	when "coords2baits"
		coords2baits
	when "aln2baits"
		aln2baits
	end
end
