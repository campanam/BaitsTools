#!/usr/bin/ruby
#-----------------------------------------------------------------------------------------------
# baitstools
BAITSTOOLSVER = "0.6"
# Michael G. Campana, 2017
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

#-----------------------------------------------------------------------------------------------
class Fa_Seq #container for fasta/fastq sequences
	attr_accessor :header, :circular, :fasta, :seq, :qual, :desc, :qual_array
	def initialize(header = "", circular = false, fasta = false, seq = "", qual = "")
		@header = header # Sequence header
		@circular = circular # Circular sequence flag
		@seq = seq # DNA sequence
		@qual = qual # Original quality scores
		@qual_array = [] # Array of numeric quality scores
		@fasta = fasta # FASTA format flag
	end
	def calc_quality # Convert quality scores to numeric values so only needed once
		for i in 0...@qual.length
			@qual_array.push($fq_hash[@qual[i]])
		end
	end
	def numeric_quality
		return @qual_array
	end
end
#-----------------------------------------------------------------------------------------------
class Chromo_SNP # Container for SNPs
	attr_accessor :chromo, :snp, :popvar_data, :ref, :alt, :qual
	def initialize(chromo, snp, popvar_data =[], ref = nil, alt = nil, qual = nil, line = nil)
		@chromo = chromo # Chromosome
		@snp = snp # SNP index in 1-based indexing
		@popvar_data = popvar_data # Array holding stacks population-specific allele data
		@ref = ref # Reference allele (vcf)
		@alt = alt # Alternate alleles (vcf)
		@qual = qual # Allele quality for vcf
		@line = line # Original datalines for recall
	end
	def within_pops? # Determine whether SNP variable within populations
		within_pops = false
		for pop in @popvar_data
			if !pop.monomorphic?
				within_pops = true
				break
			end
		end
		return within_pops
	end
	def line
		if @line.nil? # Combine data from stacks2baits Popvars if line not previously defined
			@line = ""
			for pop in @popvar_data
				@line += pop.line
			end
		end
		return @line
	end
	def alt_alleles 
		if @alt.nil? # Get alternate alleles for stacks2baits (once)
			@alt = []
			for pop in @popvar_data
				@alt += pop.alleles
			end
			@alt.uniq!
		end
		return @alt
	end
end
#-----------------------------------------------------------------------------------------------
def mean(val = [])
	mean = 0.0
	for i in 0.. val.size-1
		mean += val[i]
	end
	mean /= val.size
	return mean
end
#-----------------------------------------------------------------------------------------------
def filter_baits(bait, qual = [0])

	# To be implemented:
	# Sequence complexity filter
	# Self-complementarity filter

	keep = true
	keep = false if (bait.length < $options.baitlength && $options.completebait)
	keep = false if (bait.upcase.include?("N") and $options.no_Ns)
	keep = false if (bait.upcase.include?("-") and $options.no_gaps)
	gc = 0.0
	for i in 0 ... bait.length
		if bait[i].chr.upcase == "G" or bait[i].chr.upcase == "C"
			gc += 1.0
		end
	end
	gccont = gc/bait.length.to_f
	case $options.bait_type
	when "RNA-RNA"
		melt = 79.8 + (58.4 * gccont) + (11.8 * (gccont**2.0)) - 820.0/bait.length.to_f + 18.5 * Math::log($options.na) - 0.35 * $options.formamide
	when "RNA-DNA"
		melt = 79.8 + (58.4 * gccont) + (11.8 * (gccont**2.0)) - 820.0/bait.length.to_f + 18.5 * Math::log($options.na) - 0.50 * $options.formamide
	when "DNA-DNA"
		melt = 81.5 + (41.0 * gccont) - 500.0/bait.length.to_f + 16.6 * Math::log($options.na) - 0.62 * $options.formamide
	end	
	meanqual = mean(qual)
	minqual = qual.min
	keep = false if (gccont * 100.0 < $options.mingc && $options.mingc_filter)
	keep = false if (gccont * 100.0 > $options.maxgc && $options.maxgc_filter)
	keep = false if (melt < $options.mint && $options.mint_filter)
	keep = false if (melt > $options.maxt && $options.maxt_filter)
	keep = false if (meanqual < $options.meanqual && $options.meanqual_filter)
	keep = false if (minqual < $options.minqual && $options.minqual_filter)
	return [keep, bait.length.to_s + "\t" + (gccont * 100.0).to_s + "\t" + melt.to_s + "\t" + meanqual.to_s + "\t" + minqual.to_s + "\t" + keep.to_s + "\n"]
end
#-----------------------------------------------------------------------------------------------
def write_baits(baitsout = "", outfilter = "", paramline = "", coordline = "", filtercoordline = "", filestem = $options.infile)
	unless $options.algorithm == "checkbaits"
		File.open(filestem + "-baits.fa", 'w') do |write|
			write.puts baitsout
		end
	end
	if $options.coords
		File.open(filestem + "-baits.bed", 'w') do |write|
			write.puts coordline
		end
	end
	if $options.filter
		File.open(filestem + "-filtered-baits.fa", 'w') do |write|
			write.puts outfilter
		end
		if $options.params
			File.open(filestem + "-filtered-params.txt", 'w') do |write|
				write.puts paramline
			end
		end
		if $options.coords
			File.open(filestem + "-filtered-baits.bed", 'w') do |write|
				write.puts filtercoordline
			end
		end
	end
end
#-----------------------------------------------------------------------------------------------
def read_fasta(file) # Read fasta and fastq files
	seq_array = []
	faseq = nil # Dummy value
	qual = false # Denotes whether sequence or quality string
	File.open(file, 'r') do |seq|
		while line = seq.gets
			unless line == "\n" # Remove extraneous line breaks
				if line[0].chr == ">" and qual == false
					if !faseq.nil? # push previously completed sequence into array
						faseq.calc_quality unless faseq.fasta # Calculate quality scores
						seq_array.push(faseq)
					end
					header = line[1...-1] # Remove final line break and beginning >
					if header[-5..-1] == "#circ" #Adding this to the end of a sequence denotes as circular
						circular = true
						header = header[0...-5] # Remove circ hashtag
					else
						circular = false
					end
					if $options.ncbi # Remove sequence description from sequence name for NCBI-formatted files
						header_array = header.split(" ")
						header = header_array[0]
					end
					faseq = Fa_Seq.new(header, circular, true)
				elsif qual and faseq.qual.length < faseq.seq.length # Test if quality line is complete
					faseq.qual += line[0...-1] #Exclude final line break, allow multi-line fastq
			 	elsif line[0].chr == "@" # Start new sequence 	
					qual = false
					if !faseq.nil? # push previously completed sequence into array
						faseq.calc_quality unless faseq.fasta # Calculate quality scores
						seq_array.push(faseq)
					end 
					header = line[1...-1] # Remove final line break and beginning >
					if header[-5..-1] == "#circ" #Adding this to the end of a sequence denotes as circular
						circular = true
						header = header[0...-5] # Remove circ hashtag
					else
						circular = false
					end
					if $options.ncbi # Remove sequence description from sequence name for NCBI-formatted files
						header_array = header.split(" ")
						header = header_array[0]
					end
					faseq = Fa_Seq.new(header, circular, false)
				elsif line[0].chr == "+"
					qual = true
				else
					faseq.seq += line[0...-1] #Exclude final line break, allow multi-line fasta/fastq
				end
			end
		end
	end
	faseq.calc_quality unless faseq.fasta # Calculate quality scores
	seq_array.push(faseq) # Put last sequence into fasta array
	return seq_array
end
#-----------------------------------------------------------------------------------------------
def selectsnps(snp_hash) # Choose SNPs based on input group of SNPSs
	temp_snps = snp_hash.dup #Avoid messing with the original hash.
	selectsnps = {}
	if !$options.every
		for i in 1..$options.totalsnps
			break if temp_snps.size == 0 # Stop selecting SNPs when all contigs deleted from consideration
			selected_contig = temp_snps.keys[rand(temp_snps.size)] # Get name of contig
			selected_snp = temp_snps[selected_contig][rand(temp_snps[selected_contig].size)]
			# Delete SNPs that are too close
			tmp = temp_snps[selected_contig].dup # Duplicate this subsection so that during deletion there are not errors
			for i in 0 .. tmp.size - 1
				# 0 slot of selected_snp and temp_snp is SNP id, other slots in array are metadata
				if tmp[i].snp < selected_snp.snp && selected_snp.snp - tmp[i].snp < $options.distance
					temp_snps[selected_contig].delete(tmp[i])
				elsif tmp[i].snp > selected_snp.snp && tmp[i].snp - selected_snp.snp < $options.distance
					temp_snps[selected_contig].delete(tmp[i])
				end
			end
			# Add SNP to selected pool and delete contigs if maximum number of SNPs reached or no remaining SNPs on contig
			selectsnps[selected_contig] = [] if selectsnps[selected_contig].nil?
			selectsnps[selected_contig].push(selected_snp)
			temp_snps[selected_contig].delete(selected_snp) # So it cannot be reselected
			if $options.scale
				maxsize = $options.scalehash[selected_contig]
			else
				maxsize = $options.maxsnps
			end
			if selectsnps[selected_contig].size == maxsize or temp_snps[selected_contig].size == 0
				temp_snps.delete_if {|key, value | key == selected_contig}
			end
		end
	else
		selectsnps = snp_hash
	end
	return selectsnps
end
#-----------------------------------------------------------------------------------------------
def snp_to_baits(selectedsnps, refseq)
	baitsout = ""
	outfilter = ""
	paramline = "Chromosome:Coordinates\tSNP\tBaitLength\t%GC\tTm\tMeanQuality\tMinQuality\tKept\n"
	coordline = ""
	filtercoordline = ""
	filteredsnps = selectedsnps.dup
	for rseq in refseq
		if selectedsnps.keys.include?(rseq.header)
			for snp in selectedsnps[rseq.header]
				if $options.tiling
					tiling = $options.tiledepth
				else
					tiling = 1
					$options.tileoffset = 1
				end
				rseq_vars = [rseq]
				if $options.alt_alleles
					for altvar in snp.alt_alleles
						tmpseq = rseq.seq.dup # Need to duplicate to prevent changes to original sequence
						if snp.ref.nil?
							tmpseq[snp.snp-1]=altvar
						else
							tmpseq[snp.snp-1..snp.snp+snp.ref.length-2]=altvar
						end
						tmpheader = rseq.header.dup + "-" + altvar + "-variant"
						tmp = Fa_Seq.new(tmpheader, rseq.circular, rseq.fasta, tmpseq, rseq.qual)
						unless rseq.fasta
							qual_arr = []
							tmpqual = rseq.qual_array.dup
							if $options.algorithm == "vcf2baits"
								qual = snp.qual
								qual = 40 if qual > 40 # Adjust for quality scores beyond typical sequence capability
							else
								qual = rseq.numeric_quality[snp.snp-1] # Assume quality is equal across the variants since no quality information
							end
							for i in 1..altvar.length # Apply quality score to all bases in alternate allele
								qual_arr.push(qual)
							end
							tmpqual[snp.snp-1]=qual_arr
							tmp.qual_array = tmpqual.flatten
						end
						rseq_vars.push(tmp) # Add alternate variants to the sequence
					end
				end
				all_removed = 0 # Determine whether to remove whole site from filtered SNPs
				for rseq_var in rseq_vars
					for tile in 0...tiling
						be4 = snp.snp - $options.lenbef + tile * $options.tileoffset
						after = snp.snp + $options.lenaft + tile * $options.tileoffset
						be4 = 1 if (be4 < 1 and !rseq_var.circular)
						qual = [$options.fasta_score]
						if after > rseq_var.seq.length and !rseq_var.circular
							after = rseq_var.seq.length
						if be4 < 1
							prb = rseq_var.seq[be4-1..-1] + rseq_var.seq[0..after-1]  #Correct for 0-based counting
							qual = rseq_var.numeric_quality[be4-1..-1] + rseq_var.numeric_quality[0..after-1] unless rseq_var.fasta
						else
							prb = rseq_var.seq[be4-1..after-1]  #Correct for 0-based counting
							qual = rseq_var.numeric_quality[be4-1..after-1] unless rseq_var.fasta
						end 
						elsif after > rseq_var.seq.length and rseq_var.circular
							after -= rseq_var.seq.length
							if be4 < 1
								prb = rseq_var.seq[be4-1..-1] + rseq_var.seq[be4-1..rseq_var.seq.length-1]+rseq_var.seq[0..after-1] #Correct for 0-based counting and end of sequence
								qual = rseq_var.numeric_quality[be4-1..-1] + rseq_var.numeric_quality[be4-1..rseq_var.seq.length-1]+rseq_var.numeric_quality[0..after-1] unless rseq_var.fasta
							else
								prb = rseq_var.seq[be4-1..rseq_var.seq.length-1]+rseq_var.seq[0..after-1]  #Correct for 0-based counting and end of sequence
								qual = rseq_var.numeric_quality[be4-1..rseq_var.seq.length-1]+rseq.numeric_quality[0..after-1] unless rseq_var.fasta
							end
						else
							if be4 < 1
								prb = rseq_var.seq[be4-1..-1] + rseq_var.seq[0..after-1]
								qual = rseq_var.numeric_quality[be4-1..-1] + rseq_var.numeric_quality[0..after-1] unless rseq_var.fasta
							else
								prb = rseq_var.seq[be4-1..after-1]  #Correct for 0-based counting
								qual = rseq_var.numeric_quality[be4-1..after-1] unless rseq_var.fasta
							end
						end
						seq = ">" + rseq_var.header + "\t" + snp.snp.to_s + "\n" + prb + "\n"
						baitsout += seq
						be4 = rseq_var.seq.length + be4 if be4 < 1
						coord = rseq_var.header + "\t" + be4.to_s + "\t" + after.to_s + "\n"
						coordline += coord
						if $options.filter
							parameters = filter_baits(prb, qual)
							if parameters[0]
								outfilter += seq
								filtercoordline += coord 
							else 
								all_removed += 1
								filteredsnps[rseq.header].delete(snp) if all_removed == rseq_vars.size # Only remove completely if NO variant survives
							end
							if $options.params
								paramline += rseq_var.header + ":" + be4.to_s + "-" + after.to_s + "\t" + snp.snp.to_s + "\t" + parameters[1]
							end
						end
					end
				end
			end
		end
	end
	return [baitsout, outfilter, paramline, coordline, filtercoordline, filteredsnps]
end
#-----------------------------------------------------------------------------------------------
def get_command_line # Get command line for summary output
	# Generate basic command line
	cmdline = $options.algorithm + " -i " + $options.infile
	case $options.algorithm
	when "vcf2baits", "stacks2baits"
		if $options.every
			cmdline += " -e"
			if $options.tiling
				cmdline += " -u -L " + $options.baitlength.to_s + " -O " + $options.tileoffset.to_s + " -k " + $options.tiledepth.to_s
			else
				cmdline += " -b " + $options.lenbef.to_s + " -a " + $options.lenaft.to_s
			end
		else
			cmdline += " -t " + $options.totalsnps.to_s + " -d " + $options.distance.to_s
			if $options.scale
				cmdline += " -j"
			else
				cmdline += " -m " + $options.maxsnps.to_s
			end
			if $options.no_baits
				cmdline += " -p"
			else
				if $options.tiling
					cmdline += " -u -L " + $options.baitlength.to_s + " -O " + $options.tileoffset.to_s + " -k " + $options.tiledepth.to_s
				else
					cmdline += " -b " + $options.lenbef.to_s + " -a " + $options.lenaft.to_s
				end
			end
		end
		cmdline += " -r " + $options.refseq unless $options.no_baits
		cmdline += " -R" if $options.alt_alleles
		cmdline += " -V " + $options.varqual.to_s if $options.varqual_filter
		cmdline += " -S" if $options.sort
		cmdline += " -H -A " + $options.alpha.to_s if $options.hwe
	else
		unless $options.algorithm == "tilebaits" or $options.algorithm == "checkbaits" or $options.algorithm == "aln2baits"
			cmdline += " -r " + $options.refseq
		end
		cmdline += " -L " + $options.baitlength.to_s
		cmdline += " -O " + $options.tileoffset.to_s unless $options.algorithm == "checkbaits"
		if $options.algorithm == "aln2baits"
			cmdline += " -H " + $options.haplodef
		elsif $options.algorithm == "annot2baits"
			cmdline += " -U "
			for feature in $options.features
				cmdline += feature + ","
			end
			cmdline = cmdline[0...-1] # Remove final , from feature list
		end		
	end
	cmdline += " -B" if $options.coords
	cmdline +=	" -D" if $options.ncbi		
	# Generate filtration options
	fltline = ""
	fltline += " -w" if $options.params
	fltline += " -c" if $options.completebait
	fltline += " -G" if $options.no_gaps
	fltline += " -N" if $options.no_Ns
	fltline += " -n " + $options.mingc.to_s if $options.mingc_filter
	fltline += " -x " + $options.maxgc.to_s if $options.maxgc_filter
	fltline += " -q " + $options.mint.to_s if $options.mint_filter
	fltline += " -z " + $options.maxt.to_s if $options.maxt_filter
	if $options.mint_filter or $options.maxt_filter
		fltline += " -T " + $options.bait_type + " -s " + $options.na.to_s + " -f " + $options.formamide.to_s
	end
	fltline += " -Q " + $options.meanqual.to_s if $options.meanqual_filter
	fltline += " -M " + $options.minqual.to_s if $options.minqual_filter
	fltline += " -F " + $options.fasta_score.to_s if ($options.meanqual_filter or $options.minqual_filter)
	return [cmdline,fltline]
end
#-----------------------------------------------------------------------------------------------
class Parser
	def self.parse(options)
		# Set defaults
		args = OpenStruct.new
		algorithms = ["aln2baits", "annot2baits", "bed2baits", "checkbaits", "stacks2baits", "tilebaits", "vcf2baits"]
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
		args.no_gaps = false # Flag to omit bait sequences with gaps
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
					opts.on("-t", "--totalvars [VALUE]", Integer, "Total requested variants (Default = 30,000)") do |tsnps|
						args.totalsnps = tsnps if tsnps != nil
					end
				when "stacks2baits" #stacks2baits only options
					opts.on("-i","--input [FILE]", String, "Input Stacks summary tsv") do |fa|
						args.infile = fa
					end
					opts.on("-S", "--sort", "Sort variants according to variation between/within populations") do
						args.sort = true
					end
					opts.on("-H", "--hwe", "Sort variants within populations according to Hardy-Weinberg Equilibrium (Implies -S)") do
						args.sort = true
						args.hwe = true
					end
					opts.on("-A", "--alpha [VALUE]", Float, "Set alpha value for HWE test (Either 0.10, 0.05, 0.025, 0.01) (Default = 0.05)") do |alph|
						args.alpha = alph if alph != nil
					end
					opts.on("-t", "--totalvars [VALUE]", Integer, "Total requested variants per category (Default = 30,000)") do |tsnps|
						args.totalsnps = tsnps if tsnps != nil
					end
				else
					case args.algorithm
					when "bed2baits"
						inputstr = "Input BED file"
					when "annot2baits"
						inputstr = "Input GFF/GTF table"
					else 
						inputstr = "Input FASTA/FASTQ file"
					end
					opts.on("-i","--input [FILE]", String, inputstr) do |fa|
						args.infile = fa
					end
					if args.algorithm == "bed2baits" or args.algorithm == "annot2baits"
						opts.on("-r","--refseq [FILE]", String, "Reference FASTA/FASTQ sequence file") do |ref|
							args.refseq = ref
						end
					end
					opts.on("-L", "--length [VALUE]", Integer, "Requested bait length (Default = 120)") do |prblength|
						args.baitlength = prblength if prblength != nil
					end
					unless args.algorithm == "checkbaits"
						opts.on("-O", "--offset [VALUE]", Integer, "Base pair offset between tiled baits (Default = 20)") do |toff|
							args.tileoffset = toff if toff != nil
						end
					end
					if args.algorithm == "annot2baits"
						opts.on("-U", "--features [FEATURE]", String, "Comma-separated list of feature types") do |feat|
							args.features = feat.upcase.split(",")
						end
					elsif args.algorithm == "aln2baits"
						opts.on("-H","--haplo [VALUE]", String, "Window haplotype definition (haplotype or variant)") do |fa|
							args.haplodef = fa
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
					 opts.on("-e","--every", "Output bait sequences for every variant in the input file (Overrides -t, -j, -d, -m , -p)") do
						args.every = true
					end
					 opts.on("-R", "--alt", "Generate baits for alternate alleles (Overrides -p)") do
 						args.alt_alleles = true
 					end
					opts.on("-r","--refseq [FILE]", String, "If outputting baits, the reference FASTA/FASTQ sequence file") do |ref|
						args.refseq= ref
					end
					opts.on("-L", "--length [VALUE]", Integer, "Requested total bait length (Default = 120)") do |prblength|
						args.baitlength = prblength if prblength != nil
						args.lenbef =	prblength - 1 #Sets these values to the appropriate position
						args.lenaft =  0 #Sets these values to the appropriate position
					end
					opts.on("-b", "--lenbef [VALUE]", Integer, "If outputting baits, the number of bases before the variant (Default = 60) (Overrides -L)") do |b4|
						args.lenbef = b4 if b4 != nil
					end
					opts.on("-a", "--lenaft [VALUE]", Integer, "If outputting baits, the number of bases after the variant (Default = 59) (Overrides -L)") do |af|
						args.lenaft = af if af != nil
					end
					opts.on("-u", "--tiling", "If outputting baits, tile bait sequences") do
						args.tiling = true
					end
					opts.on("-O", "--offset [VALUE]", Integer, "If tiling baits, base pair offset between tiled baits (Default = 20) (Overrides -a and -b)") do |toff|
						args.tileoffset = toff if toff != nil
					end
					opts.on("-k", "--depth [VALUE]", Integer, "If tiling baits, requested baits per variant") do |tdep|
						args.tiledepth = tdep if tdep != nil
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
				opts.on("-G", "--nogaps", "Exclude bait sequences with gaps (-)") do
					args.no_gaps = true
				end	
				opts.on("-N", "--noNs", "Exclude bait sequences with Ns") do
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
					opts.on("-B", "--bed", "Output BED file for the candidate baits") do
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
					end
					exit
				end
			else
				opts.banner = "Common options:"
				print "Welcome to baitstools " + BAITSTOOLSVER + "\n\nTo use the interactive interface, enter <ruby baitstools.rb [subcommand]> without command-line options.\nCommand-line usage: ruby baitstools.rb [subcommand] [options]"
				print "\nAdd '-h' or '--help' to subcommands (without other options) to see their relevant options.\n\nAvailable subcommands:\n\n"
				print "    aln2baits\t\t\t\tGenerate weighted baits from a FASTA/FASTQ alignment\n"
				print "    annot2baits\t\t\t\tGenerate tiled baits from a GFF/GTF file and a reference sequence\n"				
				print "    bed2baits\t\t\t\tGenerate tiled baits from BED file and a reference sequence\n"
				print "    checkbaits\t\t\t\tFilter a FASTA/FASTQ of candidate baits by quality\n"
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
				print "Sort variants according to variation between/within populations? (y/n)\n"
				t = gets.chomp.upcase
				if t == "Y" or t == "YES"
					$options.sort = true
					print "Sort variants within populations according to Hardy-Weinberg Equilibrium? (y/n)\n"
					t = gets.chomp.upcase
					if t == "Y" or t == "YES"
						$options.hwe = true
						print "Alpha value for HWE (alpha = 0.10, 0.05, 0.025, 0.01)?\n"
						$options.alpha = gets.chomp.to_f
					end
				end
			end
		end
		if $options.algorithm == "stacks2baits"
			alphas = [0.10, 0.05, 0.025, 0.01]
			while !alphas.include?($options.alpha)
				print "Alpha value must be 0.10, 0.05, 0.025, 0.01. Re-enter.\n"
				$options.alpha = gets.chomp.to_f
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
		if $options.algorithm == "bed2baits" or $options.algorithm == "annot2baits"
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
				print "Haplotype definition? (haplotype or variant)\n"
				$options.haplodef = gets.chomp
			end
			while $options.haplodef != 'haplotype' and $options.haplodef != 'variant'
				print "Please choose a haplotype definition (haplotype or variant)\n"
				$options.haplodef = gets.chomp
			end
		end	
	end
	if $options.interact and !$options.no_baits
		print "Do FASTA/FASTQ sequence headers include NCBI-style descriptions? (y/n)\n"
		t = gets.chomp.upcase
		if t == "Y" or t == "YES"
			$options.ncbi = true
		end
		if $options.algorithm != "checkbaits"
			print "Output BED file(s) for candidate baits? (y/n)\n"
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
		print "Exclude baits with gaps (-)? (y/n)\n"
		t = gets.chomp.upcase
		if t == "Y" or t == "YES"
			$options.no_gaps = true
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
	$options.filter = true if ($options.completebait or $options.params or $options.algorithm == "checkbaits" or $options.mingc_filter or $options.maxgc_filter or $options.mint_filter or $options.maxt_filter or $options.meanqual_filter or $options.minqual_filter or $options.no_gaps or $options.no_Ns) # Force filtration as necessary
	cmdline = get_command_line
	print "** Starting program with the following options: **\n"
	print "** Basic command: " + cmdline[0] + " **\n"
	print "** Filtration options:" + cmdline[1] + " **\n" # filtered line always starts with a space if present
	case $options.algorithm
	when "aln2baits"
		aln2baits
	when "annot2baits"
		annot2baits
	when "bed2baits"
		bed2baits
	when "checkbaits"
		checkbaits
	when "stacks2baits"
		stacks2baits
	when "tilebaits"
		tilebaits($options.infile)
	when "vcf2baits"
		vcf2baits
	end
	print "** Program complete **\n"
end
