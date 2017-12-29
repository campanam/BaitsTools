#!/usr/bin/env ruby
#-----------------------------------------------------------------------------------------------
# baitstools
BAITSTOOLSVER = "0.9"
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
require_relative 'pyrad2baits'

#-----------------------------------------------------------------------------------------------
class Fa_Seq #container for fasta/fastq sequences
	attr_accessor :header, :circular, :fasta, :seq, :qual, :qual_array, :bedstart, :locus
	def initialize(header = "", circular = false, fasta = false, seq = "", qual = "", bedstart = 0, locus = "")
		@header = header # Sequence header
		@circular = circular # Circular sequence flag
		@seq = seq # DNA sequence
		@qual = qual # Original quality scores
		@qual_array = [] # Array of numeric quality scores
		@fasta = fasta # FASTA format flag
		@bedstart = bedstart # Sequence start in absolute BED coordinates
		@locus = locus # locus id for alignment
	end
	def make_dna # Replace uracils with thymines for internal consistency
		@seq.gsub!("u", "t")
		@seq.gsub!("U", "T")
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
def get_ambiguity(base)
	vars = []
	case base
	when "A","a","T","t","G","g","C","c","-"
		vars.push(base)
	when "R"
		vars.push("A","G")
	when "Y"
		vars.push("C","T")
	when "M"
		vars.push("A","C")
	when "K"
		vars.push("G","T")
	when "S"
		vars.push("C","G")
	when "W"
		vars.push("A","T")
	when "H"
		vars.push("A","C","T")
	when "B"
		vars.push("C","G","T")
	when "V"
		vars.push("A","C","G")
	when "D"
		vars.push("A","G","T")
	when "N"
		vars.push("A","C","G","T")
	when "r"
		vars.push("a","g")
	when "y"
		vars.push("c","t")
	when "m"
		vars.push("a","c")
	when "k"
		vars.push("g","t")
	when "s"
		vars.push("c","g")
	when "w"
		vars.push("a","t")
	when "h"
		vars.push("a","c","t")
	when "b"
		vars.push("c","g","t")
	when "v"
		vars.push("a","c","g")
	when "d"
		vars.push("a","g","t")
	when "n"
		vars.push("a","c","g","t")
	end
	return vars
end
#-----------------------------------------------------------------------------------------------
def collapse_ambiguity(bait) # Ambiguity handling
	r = ["A","G"]
	rl = ["a","g"]
	y = ["C","T"]
	yl = ["c","t"]
	m = ["A","C"]
	ml = ["a","c"]
	k = ["G","T"]
	kl = ["g","t"]
	s = ["C","G"]
	sl = ["c","g"]
	w = ["A","T"]
	wl = ["a","t"]
	h = ["A","C","T"]
	hl = ["a","c","t"]
	b = ["C","G","T"]
	bl = ["c","g","t"]
	v = ["A","C","G"]
	vl = ["a","c","g"]
	d = ["A","G","T"]
	dl = ["a","g","t"]
	n = ["A","C","G","T"]
	nl = ["a","c","g","t"]
	bait.gsub!("R",r[rand(2)])
	bait.gsub!("Y",y[rand(2)])
	bait.gsub!("M",m[rand(2)])
	bait.gsub!("K",k[rand(2)])
	bait.gsub!("S",s[rand(2)])
	bait.gsub!("W",w[rand(2)])
	bait.gsub!("H",h[rand(3)])
	bait.gsub!("B",b[rand(3)])
	bait.gsub!("V",v[rand(3)])
	bait.gsub!("D",d[rand(3)])
	bait.gsub!("N",n[rand(4)]) unless $options.no_Ns
	bait.gsub!("T","U") if $options.rna # Remove any residual Ts after ambiguity collapsing
	bait.gsub!("r",rl[rand(2)])
	bait.gsub!("y",yl[rand(2)])
	bait.gsub!("m",ml[rand(2)])
	bait.gsub!("k",kl[rand(2)])
	bait.gsub!("s",sl[rand(2)])
	bait.gsub!("w",wl[rand(2)])
	bait.gsub!("h",hl[rand(3)])
	bait.gsub!("b",bl[rand(3)])
	bait.gsub!("v",vl[rand(3)])
	bait.gsub!("d",dl[rand(3)])
	bait.gsub!("n",nl[rand(4)]) unless $options.no_Ns
	bait.gsub!("t","u") if $options.rna # Remove any residual Ts after ambiguity collapsing
	return bait
end
#-----------------------------------------------------------------------------------------------
def extend_baits(bait, reference, seqst, seqend) # Extend baits with gap characters
	bait.delete!("-")
	while bait.length < $options.baitlength
		seqend += 1
		break if seqend >= reference.length
		bait += reference[seqend].chr
		bait.delete!("-")
	end
	while bait.length < $options.baitlength
		seqst -=1
		break if seqst < 0
		bait.insert(0,reference[seqst].chr)
		bait.delete!("-")
	end
	return bait
end
#-----------------------------------------------------------------------------------------------
def filter_baits(bait, qual = [0])

	# To be implemented:
	# Self-complementarity filter
	
	bait = collapse_ambiguity(bait) if $options.collapse_ambiguities # Ambiguity handling
	keep = true
	keep = false if (bait.length < $options.baitlength && $options.completebait)
	keep = false if (bait.include?("N") and $options.no_Ns)
	keep = false if (bait.include?("-") and $options.gaps == "exclude")
	gc = 0.0
	maskbases = ["a","t","g","c","u","r","y","m","k","s","w","h","b","v","d","n"] # Array of masked base possibilities
	mask = 0.0
	for i in 0 ... bait.length
		if bait[i].chr.upcase == "G" or bait[i].chr.upcase == "C"
			gc += 1.0
		end
		mask += 1.0 if maskbases.include?(bait[i].chr)
	end
	gccont = gc/bait.length.to_f
	maskcont = mask/bait.length.to_f
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
	keep = false if (maskcont * 100.0 > $options.maxmask && $options.maxmask_filter)
	keep = false if (meanqual < $options.meanqual && $options.meanqual_filter)
	keep = false if (minqual < $options.minqual && $options.minqual_filter)
	return [keep, bait.length.to_s + "\t" + (gccont * 100.0).to_s + "\t" + melt.to_s + "\t" + maskcont.to_s + "\t" +  meanqual.to_s + "\t" + minqual.to_s + "\t" + bait.include?("N").to_s + "\t" + bait.include?("-").to_s + "\t" + keep.to_s + "\n"]
end
#-----------------------------------------------------------------------------------------------
def write_baits(baitsout = [""], outfilter = [""], paramline = [""], coordline = [""], filtercoordline = [""], filestem = $options.outdir + "/" + $options.outprefix, rbedline = [""], rfilterline = [""])
	[baitsout,outfilter,paramline,coordline,filtercoordline, rbedline, rfilterline].each do |output|
		output.flatten!
		output.delete(nil)
		output.join("\n")
	end
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
	if $options.rbed
		File.open(filestem + "-baits-relative.bed", 'w') do |write|
			write.puts rbedline
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
		if $options.rbed
			File.open(filestem + "-filtered-baits-relative.bed", 'w') do |write|
				write.puts rfilterline
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
						seq_array.push(faseq)
					end
					header = line[1...-1] # Remove final line break and beginning >
					header_array = header.split("#")
					if header_array.include?("circ")  # Determine circularity
						circular = true
						header_array.delete("circ") # Delete array value so bed value always at same index
					else
						circular = false
					end
					# Get locus information
					locus = header_array.find { |gene| /^loc/ =~ gene }
					header_array.delete(locus) unless locus.nil?
					if header_array[1].nil?
						bedstart = 0 # Set default bedstart to 0
					else
						bedstart = header_array[1][3..-1].to_i
					end
					header_array = header_array[0].split(" ") if $options.ncbi # Remove sequence description from sequence name for NCBI-formatted files
					header = header_array[0]
					faseq = Fa_Seq.new(header, circular, true)
					faseq.bedstart = bedstart
					faseq.locus = locus[3..-1] unless locus.nil?
				elsif qual and faseq.qual.length < faseq.seq.length # Test if quality line is complete
					faseq.qual += line[0...-1] #Exclude final line break, allow multi-line fastq
			 	elsif line[0].chr == "@" # Start new sequence
					qual = false
					if !faseq.nil? # push previously completed sequence into array
						seq_array.push(faseq)
					end
					header = line[1...-1] # Remove final line break and beginning >
					header_array = header.split("#")
					if header_array.include?("circ")  # Determine circularity
						circular = true
						header_array.delete!("circ") # Delete array value so bed value always at same index
					else
						circular = false
					end
					# Get locus information
					locus = header_array.find { |gene| /^loc/ =~ gene }
					header_array.delete(locus) unless locus.nil?
					if header_array[1].nil?
						bedstart = 0 # Set default bedstart to 0
					else
						bedstart = header_array[1][3..-1].to_i
					end
					header_array = header_array[0].split(" ") if $options.ncbi # Remove sequence description from sequence name for NCBI-formatted files
					header = header_array[0]
					faseq = Fa_Seq.new(header, circular, false)
					faseq.bedstart = bedstart
					faseq.locus = locus[3..-1] unless locus.nil?
				elsif line[0].chr == "+"
					qual = true
				else
					faseq.seq += line[0...-1] #Exclude final line break, allow multi-line fasta/fastq
				end
			end
		end
	end
	seq_array.push(faseq) # Put last sequence into fasta array
	threads = [] # Array to hold threads
	$options.threads.times do |i|
		threads[i] = Thread.new {
			for Thread.current[:j] in 0 ... seq_array.size
				if Thread.current[:j] % $options.threads == i
					seq_array[Thread.current[:j]].calc_quality
					seq_array[Thread.current[:j]].make_dna
				end
			end
		}
	end
	threads.each { |thr| thr.join }
	if $options.log
		$options.logtext += "SequencesRead\nType\tNumberCircular\tNumberLinear\tTotalNumber\tCircularBp\tLinearBp\tTotalBp\n"
		facirc = 0
		falin = 0
		facircbp = 0
		falinbp = 0
		fqcirc = 0
		fqlin = 0
		fqcircbp = 0
		fqlinbp = 0
		for seq in seq_array
			if seq.fasta && seq.circular
				facirc += 1
				facircbp += seq.seq.length
			elsif seq.fasta && !seq.circular
				falin += 1
				falinbp += seq.seq.length
			elsif !seq.fasta && seq.circular
				fqcirc += 1
				fqcircbp += seq.seq.length
			else
				fqlin += 1
				fqlinbp += seq.seq.length
			end
		end
		$options.logtext += "FASTA\t" + facirc.to_s + "\t" + falin.to_s + "\t" + (facirc + falin).to_s + "\t" + facircbp.to_s + "\t" + falinbp.to_s + "\t" + (facircbp + falinbp).to_s + "\n"
		$options.logtext += "FASTQ\t" + fqcirc.to_s + "\t" + fqlin.to_s + "\t" + (fqcirc + fqlin).to_s + "\t" + fqcircbp.to_s + "\t" + fqlinbp.to_s + "\t" + (fqcircbp + fqlinbp).to_s + "\n"
		$options.logtext += "Total\t" + (facirc + fqcirc).to_s + "\t" + (falin + fqlin).to_s + "\t" + (facirc + falin + fqcirc + fqlin).to_s + "\t" + (facircbp + fqcircbp).to_s + "\t" + (falinbp + fqlinbp).to_s + "\t" + (facircbp + falinbp + fqcircbp + fqlinbp).to_s + "\n\n"
	end
	return seq_array
end
#-----------------------------------------------------------------------------------------------
def selectsnps(snp_hash) # Choose SNPs based on input group of SNPSs
	# Sort chromosomal SNPs in case unsorted
	for chromo in snp_hash.keys
		snp_hash[chromo].sort_by! { |snp| snp.snp }
	end
	temp_snps = snp_hash.dup #Avoid messing with the original hash.
	selectsnps = {}
	if !$options.every
		for i in 1..$options.totalsnps
			selected_contig = temp_snps.keys[rand(temp_snps.size)] # Get name of contig
			snpindex = rand(temp_snps[selected_contig].size)
			selected_snp = temp_snps[selected_contig][snpindex]
			# Delete SNPs that are too close
			if snpindex + 1 < temp_snps[selected_contig].size
				while (temp_snps[selected_contig][snpindex + 1].snp - selected_snp.snp).abs < $options.distance
					temp_snps[selected_contig].delete(temp_snps[selected_contig][snpindex + 1])
					break if snpindex + 1 >= temp_snps[selected_contig].size
				end
			end
			if snpindex > 0
				while (temp_snps[selected_contig][snpindex - 1].snp - selected_snp.snp).abs < $options.distance
					temp_snps[selected_contig].delete(temp_snps[selected_contig][snpindex - 1])
					snpindex -= 1
					break if snpindex < 1
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
			break if temp_snps.size == 0 # Stop selecting SNPs when all contigs deleted from consideration
		end
	else
		selectsnps = snp_hash
	end
	if $options.log
		totalvar = 0
		selectvar = 0
		for chromo in snp_hash.keys
			totalvar += snp_hash[chromo].size
		end
		$options.logtext += "Chromosome\tSelectedVariants\n"
	end
	for chromo in selectsnps.keys
		selectsnps[chromo].sort_by! { |snp| snp.snp }
		if $options.log
			selectvar += selectsnps[chromo].size
			$options.logtext += chromo + "\t"
			for snp in selectsnps[chromo]
				$options.logtext += snp.snp.to_s + ","
			end
			$options.logtext = $options.logtext[0...-1] # Remove final comma
		end
	end
	$options.logtext += "\n\nNumberTotalVariants\tNumberSelectedVariants\n" if $options.log
	$options.logtext += totalvar.to_s + "\t" + selectvar.to_s + "\n\n" if $options.log
	return selectsnps
end
#-----------------------------------------------------------------------------------------------
def snp_to_baits(selectedsnps, refseq)
	baitsout = []
	outfilter = []
	paramline = ["Chromosome:Coordinates\tSNP\tBaitLength\tGC%\tTm\tMasked%\tMeanQuality\tMinQuality\tNs\tGaps\tKept\n"]
	coordline = []
	filtercoordline = []
	refseq.size.times do
		baitsout.push([])
		outfilter.push([])
		paramline.push([])
		coordline.push([])
		filtercoordline.push([])
	end
	filteredsnps = {}
	threads = []
	if $options.log
		logs = []
		refseq.size.times { logs.push([]) }
		$options.logtext += "Chromosome\tVariant\tNumberBaits\tRetainedBaits\tExcludedBaits\n"
	end
	$options.threads.times do |i|
		threads[i] = Thread.new {
			for Thread.current[:j] in 0...refseq.size
				if Thread.current[:j] % $options.threads == i
					if selectedsnps.keys.include?(refseq[Thread.current[:j]].header)
						for Thread.current[:snp] in selectedsnps[refseq[Thread.current[:j]].header]
							Thread.current[:totalbaits] = 0 if $options.log
							Thread.current[:retainedbaits] = 0 if $options.log
							Thread.current[:tiling] = $options.tiledepth
							Thread.current[:rseq_vars] = [refseq[Thread.current[:j]]]
							if $options.alt_alleles
								for Thread.current[:altvar] in Thread.current[:snp].alt_alleles
									Thread.current[:tmpseq] = refseq[Thread.current[:j]].seq.dup # Need to duplicate to prevent changes to original sequence
									if Thread.current[:snp].ref.nil?
										Thread.current[:tmpseq][Thread.current[:snp].snp-1]=Thread.current[:altvar]
									else
										Thread.current[:tmpseq][Thread.current[:snp].snp-1..Thread.current[:snp].snp+Thread.current[:snp].ref.length-2]=Thread.current[:altvar]
									end
									Thread.current[:tmpheader] = refseq[Thread.current[:j]].header.dup + "_" + Thread.current[:altvar] + "-variant"
									Thread.current[:tmp] = Fa_Seq.new(Thread.current[:tmpheader], refseq[Thread.current[:j]].circular, refseq[Thread.current[:j]].fasta, Thread.current[:tmpseq], refseq[Thread.current[:j]].qual)
									unless refseq[Thread.current[:j]].fasta
										Thread.current[:qual_arr] = []
										Thread.current[:tmpqual] = refseq[Thread.current[:j]].qual_array.dup
										if $options.algorithm == "vcf2baits"
											Thread.current[:qual] = Thread.current[:snp].qual
											Thread.current[:qual] = 40 if Thread.current[:qual] > 40 # Adjust for quality scores beyond typical sequence capability
										else
											Thread.current[:qual] = refseq[Thread.current[:j]].numeric_quality[Thread.current[:snp].snp-1] # Assume quality is equal across the variants since no quality information
										end
										for Thread.current[:i] in 1..Thread.current[:altvar].length # Apply quality score to all bases in alternate allele
											Thread.current[:qual_arr].push(Thread.current[:qual])
										end
										Thread.current[:tmpqual][Thread.current[:snp].snp-1]=Thread.current[:qual_arr]
										Thread.current[:tmp].qual_array = Thread.current[:tmpqual].flatten
									end
									Thread.current[:rseq_vars].push(Thread.current[:tmp]) # Add alternate variants to the sequence
								end
							end
							for Thread.current[:rseq_var] in Thread.current[:rseq_vars]
								for Thread.current[:tile] in 0...Thread.current[:tiling]
									Thread.current[:be4] = Thread.current[:snp].snp - $options.lenbef + Thread.current[:tile] * $options.tileoffset
									Thread.current[:after] = Thread.current[:snp].snp + $options.baitlength + 1 - $options.lenbef + Thread.current[:tile] * $options.tileoffset
									Thread.current[:be4] = 1 if (Thread.current[:be4] < 1 and !Thread.current[:rseq_var].circular)
									Thread.current[:qual] = [$options.fasta_score]
									if Thread.current[:after] > Thread.current[:rseq_var].seq.length and !Thread.current[:rseq_var].circular
										Thread.current[:after] = Thread.current[:rseq_var].seq.length
										if Thread.current[:be4] < 1
											Thread.current[:prb] = Thread.current[:rseq_var].seq[Thread.current[:be4]-1..-1] + Thread.current[:rseq_var].seq[0..Thread.current[:after]-1]  #Correct for 0-based counting
											Thread.current[:qual] = Thread.current[:rseq_var].numeric_quality[Thread.current[:be4]-1..-1] + Thread.current[:rseq_var].numeric_quality[0..Thread.current[:after]-1] unless Thread.current[:rseq_var].fasta
										else
											Thread.current[:prb] = Thread.current[:rseq_var].seq[Thread.current[:be4]-1..Thread.current[:after]-1]  #Correct for 0-based counting
											Thread.current[:qual] = Thread.current[:rseq_var].numeric_quality[Thread.current[:be4]-1..Thread.current[:after]-1] unless Thread.current[:rseq_var].fasta
										end
									elsif Thread.current[:after] > Thread.current[:rseq_var].seq.length and Thread.current[:rseq_var].circular
										Thread.current[:after] -= Thread.current[:rseq_var].seq.length
										if Thread.current[:be4] < 1
											Thread.current[:prb] = Thread.current[:rseq_var].seq[Thread.current[:be4]-1..-1] + Thread.current[:rseq_var].seq[Thread.current[:be4]-1..Thread.current[:rseq_var].seq.length-1]+Thread.current[:rseq_var].seq[0..Thread.current[:after]-1] #Correct for 0-based counting and end of sequence
											Thread.current[:qual] = Thread.current[:rseq_var].numeric_quality[Thread.current[:be4]-1..-1] + Thread.current[:rseq_var].numeric_quality[Thread.current[:be4]-1..Thread.current[:rseq_var].seq.length-1]+Thread.current[:rseq_var].numeric_quality[0..Thread.current[:after]-1] unless Thread.current[:rseq_var].fasta
										else
											Thread.current[:prb] = Thread.current[:rseq_var].seq[Thread.current[:be4]-1..Thread.current[:rseq_var].seq.length-1]+Thread.current[:rseq_var].seq[0..Thread.current[:after]-1]  #Correct for 0-based counting and end of sequence
											Thread.current[:qual] = Thread.current[:rseq_var].numeric_quality[Thread.current[:be4]-1..Thread.current[:rseq_var].seq.length-1]+rseq.numeric_quality[0..Thread.current[:after]-1] unless Thread.current[:rseq_var].fasta
										end
									else
										if Thread.current[:be4] < 1
											Thread.current[:prb] = Thread.current[:rseq_var].seq[Thread.current[:be4]-1..-1] + Thread.current[:rseq_var].seq[0..Thread.current[:after]-1]
											Thread.current[:qual] = Thread.current[:rseq_var].numeric_quality[Thread.current[:be4]-1..-1] + Thread.current[:rseq_var].numeric_quality[0..Thread.current[:after]-1] unless Thread.current[:rseq_var].fasta
										else
											Thread.current[:prb] = Thread.current[:rseq_var].seq[Thread.current[:be4]-1..Thread.current[:after]-1]  #Correct for 0-based counting
											Thread.current[:qual] = Thread.current[:rseq_var].numeric_quality[Thread.current[:be4]-1..Thread.current[:after]-1] unless Thread.current[:rseq_var].fasta
										end
									end
									Thread.current[:prb].gsub!("T","U") if $options.rna # RNA output handling
									Thread.current[:prb].gsub!("t","u") if $options.rna # RNA output handling
									Thread.current[:prb] = extend_baits(Thread.current[:prb], Thread.current[:rseq_var].seq, Thread.current[:be4]-1, Thread.current[:after]-1) if $options.gaps == "extend" # Basic gap extension
									Thread.current[:seq] = ">" + Thread.current[:rseq_var].header + "_site" + Thread.current[:snp].snp.to_s + "\n" + Thread.current[:prb] + "\n"
									Thread.current[:be4] = Thread.current[:rseq_var].seq.length + Thread.current[:be4] if Thread.current[:be4] < 1
									Thread.current[:coord] = Thread.current[:rseq_var].header + "\t" + (Thread.current[:be4]-1).to_s + "\t" + Thread.current[:after].to_s + "\n"
									baitsout[Thread.current[:j]].push(Thread.current[:seq])
									Thread.current[:totalbaits] += 1 if $options.log
									coordline[Thread.current[:j]].push(Thread.current[:coord])
									if $options.filter
										Thread.current[:parameters] = filter_baits(Thread.current[:prb], Thread.current[:qual]) # U should not affect filtration
										if Thread.current[:parameters][0]
											outfilter[Thread.current[:j]].push(Thread.current[:seq])
											filtercoordline[Thread.current[:j]].push(Thread.current[:coord])
											Thread.current[:retainedbaits] += 1 if $options.log
											if filteredsnps[refseq[Thread.current[:j]].header].nil?
												filteredsnps[refseq[Thread.current[:j]].header] = [Thread.current[:snp]]
											else
												filteredsnps[refseq[Thread.current[:j]].header].push(Thread.current[:snp])
												filteredsnps[refseq[Thread.current[:j]].header].uniq!
											end
										end
										if $options.params
											Thread.current[:param] = Thread.current[:rseq_var].header + ":" + Thread.current[:be4].to_s + "-" + Thread.current[:after].to_s + "\t" + Thread.current[:snp].snp.to_s + "\t" + Thread.current[:parameters][1]
											paramline[Thread.current[:j]+1].push(Thread.current[:param])
										end
									end
								end
							end
							if $options.log
								Thread.current[:vals] = [refseq[Thread.current[:j]].header, Thread.current[:snp].snp, Thread.current[:totalbaits]]
								if $options.filter
									Thread.current[:vals].push(Thread.current[:retainedbaits], Thread.current[:totalbaits] - Thread.current[:retainedbaits])
								else
									Thread.current[:vals].push(Thread.current[:totalbaits], "NA")
								end
								logs[Thread.current[:j]].push(Thread.current[:vals])
							end
						end
					end
				end
			end
		}
	end
	threads.each { |thr| thr.join }
	if $options.log
		vlogs = [[],[]]
		for i in 0 ... logs.size
			for l in 0 ... logs[i].size
				$options.logtext += logs[i][l].join("\t") + "\n"
				vlogs[0].push(logs[i][l][2])
				vlogs[1].push(logs[i][l][3])
			end
		end
		$options.logtext += "\nTotalBaitCoverage(√ó)\tFilteredBaitCoverage(√ó)\n"
		if $options.filter
			$options.logtext += mean(vlogs[0]).to_s + "\t" + mean(vlogs[1]).to_s + "\n"
		else
			$options.logtext += mean(vlogs[0]).to_s + "\tNA\n"
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
			cmdline += " -L" + $options.baitlength.to_s + " -O" + $options.tileoffset.to_s + " -b" + $options.lenbef.to_s + " -k" + $options.tiledepth.to_s
		else
			cmdline += " -t" + $options.totalsnps.to_s + " -d" + $options.distance.to_s
			if $options.scale
				cmdline += " -j"
			else
				cmdline += " -m" + $options.maxsnps.to_s
			end
			if $options.no_baits
				cmdline += " -p"
			else
				cmdline += " -L" + $options.baitlength.to_s + " -O" + $options.tileoffset.to_s + " -b" + $options.lenbef.to_s + " -k" + $options.tiledepth.to_s
			end
		end
		cmdline += " -r " + $options.refseq unless $options.no_baits
		cmdline += " -a" if $options.alt_alleles
		cmdline += " -V" + $options.varqual.to_s if $options.varqual_filter
		cmdline += " -S" if $options.sort
		cmdline += " -H -A" + $options.alpha.to_s if $options.hwe
	else
		if $options.algorithm == "annot2baits" or $options.algorithm == "bed2baits"
			cmdline += " -r " + $options.refseq
			cmdline += " -P" + $options.pad.to_s
		end
		cmdline += " -L" + $options.baitlength.to_s
		cmdline += " -O" + $options.tileoffset.to_s unless $options.algorithm == "checkbaits"
		if $options.algorithm == "pyrad2baits"
			cmdline += " -I" + $options.minind.to_s
			cmdline += " -W " + $options.strategy
		end
		if $options.algorithm == "aln2baits" or ($options.algorithm == "pyrad2baits" && $options.strategy == "alignment")
			cmdline += " -H " + $options.haplodef
		elsif $options.algorithm == "annot2baits"
			cmdline += " -U "
			for feature in $options.features
				cmdline += feature + ","
			end
			cmdline = cmdline[0...-1] # Remove final , from feature list
		end
		if $options.algorithm == "pyrad2baits" && $options.strategy != "alignment"
			cmdline += " -t" + $options.totalsnps.to_s + " -m" + $options.maxsnps.to_s + " -d" + $options.distance.to_s + " -k" + $options.tiledepth.to_s
			cmdline += " -a" if $options.alt_alleles
		end
	end
	cmdline += " -o " + $options.outprefix
	cmdline += " -Z " + $options.outdir
	cmdline += " -l" if $options.log
	cmdline += " -B" if $options.coords
	cmdline += " -E" if $options.rbed
	cmdline += " -D" if $options.ncbi
	cmdline += " -Y" if $options.rna
	cmdline += " -G " + $options.gaps
	cmdline += " -X" + $options.threads.to_s
	# Generate filtration options
	fltline = ""
	fltline += " -w" if $options.params
	fltline += " -c" if $options.completebait
	fltline += " -N" if $options.no_Ns
	fltline += " -C" if $options.collapse_ambiguities
	fltline += " -n" + $options.mingc.to_s if $options.mingc_filter
	fltline += " -x" + $options.maxgc.to_s if $options.maxgc_filter
	fltline += " -q" + $options.mint.to_s if $options.mint_filter
	fltline += " -z" + $options.maxt.to_s if $options.maxt_filter
	if $options.mint_filter or $options.maxt_filter
		fltline += " -T " + $options.bait_type + " -s" + $options.na.to_s + " -f" + $options.formamide.to_s
	end
	fltline += " -K" + $options.maxmask.to_s if $options.maxmask_filter
	fltline += " -Q" + $options.meanqual.to_s if $options.meanqual_filter
	fltline += " -M" + $options.minqual.to_s if $options.minqual_filter
	fltline += " -F" + $options.fasta_score.to_s if ($options.meanqual_filter or $options.minqual_filter)
	return [cmdline,fltline]
end
#-----------------------------------------------------------------------------------------------
class Parser
	def self.parse(options)
		# Set defaults
		args = OpenStruct.new
		algorithms = ["aln2baits", "annot2baits", "bed2baits", "checkbaits", "pyrad2baits", "stacks2baits", "tilebaits", "vcf2baits"]
		args.algorithm = ARGV[0] # Define subcommands
		if algorithms.include?(args.algorithm) #Set interactive mode or whether to kill
			args.interact = true if ARGV.size == 1
			args.kill = false
		else
			args.kill = true # kill the program if neither interactive nor command line
		end
		ARGV.delete_at(0) # Remove subcommand from input otherwise crashes
		args.threads = 1 # Number of threads
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
		args.coords = false # Flag to output BED file of baits (absolute coordinates)
		args.rbed = false # Flag to output BED file of baits (relative coordinates)
		args.gaps = "include" # Gapped bait handling
		args.no_Ns = false # Flag to omit bait sequences with Ns
		args.collapse_ambiguities = false # Flag to collapse ambiguities to a single nucleotide
		args.haplodef = "haplotype" # Haplotype definition for aln2baits
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
		args.features = [] # Array holding desired features
		args.pad = 0 # BP to pad ends of extracted regions
		args.ncbi = false # Flag whether FASTA/FASTQ headers include NCBI-style descriptors
		args.rna = false # Flag whether baits are output as RNA
		args.alt_alleles = false # Flag to apply alternate alleles
		args.varqual_filter = false # Flag to determine whether to filter vcf variants by QUAL scores
		args.varqual = 30 # Minimum vcf variant QUAL score
		args.outdir = File.expand_path("./") # Output directory
		args.outprefix = "out" # Output prefix
		args.log = false # Flag to output detailed log
		args.logtext = "Parsed commands: " # Variable to store log information
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
					opts.on("-A", "--alpha [VALUE]", Float, "Set alpha value for HWE test (Either 0.10, 0.05, 0.025, 0.01) (Default = 0.05)") do |alph|
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
						args.strategy = strat if strat != nil
					end
					opts.on("-H","--haplo [VALUE]", String, "If using alignment strategy, window haplotype definition (haplotype or variant) (Default = haplotype)") do |fa|
						args.haplodef = fa if fa != nil
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
						opts.on("-P", "--pad [VALUE]", Integer, "Length in bp to pad beginning and ending of extracted regions (Default = 0)") do |pad|
							args.pad = pad if pad != nil
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
							args.haplodef = fa if fa != nil
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
					opts.on("-k", "--depth [VALUE]", Integer, "Requested baits per variant (Default = 1)") do |tdep|
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
				opts.on("-N", "--noNs", "Exclude bait sequences with Ns") do
					args.no_Ns = true
				end
				opts.on("-C", "--collapse", "Collapse ambiguities to a single nucleotide") do
					args.collapse_ambiguities = true
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
					args.bait_type = hyb if hyb != nil
				end
				opts.on("-s","--na [VALUE]", Float, "Melting temperature sodium concentration in M (Default = 0.9)") do |na|
					args.na = na if na != nil
				end
				opts.on("-f", "--formamide [VALUE]", Float, "Melting temperature formamide concentration in M (Default = 0.0)") do |form|
					args.formamide = form if form != nil
				end
				opts.on("-K", "--maxmask [VALUE]", Float, "Maximum % sequence consists of repetitive/low-complexity elements (Default 25.0)") do |maxrep|
					args.maxmask = maxrep if maxrep != nil
					args.maxmask_filter = true # Force bait filtration if called
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
					if args.algorithm == "aln2baits" or args.algorithm == "annot2baits" or args.algorithm == "bed2baits" or args.algorithm == "tilebaits"
						opts.on("-E", "--rbed", "Output BED file for the baits relative to extracted sequences") do
							args.rbed = true
						end
					end
				end
				unless args.algorithm == "pyrad2baits"
					opts.on("-D", "--ncbi", "FASTA/FASTQ file headers have NCBI-style descriptions") do
						args.ncbi = true
					end
				end
				opts.on("-Y","--rna", "Output baits as RNA rather than DNA") do
					args.rna = true
				end
				opts.on("-G", "--gaps [VALUE]", String, "Strategy to handle sequence gaps (-) (include, exclude, or extend) (Default = include)") do |gap|
					args.gaps = gap if gap != nil
				end
				opts.on("-X", "--threads [VALUE]", Integer, "Number of threads (Default = 1)") do |thr|
					args.threads = thr if thr != nil
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
	if $options.interact
		print "Enter output file prefix.\n"
		$options.outprefix = gets.chomp
		print "Enter output file directory.\n"
		$options.outdir = File.expand_path(gets.chomp)
		print "Output detailed log? (y/n)\n"
		t = gets.chomp.upcase
		if t == "Y" or t == "YES"
			$options.log = true
		end
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
			while $options.tileoffset > $options.baitlength or $options.tileoffset < 1
				print "Tiling offset cannot be less than 1 or greater than bait length. Please re-enter.\n"
				$options.tileoffset = gets.chomp.to_i
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
		if $options.interact
			print "Enter tiling bp offset.\n" unless $options.algorithm == "checkbaits"
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
				$options.strategy = gets.chomp
			end
			while $options.strategy != 'alignment' and $options.strategy != 'SNPs' and $options.strategy != 'informative'
				print "Please choose a strategy (alignment, SNPs, or informative)\n"
				$options.strategy = gets.chomp
			end
			if $options.strategy != "alignment"
				if $options.interact
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
				$options.haplodef = gets.chomp
			end
			while $options.haplodef != 'haplotype' and $options.haplodef != 'variant'
				print "Please choose a haplotype definition (haplotype or variant)\n"
				$options.haplodef = gets.chomp
			end
		end
	end
	if $options.interact and !$options.no_baits
		if $options.algorithm != "pyrad2baits"
			print "Do FASTA/FASTQ sequence headers include NCBI-style descriptions? (y/n)\n"
			t = gets.chomp.upcase
			if t == "Y" or t == "YES"
				$options.ncbi = true
			end
		end
		if $options.algorithm != "checkbaits"
			print "Output BED file(s) for candidate baits? (y/n)\n"
			t = gets.chomp.upcase
			if t == "Y" or t == "YES"
				$options.coords = true
			end
			if $options.algorithm == "aln2baits" or $options.algorithm == "annot2baits" or $options.algorithm == "bed2baits" or $options.algorithm == "tilebaits"
				print "Output BED file(s) for the baits relative to extracted sequences? (y/n)\n"
				t = gets.chomp.upcase
				if t == "Y" or t == "YES"
					$options.rbed = true
				end
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
		$options.gaps = "" # Turn off default to allow value entry entry
	end
	while $options.gaps != 'include' and $options.gaps != 'exclude' and $options.gaps != 'extend'
		print "Strategy to handle sequence gaps (-) ('include', 'exclude', or 'extend')\n"
		$options.gaps = gets.chomp.downcase
	end
	if $options.interact and !$options.no_baits
		print "Exclude baits with Ns? (y/n)\n"
		t = gets.chomp.upcase
		if t == "Y" or t == "YES"
			$options.no_Ns = true
		end
		print "Collapse ambiguity codes to a single nucleotide? (y/n)\n"
		t = gets.chomp.upcase
		if t == "Y" or t == "YES"
			$options.collapse_ambiguities = true
		end
		print "Output baits as RNA rather than DNA? (y/n)\n"
		t = gets.chomp.upcase
		if t == "Y" or t == "YES"
			$options.rna = true
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
			$options.bait_type = gets.chomp
		end
	end
	if $options.interact
		print "Filter by maximum percent masked sequence? (y/n)\n"
		t = gets.chomp.upcase
		if t == "Y" or t == "YES"
			$options.maxmask_filter = true
			print "Enter maximum percent masked sequence\n"
			$options.maxmask = gets.chomp.to_f
		end
	end
	while $options.maxmask < 0.0 or $options.maxmask > 100.0
		print "Percent masked sequence must be between 0.0 and 100.0. Re-enter.\n"
		$options.maxmask = gets.chomp.to_f
	end
	if $options.interact and $options.algorithm != "aln2baits" and !$options.no_baits
		print "Filter by minimum mean base quality? (y/n)\n"
		t = gets.chomp.upcase
		if t == "Y" or t == "YES"
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
		t = gets.chomp.upcase
		if t == "Y" or t == "YES"
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
	$options.no_baits = false if ($options.every or $options.alt_alleles) # Override -p when needed
	$options.filter = true if ($options.completebait or $options.params or $options.algorithm == "checkbaits" or $options.mingc_filter or $options.maxgc_filter or $options.mint_filter or $options.maxt_filter or $options.maxmask_filter or $options.meanqual_filter or $options.minqual_filter or $options.gaps == "exclude" or $options.no_Ns or $options.collapse_ambiguities) # Force filtration as necessary
	cmdline = get_command_line
	print "** Starting program with the following options: **\n"
	print "** Basic command: " + cmdline[0] + " **\n"
	print "** Filtration options:" + cmdline[1] + " **\n" # filtered line always starts with a space if present
	$options.logtext += cmdline[0] + cmdline[1] + "\n\n" if $options.log
	case $options.algorithm
	when "aln2baits"
		aln2baits($options.infile)
	when "annot2baits"
		annot2baits
	when "bed2baits"
		bed2baits
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
	if $options.log
		File.open($options.outdir + "/" + $options.outprefix + ".log.txt", 'w') do |write|
			write.puts $options.logtext
		end
	end
	print "** Program complete **\n"
end
