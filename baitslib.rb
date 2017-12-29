#!/usr/bin/ruby
#-----------------------------------------------------------------------------------------------
# baitslib 0.4
# Michael G. Campana, 2016
# Smithsonian Conservation Biology Institute
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
	attr_accessor :chromo, :snp, :popvar_data, :alt, :qual
	def initialize(chromo, snp, popvar_data =[], alt = nil, qual = nil, line = nil)
		@chromo = chromo # Chromosome
		@snp = snp # SNP index in 1-based indexing
		@popvar_data = popvar_data # Array holding stacks population-specific allele data
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
		if @line.nil? # Combine data from Stacks2baits Popvars if line not previously defined
			@line = ""
			for pop in @popvar_data
				@line += pop.line
			end
		end
		return @line
	end
	def alt_alleles 
		if @alt.nil? # Get alternate alleles for Stacks2Baits (once)
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
	keep = false if (bait.upcase.include?("-") and $options.no_indels)
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
		File.open(filestem + "-baits-coords.txt", 'w') do |write|
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
			File.open(filestem + "-filtered-baits-coords.txt", 'w') do |write|
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
						tmpseq[snp.snp-1]=altvar
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
								prb = rseq_var.seq[be4-1..-1] + rseq_var.seq[be4-1..rseq_var.seq.length-1]+rseq_var.seq[0..after-1]
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
						coord = rseq_var.header + ":" + be4.to_s + "-" + after.to_s + "\n"
						coordline += coord
						if $options.filter
							parameters = filter_baits(prb, qual) #Correct for 0-based counting and end of sequence
							if parameters[0]
								outfilter += seq
								filtercoordline += coord 
							elsif 
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
			cmdline += " -G "
			for feature in $options.features
				cmdline += feature + ","
			end
			cmdline = cmdline[0...-1] # Remove final , from feature list
		end		
	end
	cmdline += " -o" if $options.coords
	cmdline +=	" -D" if $options.ncbi		
	# Generate filtration options
	fltline = ""
	fltline += " -w" if $options.params
	fltline += " -c" if $options.completebait
	fltline += " -I" if $options.no_indels
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
