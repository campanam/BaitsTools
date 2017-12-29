#!/usr/bin/ruby
#-----------------------------------------------------------------------------------------------
# baitslib 0.3
# Michael G. Campana, 2016
# Smithsonian Conservation Biology Institute
#-----------------------------------------------------------------------------------------------

# To be implemented
# Sequence complexity filter
# Minimum PIC/Calculate PIC
# Self-complementarity filter

#-----------------------------------------------------------------------------------------------
class Fa_Seq #container for fasta sequences
	attr_accessor :header, :circular, :seq
	def initialize(header = "", circular = false, seq = "")
		@header = header
		@circular = circular
		@seq = seq
	end 
end
#-----------------------------------------------------------------------------------------------
class Chromo_SNP # Container for SNPs
	attr_accessor :chromo, :snp, :popvar_data # Chromosome, SNP index in 1-based indexing, array containing variant data
	def initialize(chromo, snp, popvar_data =[])
		@chromo = chromo
		@snp = snp
		@popvar_data = popvar_data
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
end
#-----------------------------------------------------------------------------------------------
def filter_baits(bait)
	keep = true
	keep = false if bait.length < $options.baitlength
	gc = 0.0
	for i in 0 ... bait.length
		if bait[i].chr.upcase == "G" or bait[i].chr.upcase == "C"
			gc += 1.0
		end
	end
	gccont = gc/bait.length.to_f
	melt = 79.8 + (58.4 * gccont) + (11.8 * (gccont**2.0)) - 820.0/bait.length.to_f + 18.5 * Math::log($options.na)
	if $options.gc
		keep = false if gccont * 100.0 < $options.mingc
		keep = false if gccont * 100.0 > $options.maxgc
	end
	if $options.melt
		keep = false if melt < $options.mint
		keep = false if melt > $options.maxt
	end
	return [keep, bait.length.to_s + "\t" + gccont.to_s + "\t" + melt.to_s + "\t" + keep.to_s + "\n"]
end
#-----------------------------------------------------------------------------------------------
def write_baits(baitsout = "", outfilter = "", paramline = "", coordline = "", filtercoordline = "", filestem = $options.infile)
	File.open(filestem + "-baits.fa", 'w') do |write|
		write.puts baitsout
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
def read_fasta(file)
	seq_array = []
	faseq = nil # Dummy value
	File.open(file, 'r') do |seq|
		while line = seq.gets
			if line[0].chr == ">"
				seq_array.push(faseq) if !faseq.nil? # push previously completed sequence into array
				header = line[1...-1] # Remove final line break and beginning >
				if header[-5..-1] == "#circ" #Adding this to the end of a sequence denotes as circular
					circular = true
					header = header[0...-5] # Remove circ hashtag
				else
					circular = false
				end
				faseq = Fa_Seq.new(header, circular)
			else
				faseq.seq += line[0...-1] #Exclude final line break, allow multi-line fasta
			end
		end
	end
	seq_array.push(faseq) # Put last sequence into fasta array
	return seq_array
end
#-----------------------------------------------------------------------------------------------
def selectsnps(snp_hash) # Choose SNPs based on input group of SNPS
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
				if tmp[i] < selected_snp && selected_snp - tmp[i] < $options.distance
					temp_snps[selected_contig].delete(tmp[i])
				elsif tmp[i] > selected_snp && tmp[i] - selected_snp < $options.distance
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
def snp_to_baits(selectedsnps, filter_search_keys = [])
	baitsout = ""
	outfilter = ""
	paramline = "Chromosome:Coordinates\tSNP\tBaitLength\t%GC\tTm\tKept\n"
	coordline = ""
	filtercoordline = ""
	refseq = read_fasta($options.refseq)
	for rseq in refseq
		if selectedsnps.keys.include?(rseq.header)
			for snp in selectedsnps[rseq.header]
				if $options.tiling
					tiling = $options.tiledepth
				else
					tiling = 1
					$options.tileoffset = 1
				end
				for tile in 0...tiling
					be4 = snp - $options.lenbef + tile * $options.tileoffset
					after = snp + $options.lenaft + tile * $options.tileoffset
					be4 = 1 if (be4 < 1 and !rseq.circular)
					if after > rseq.seq.length and !rseq.circular
						after = rseq.seq.length
						if be4 < 1
							prb = rseq.seq[be4-1..-1] + rseq.seq[0..after-1]  #Correct for 0-based counting
						else
							prb = rseq.seq[be4-1..after-1]  #Correct for 0-based counting
						end 
					elsif after > rseq.seq.length and rseq.circular
						after -= rseq.seq.length
						if be4 < 1
							prb = rseq.seq[be4-1..-1] + rseq.seq[be4-1..rseq.seq.length-1]+rseq.seq[0..after-1]
						else
							prb = rseq.seq[be4-1..rseq.seq.length-1]+rseq.seq[0..after-1]  #Correct for 0-based counting and end of sequence
						end
					else
						if be4 < 1
							prb = rseq.seq[be4-1..-1] + rseq.seq[0..after-1]
						else
							prb = rseq.seq[be4-1..after-1]  #Correct for 0-based counting
						end
					end
					seq = ">" + rseq.header + "\t" + snp.to_s + "\n" + prb + "\n"
					baitsout += seq
					be4 = rseq.seq.length + be4 if be4 < 1
					coord = rseq.header + ":" + be4.to_s + "-" + after.to_s + "\n"
					coordline += coord
					if $options.filter
						parameters = filter_baits(prb) #Correct for 0-based counting and end of sequence
						if parameters[0]
							outfilter += seq
							filtercoordline += coord 
						else
							search = rseq.header + "\t" + snp.to_s + "\t"
							filter_search_keys.delete(search)
						end
						if $options.params
							paramline += rseq.header + ":" + be4.to_s + "-" + after.to_s + "\t" + snp.to_s + "\t" + parameters[1]
						end
					end
				end
			end
		end
	end
	return [baitsout, outfilter, paramline, coordline, filtercoordline, filter_search_keys]
end