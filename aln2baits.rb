#!/usr/bin/ruby
#-----------------------------------------------------------------------------------------------
# aln2baits 0.2
# Michael G. Campana, 2016
# Smithsonian Conservation Biology Institute
#-----------------------------------------------------------------------------------------------

class Hap_Window # Object defining a haplotype window
	attr_accessor :header, :seqstart, :seqend, :haplotypes
	def initialize(header, seqstart, seqend, haplotypes = [])
		@header = header
		@seqstart = seqstart
		@seqend = seqend
		@haplotypes = haplotypes
	end
	def var_permutations # Get number of possible variant permutations
		n_vars = 1 # This code ignores wobble bases completely
		for i in 0...@haplotypes[0].length
			varA = 0 #A allele exists?
			varT = 0 #T allele exists?
			varG = 0 #G allele exits?
			varC = 0 #C allele exists?
			indel = 0 #Indel allele exists?
			for seq in @haplotypes
				case seq[i].upcase
				when "A"
					varA = 1
				when "T"
					varT = 1
				when "G"
					varG = 1
				when "C"
					varC = 1
				when "-"
					indel = 1
				end
				break if varA + varT + varG +varC +indel == 5 #End loop early if all vars present
			end
			n_vars *= (varA+varT+varG+varC+indel)
		end
		return n_vars
	end
end
#-----------------------------------------------------------------------------------------------
def aln2baits
	aln = read_fasta($options.infile)
	# Get haplotypes
	seqstart = 0
	windows = []
	while seqstart < aln[0].seq.length
		seqend = seqstart + $options.baitlength - 1
		seqend = aln[0].seq.length - 1 if seqend > aln[0].seq.length-1 # Correct for circular sequences later
		window = Hap_Window.new([], seqstart, seqend)
		for seq in aln
			tmp_seq = seq.seq[seqstart..seqend]
			unless tmp_seq.include?("-") and $options.no_indels # Exclude sequences with indels if not permitted
				unless window.haplotypes.include?(tmp_seq) # Add new sequences to haplotype list
					window.haplotypes.push(tmp_seq) 
					window.header.push(seq.header)
				end
			end
		end
		windows.push(window)
		seqstart += $options.tileoffset # Control the window by tiling density
	end
	# Get bait candidates
	baitsout = ""
	outfilter = ""
	paramline = "Sequence:Coordinates:Haplotype\tBaitLength\t%GC\tTm\tKept\n"
	coordline = ""
	filtercoordline = ""
	for window in windows
		case $options.haplodef
		when "haplotype"
			for hapno in 1..window.haplotypes.size
				rng = (window.seqstart+1).to_s+"-"+(window.seqend+1).to_s #Adjust for 1-based indexing
				baitsout += ">" + window.header[hapno-1] + "_" + rng + "_haplotype" + hapno.to_s + "\n" + window.haplotypes[hapno-1] + "\n"
				coordline += window.header[hapno-1] + ":" + rng + "\n"
				if $options.filter
					flt = filter_baits(window.haplotypes[hapno-1])
					if flt[0]
						outfilter += ">" + window.header[hapno-1] + "_" + rng + "_haplotype" + hapno.to_s+ "\n" + window.haplotypes[hapno-1] + "\n"
						filtercoordline += window.header[hapno-1] + ":" + rng + "\n"
					end
					if $options.params
						paramline += window.header[hapno-1] + ":" + rng + ":" + hapno.to_s + "\t" + flt[1]
					end
				end	
			end
		when "variants"
			nvars = window.var_permutations # Get maximum number of sequences
			for hapno in 1..nvars
				chosen = rand(window.haplotypes.size) # Choose a random haplotype
				haplo = window.haplotypes[chosen] 
				rng = (window.seqstart+1).to_s+"-"+(window.seqend+1).to_s #Adjust for 1-based indexing
				baitsout += ">" + window.header[chosen] + "_" + rng + "_haplotype" + hapno.to_s + "\n" + haplo + "\n"
				coordline += window.header[chosen] + ":" + rng + "\n"
				if $options.filter
					flt = filter_baits(haplo)
					if flt[0]
						outfilter += ">" + window.header[chosen] + "_" + rng + "_haplotype" + hapno.to_s+ "\n" + haplo + "\n"
						filtercoordline += window.header[chosen] + ":" + rng + "\n"
					end
					if $options.params
						paramline += window.header[chosen] + ":" + rng + ":" + hapno.to_s + "\t" + flt[1]
					end
				end	
			end
		end
	end
	write_baits(baitsout, outfilter, paramline, coordline, filtercoordline)
end