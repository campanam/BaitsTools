#!/usr/bin/ruby
#-----------------------------------------------------------------------------------------------
# aln2baits 0.1
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
		window = Hap_Window.new(aln[0].header, seqstart, seqend)
		for seq in aln
			tmp_seq = seq.seq[seqstart..seqend]
			unless tmp_seq.include?("-") and $options.no_indels # Exclude sequences with indels if not permitted
				window.haplotypes.push(tmp_seq) unless window.haplotypes.include?(tmp_seq) # Add new sequences to haplotype list
			end
		end
		windows.push(window)
		seqstart += $options.tileoffset # Control the window by tiling density
	end
	# Get bait candidates
	baitsout = ""
	outfilter = ""
	paramline = "AlignmentReference:Coordinates:Haplotype\tBaitLength\t%GC\tTm\tKept\n"
	coordline = ""
	filtercoordline = ""
	for window in windows
		case $options.haplodef
		when "haplotype"
			hapno = 1
			for haplo in window.haplotypes
				rng = (window.seqstart+1).to_s+"-"+(window.seqend+1).to_s #Adjust for 1-based indexing
				baitsout += ">" + window.header + "_" + rng + "_haplotype" + hapno.to_s + "\n" + haplo + "\n"
				coordline += window.header + ":" + rng + "\n"
				if $options.filter
					flt = filter_baits(haplo)
					if flt[0]
						outfilter += ">" + window.header + "_" + rng + "_haplotype" + hapno.to_s+ "\n" + haplo + "\n"
						filtercoordline += window.header + ":" + rng + "\n"
					end
					if $options.params
						paramline += window.header + ":" + rng + ":" + hapno.to_s + "\t" + "\t" + flt[1]
					end
				end	
				hapno += 1
			end
		when "variants"
			nvars = window.var_permutations # Get maximum number of sequences
			for hapno in 1..nvars
				haplo = window.haplotypes[rand(window.haplotypes.size)] # Choose a random haplotype
				rng = (window.seqstart+1).to_s+"-"+(window.seqend+1).to_s #Adjust for 1-based indexing
				baitsout += ">" + window.header + "_" + rng + "_haplotype" + hapno.to_s + "\n" + haplo + "\n"
				coordline += window.header + ":" + rng + "\n"
				if $options.filter
					flt = filter_baits(haplo)
					if flt[0]
						outfilter += ">" + window.header + "_" + rng + "_haplotype" + hapno.to_s+ "\n" + haplo + "\n"
						filtercoordline += window.header + ":" + rng + "\n"
					end
					if $options.params
						paramline += window.header + ":" + rng + ":" + hapno.to_s + "\t" + "\t" + flt[1]
					end
				end	
			end
		end
	end
	write_baits(baitsout, outfilter, paramline, coordline, filtercoordline)
end