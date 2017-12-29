#!/usr/bin/ruby
#-----------------------------------------------------------------------------------------------
# aln2baits
ALN2BAITSVER = "0.5"
# Michael G. Campana, 2017
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
	def var_permutations # Get possible variant permutations
		varindex = 1 # Index for sequence
		variants = []
		for i in 0...@haplotypes[0].length
			variants.push([])
		end
		for i in 0...@haplotypes[0].length
			for seq in @haplotypes
				if !variants[i].include?(seq[i])
					case seq[i]
					when "A","T","G","C","-"
						variants[i].push(seq[i])
					when "R"
						variants[i].push("A","G")
					when "Y"
						variants[i].push("C","T")
					when "M"
						variants[i].push("A","C")
					when "K"
						variants[i].push("G","T")
					when "S"
						variants[i].push("C","G")
					when "W"
						variants[i].push("A","T")
					when "H"
						variants[i].push("A","C","T")
					when "B"
						variants[i].push("C","G","T")
					when "V"
						variants[i].push("A","C","G")
					when "D"
						variants[i].push("A","G","T")
					when "N"
						variants[i].push("A","C","G","T")
					end
				end
			end
			variants[i].uniq! # Remove duplicate variants due to wobble base handling
			varindex *= variants[i].size
		end
		revised_haplos = []
		for i in 1..varindex
			revised_haplos.push("")
		end
		for i in 0...@haplotypes[0].length
			for k in 0...varindex
				var = k % variants[i].size
				revised_haplos[k]+=variants[i][var]
			end
		end
		return revised_haplos
	end
end
#-----------------------------------------------------------------------------------------------
def aln2baits
	print "** Reading alignment **\n"
	aln = read_fasta($options.infile)
	print "** Identifying haplotypes **\n"
	# Get haplotypes
	seqstart = 0
	windows = []
	while seqstart < aln[0].seq.length
		seqend = seqstart + $options.baitlength - 1
		seqend = aln[0].seq.length - 1 if seqend > aln[0].seq.length-1 # Correct for circular sequences later
		window = Hap_Window.new([], seqstart, seqend)
		for seq in aln
			tmp_seq = seq.seq[seqstart..seqend]
			unless window.haplotypes.include?(tmp_seq) # Add new sequences to haplotype list
				window.haplotypes.push(tmp_seq) 
				window.header.push(seq.header)
			end
		end
		windows.push(window)
		seqstart += $options.tileoffset # Control the window by tiling density
	end
	# Get bait candidates
	print "** Generating and filtering baits **\n"
	baitsout = ""
	outfilter = ""
	paramline = "Sequence:Coordinates:Haplotype\tBaitLength\t%GC\tTm\tMeanQuality\tMinQuality\tKept\n"
	coordline = ""
	filtercoordline = ""
	for window in windows
		case $options.haplodef
		when "haplotype"
			for hapno in 1..window.haplotypes.size
				rng = (window.seqstart+1).to_s+"-"+(window.seqend+1).to_s # Adjust for 1-based indexing
				window.haplotypes[hapno-1].gsub!("T","U") if $options.rna # Will correct both raw and filtered sequences
				baitsout += ">" + window.header[hapno-1] + "_" + rng + "_haplotype" + hapno.to_s + "\n" + window.haplotypes[hapno-1] + "\n"
				coordline += window.header[hapno-1] + "\t" + (window.seqstart+1).to_s + "\t" + (window.seqend+1).to_s + "\n"
				if $options.filter
					flt = filter_baits(window.haplotypes[hapno-1]) # U won't affect filtration
					if flt[0]
						outfilter += ">" + window.header[hapno-1] + "_" + rng + "_haplotype" + hapno.to_s+ "\n" + window.haplotypes[hapno-1] + "\n"
						filtercoordline += window.header[hapno-1] + "\t" + (window.seqstart+1).to_s + "\t" + (window.seqend+1).to_s + "\n"
					end
					if $options.params
						paramline += window.header[hapno-1] + ":" + rng + ":" + hapno.to_s + "\t" + flt[1]
					end
				end	
			end
		when "variant"
			nvars = window.var_permutations # Get window_permutations
			for hapno in 1..nvars.size
				haplo = nvars[hapno-1]
				haplo.gsub!("T","U") if $options.rna # Will correct both raw and filtered sequences
				rng = (window.seqstart+1).to_s + "-" + (window.seqend+1).to_s #Adjust for 1-based indexing
				baitsout += ">Alignment_" + rng + "_haplotype" + hapno.to_s + "\n" + haplo + "\n" # Original window headers are meaningless
				coordline += "Alignment\t" + (window.seqstart+1).to_s + "\t" + (window.seqend+1).to_s + "\n" # Original window headers are meaningless
				if $options.filter
					flt = filter_baits(haplo) # U won't affect filtration
					if flt[0]
						outfilter += ">Alignment_" + rng + "_haplotype" + hapno.to_s+ "\n" + haplo + "\n"
						filtercoordline += "Alignment\t" + (window.seqstart+1).to_s + "\t" + (window.seqend+1).to_s + "\n"
					end
					if $options.params
						paramline += "Alignment:" + rng + ":" + hapno.to_s + "\t" + flt[1]
					end
				end	
			end
		end
	end
	write_baits(baitsout, outfilter, paramline, coordline, filtercoordline)
end
