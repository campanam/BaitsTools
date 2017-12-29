#!/usr/bin/env ruby
#-----------------------------------------------------------------------------------------------
# aln2baits
ALN2BAITSVER = "0.6"
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
		variants = []
		for i in 0...@haplotypes[0].length
			variants.push([])
		end
		threads = [] # Array to hold threads
		$options.threads.times do |j|
			threads[j] = Thread.new {
				for Thread.current[:i] in 0...@haplotypes[0].length
					if Thread.current[:i] % $options.threads == j
						Thread.current[:vars] = []
						for Thread.current[:seq] in @haplotypes
							case Thread.current[:seq][Thread.current[:i]]
							when "A","T","G","C","-"
								Thread.current[:vars].push(Thread.current[:seq][Thread.current[:i]])
							when "R"
								Thread.current[:vars].push("A","G")
							when "Y"
								Thread.current[:vars].push("C","T")
							when "M"
								Thread.current[:vars].push("A","C")
							when "K"
								Thread.current[:vars].push("G","T")
							when "S"
								Thread.current[:vars].push("C","G")
							when "W"
								Thread.current[:vars].push("A","T")
							when "H"
								Thread.current[:vars].push("A","C","T")
							when "B"
								Thread.current[:vars].push("C","G","T")
							when "V"
								Thread.current[:vars].push("A","C","G")
							when "D"
								Thread.current[:vars].push("A","G","T")
							when "N"
								Thread.current[:vars].push("A","C","G","T")
							end
						end
						Thread.current[:vars].uniq! # Remove duplicate variants
						variants[Thread.current[:i]] = Thread.current[:vars] # Minimize lock time
					end
				end
			}
		end
		threads.each { |thr| thr.join }
		varindex = 1 # Index for sequence
		for var in variants
			varindex *= var.size # Must be complete down here or interferes with multithreading
		end
		revised_haplos = []
		for i in 1..varindex
			revised_haplos.push("")
		end
		threads = []
		$options.threads.times do |j|
			threads[j] = Thread.new {
				for Thread.current[:k] in 0...varindex
					if Thread.current[:k] % $options.threads == j
						for Thread.current[:i] in 0...@haplotypes[0].length
							Thread.current[:var] = Thread.current[:k] % variants[Thread.current[:i]].size
							revised_haplos[Thread.current[:k]]+=variants[Thread.current[:i]][Thread.current[:var]] # Minimize lock time
						end
					end
				end
			}
		end
		threads.each { |thr| thr.join }
		@haplotypes = revised_haplos
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
		window.var_permutations if $options.haplodef == "variant" # Call here for code efficiency
		windows.push(window)
		seqstart += $options.tileoffset # Control the window by tiling density
	end
	# Get bait candidates
	print "** Generating and filtering baits **\n"
	baitsout = []
	outfilter = []
	paramline = ["Sequence:Coordinates:Haplotype\tBaitLength\t%GC\tTm\tMeanQuality\tMinQuality\tKept\n"]
	coordline = []
	filtercoordline = []
	threads = []
	windows.size.times do
		baitsout.push([])
		outfilter.push([])
		paramline.push([])
		coordline.push([])
		filtercoordline.push([])
	end
	$options.threads.times do |i|
		threads[i] = Thread.new {
			for Thread.current[:j] in 0 ... windows.size
				if Thread.current[:j] % $options.threads == i
					for Thread.current[:hapno] in 1..windows[Thread.current[:j]].haplotypes.size
						Thread.current[:rng] = (windows[Thread.current[:j]].seqstart+1).to_s+"-"+(windows[Thread.current[:j]].seqend+1).to_s # Adjust for 1-based indexing
						windows[Thread.current[:j]].haplotypes[Thread.current[:hapno]-1].gsub!("T","U") if $options.rna # Will correct both raw and filtered sequences
						$options.haplodef == "haplotype" ? Thread.current[:header] = windows[Thread.current[:j]].header[Thread.current[:hapno]-1] : Thread.current[:header] = "Alignment"
						Thread.current[:bait] = ">" + Thread.current[:header] + "_" + Thread.current[:rng] + "_haplotype" + Thread.current[:hapno].to_s + "\n" + windows[Thread.current[:j]].haplotypes[Thread.current[:hapno]-1] + "\n"
						Thread.current[:coord] = Thread.current[:header] + "\t" + (windows[Thread.current[:j]].seqstart).to_s + "\t" + (windows[Thread.current[:j]].seqend+1).to_s + "\n"
						baitsout[Thread.current[:j]].push(Thread.current[:bait])
						coordline[Thread.current[:j]].push(Thread.current[:coord])
						if $options.filter
							Thread.current[:flt] = filter_baits(windows[Thread.current[:j]].haplotypes[Thread.current[:hapno]-1]) # U won't affect filtration
							if Thread.current[:flt][0]
								outfilter[Thread.current[:j]].push(Thread.current[:bait])
								filtercoordline[Thread.current[:j]].push(Thread.current[:coord])
							end
							if $options.params
								Thread.current[:param] = Thread.current[:header] + ":" + Thread.current[:rng] + ":" + Thread.current[:hapno].to_s + "\t" + Thread.current[:flt][1]
								paramline[Thread.current[:j]+1].push(Thread.current[:param])
							end
						end
					end
				end
			end
		}
	end
	threads.each { |thr| thr.join }
	write_baits(baitsout, outfilter, paramline, coordline, filtercoordline)
end
