#!/usr/bin/env ruby
#-----------------------------------------------------------------------------------------------
# blast2baits
BLAST2BAITSVER = "1.6.0"
# Michael G. Campana, 2017-2019
# Smithsonian Conservation Biology Institute
#-----------------------------------------------------------------------------------------------

def blast2baits
	#Import reference sequence
	print "** Reading reference sequence **\n"
	refseq = read_fasta($options.refseq)
	refhash = {} # sort sequences by name
	for ref in refseq
		refhash[ref.header]=ref
	end
	# Read blast file	
	print "** Reading BLAST file **\n"
	regions = [] #Array to hold generated fasta sequences
	write_file(".log.txt", "ExtractedRegions\nRegion\tChromosome\tStart\tEnd\tLength\tStrand") if $options.log
	totallength = 0
	eval(gz_file_open($options.infile)).open($options.infile) do |blast|
		while line = blast.gets
			if line[0].chr != "#"
				line_arr = line.split("\t")
				target = line_arr[0]
				chromo = line_arr[1]
				percid = line_arr[2].to_f
				evalue = line_arr[10].to_f
				if line_arr[9].to_i > line_arr[8].to_i
					seqst = line_arr[8].to_i - 1 - $options.pad # Correct for 0-based counting and padding
					seqend = line_arr[9].to_i - 1 + $options.pad # Correct for 0-based counting and padding
					strand = "+"
				else 
					seqend = line_arr[8].to_i - 1 + $options.pad # Correct for 0-based counting and padding
					seqst = line_arr[9].to_i - 1 - $options.pad # Correct for padding going off end
					strand = "-"
				end
				if refhash[chromo].nil?
					print "** Sequence " + chromo + " not found in reference sequence file. **\n"
				else
					seqst, seqend, seqcycles = get_looped_sequence(refhash, chromo, seqst, seqend)
					nuclen = (line_arr[9].to_i - line_arr[8].to_i).abs + 1 # Get length of matching nucleotide in reference
					if nuclen >= $options.blastlen && percid >= $options.percid
						unless $options.evalue_filter && evalue > $options.evalue
							seq = get_padded_faseq(refhash, chromo, seqst, seqend, seqcycles)
							if line_arr[9].to_i < line_arr[8].to_i
								seq.seq, seq.qual = reversecomp(seq.seq, seq.qual_array)
								seq.qual.reverse!
							end
							regions.push(seq)
							if $options.log
								write_file(".log.txt", seq.header + "\t" + chromo + "\t" + (seqst+1).to_s + "\t" + (seqend+1).to_s + "\t" + seq.seq.length.to_s + "\t" + strand)
								totallength += seq.seq.length
							end
						end
					end
				end
			end
		end
	end
	tile_regions(regions, totallength)
end
