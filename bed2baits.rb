#!/usr/bin/env ruby
#-----------------------------------------------------------------------------------------------
# bed2baits
BED2BAITSVER = "1.6.7"
# Michael G. Campana, 2017-2020
# Smithsonian Conservation Biology Institute
#-----------------------------------------------------------------------------------------------

def bed2baits
	#Import reference sequence
	print "** Reading reference sequence **\n"
	refseq = read_fasta($options.refseq)
	refhash = {} # sort sequences by name
	for ref in refseq
		refhash[ref.header]=ref
	end
	#Read coordinates table
	if $options.list_format == "bed"
		print "** Reading BED file **\n" 
	else
		print "** Reading interval list **\n"
	end
	regions = [] #Array to hold generated fasta sequences
	write_file(".log.txt", "ExtractedRegions\nRegion\tChromosome\tStart\tEnd\tLength") if $options.log
	totallength = 0
	bed_header = true
	gz_file_open($options.infile).open($options.infile) do |coord|
		while line = coord.gets
			case $options.list_format
			when "bed"
				line_arr = line.split
				bed_header = false if line_arr[0] != "track" && line_arr[0] != "browser" && bed_header # Last comparison to attempt to prevent redefining bed_header with each feature
				if bed_header
					next
				else
					chromo = line_arr[0]
					seqst = line_arr[1].to_i
					seqend = line_arr[2].to_i - 1
				end
			when "gatk"
				chromo = line.split(":")[0]
				seqst = line.split(":")[1].split("-")[0].to_i
				seqend = line.split(":")[1].split("-")[1].to_i - 1
			when "picard"
				if line[0].chr == "@"
					next
				else
					line_arr = line.split("\t")
					chromo = line_arr[0]
					seqst = line_arr[1].to_i - 1
					seqend = line_arr[2].to_i - 1
				end
			end
			seqst -= $options.pad
			seqend += $options.pad
			if refhash.include?(chromo)
				seqst, seqend, seqcycles = get_looped_sequence(refhash, chromo, seqst, seqend)
				seq = get_padded_faseq(refhash, chromo, seqst, seqend, seqcycles)
				regions.push(seq)
				if $options.log
					write_file(".log.txt", seq.header + "\t" + chromo + "\t" + (seqst+1).to_s + "\t" + (seqend+1).to_s + "\t" + seq.seq.length.to_s)
					totallength += seq.seq.length
				end
			else
				print "** Chromosome " + chromo + " not found in reference sequence file. **\n"
			end
		end
	end
	tile_regions(regions, totallength)
end
