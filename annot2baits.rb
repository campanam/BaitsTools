#!/usr/bin/env ruby
#-----------------------------------------------------------------------------------------------
# annot2baits
ANNOT2BAITSVER = "1.7.3"
# Michael G. Campana, 2017-2021
# Smithsonian Conservation Biology Institute
#-----------------------------------------------------------------------------------------------

def annot2baits
	#Import reference sequence
	print "** Reading reference sequence **\n"
	refseq = read_fasta($options.refseq)
	refhash = {} # sort sequences by name
	for ref in refseq
		refhash[ref.header]=ref
	end
	# Read annotation file
	print "** Reading annotation file **\n"
	regions = [] #Array to hold generated fasta sequences
	write_file(".log.txt", "ExtractedRegions\nRegion\tChromosome\tStart\tEnd\tLength") if $options.log
	totallength = 0
	gz_file_open($options.infile) do |annot|
		while line = annot.gets
			if line[0].chr != "#"
				line_arr = line.split("\t")
				chromo = line_arr[0]
				feature = line_arr[2]
				seqst = line_arr[3].to_i - 1 - $options.pad # Correct for 0-based counting and padding
				seqend = line_arr[4].to_i - 1 + $options.pad # Correct for 0-based counting and padding
				unless feature.nil? # Skip bad feature lines such as line breaks
					if $options.features.include?(feature.upcase) # test whether included feature
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
			end
		end
	end
	tile_regions(regions, totallength)
end

