#!/usr/bin/env ruby
#-----------------------------------------------------------------------------------------------
# annot2baits
ANNOT2BAITSVER = "1.4.0"
# Michael G. Campana, 2017-2019
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
	File.open($options.infile, 'r') do |annot|
		while line = annot.gets
			if line[0].chr != "#"
				line_arr = line.split("\t")
				chromo = line_arr[0]
				feature = line_arr[2]
				seqst = line_arr[3].to_i - 1 - $options.pad # Correct for 0-based counting and padding
				seqend = line_arr[4].to_i - 1 + $options.pad # Correct for 0-based counting and padding
				seqcycles = -2 # Number of times around reference sequence
				unless feature.nil? # Skip bad feature lines such as line breaks
					if $options.features.include?(feature.upcase) # test whether included feature
						if refhash[chromo].circular
							while seqend > refhash[chromo].seq.length - 1
								seqcycles += 1
								seqend -= refhash[chromo].seq.length if seqend > refhash[chromo].seq.length - 1 # Correct for padding going off end
							end
							while seqst < 0
								seqcycles += 1
								seqst += refhash[chromo].seq.length
							end				
						else
							if seqend > refhash[chromo].seq.length - 1 # Correct for padding going off end
								print "** Chromosome " + chromo + " starting coordinate set to 1. **\n"
								seqend = refhash[chromo].seq.length - 1
							end
							if seqst < 0
								print "** Chromosome " + chromo + " final coordinate set to " + refhash[chromo].seq.length.to_s + " **\n"
								seqst = 0
							end
						end
						if refhash[chromo].fasta
							seq = Fa_Seq.new(chromo + "_" + (seqst+1).to_s + "-" + (seqend+1).to_s, false, true) #Correct for 0-based counting
						else
							seq = Fa_Seq.new(chromo + "_" + (seqst+1).to_s + "-" + (seqend+1).to_s, false, false) #Correct for 0-based counting
							if seqcycles == -2
								seq.qual = refhash[chromo].qual[seqst..seqend]
							else
								seq.qual = refhash[chromo].qual[seqst..-1]
								seqcycles.times { seq.qual << refhash[chromo].qual }
								seq.qual << refhash[chromo].qual[0..seqend]
							end
							seq.calc_quality
						end
						if seqcycles == -2
							seq.seq = refhash[chromo].seq[seqst..seqend]
						else
							seq.seq = refhash[chromo].seq[seqst..-1]
							(seqcycles+1).times { seq.seq << refhash[chromo].seq }
							seq.seq << refhash[chromo].seq[0..seqend]
						end
						seq.bedheader = chromo
						seq.bedstart = seqst 
						regions.push(seq)
						if $options.log
							write_file(".log.txt", seq.header + "\t" + chromo + "\t" + (seqst+1).to_s + "\t" + (seqend+1).to_s + "\t" + seq.seq.length.to_s)
							totallength += seq.seq.length
						end
					end
				end
			end
		end
	end
	tile_regions(regions, totallength)
end

