#!/usr/bin/env ruby
#-----------------------------------------------------------------------------------------------
# annot2baits
ANNOT2BAITSVER = "1.1.0"
# Michael G. Campana, 2017-2018
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
				seqst = 0 if seqst < 0 # Correct for padding going off end
				seqend = line_arr[4].to_i - 1 + $options.pad # Correct for 0-based counting and padding
				unless feature.nil? # Skip bad feature lines such as line breaks
					if $options.features.include?(feature.upcase) # test whether included feature
						seqend = refhash[chromo].seq.length - 1 if seqend > refhash[chromo].seq.length - 1 # Correct for padding going off end
						if refhash[chromo].fasta
							seq = Fa_Seq.new(chromo + "_" + (seqst+1).to_s + "-" + (seqend+1).to_s, false, true) #Correct for 0-based counting
						else
							seq = Fa_Seq.new(chromo + "_" + (seqst+1).to_s + "-" + (seqend+1).to_s, false, false) #Correct for 0-based counting
							seq.qual = refhash[chromo].qual[seqst..seqend] 
							seq.calc_quality
						end
						seq.seq = refhash[chromo].seq[seqst..seqend]
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
	#Write fasta sequences from the files
	for reg in regions
		write_file("-regions.fa", ">" + reg.header + "\n" + reg.seq)
	end
	write_file(".log.txt", "\nTotalRegions\tTotalRegionLength\n" + regions.size.to_s + "\t" + totallength.to_s + "\n") if $options.log
	#Generate probes using methods from tilebaits
	tilebaits(regions)
end

