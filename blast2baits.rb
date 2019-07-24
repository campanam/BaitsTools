#!/usr/bin/env ruby
#-----------------------------------------------------------------------------------------------
# aln2baits
BLAST2BAITSVER = "1.4.0"
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
	File.open($options.infile, 'r') do |blast|
		while line = blast.gets
			if line[0].chr != "#"
				line_arr = line.split("\t")
				target = line_arr[0]
				chromo = line_arr[1]
				percid = line_arr[2].to_f
				evalue = line_arr[10].to_f
				seqcycles = -2 # Number of times around reference sequence
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
					if refhash[chromo].circular
						while seqend > refhash[chromo].seq.length - 1
							seqcycles += 1
							seqend -= refhash[chromo].seq.length # Correct for padding going off end
						end
						while seqst < 0
							seqcycles += 1
							seqst += refhash[chromo].seq.length
						end				
					else
						if seqst < 0 # Correct for padding going off end
							print "** Sequence " + chromo + " starting coordinate set to 1. **\n"
							seqst = 0 
						end
						if seqend > refhash[chromo].seq.length - 1 # Correct for padding going off end
							print "** Sequence " + chromo + " final coordinate set to " + refhash[chromo].seq.length.to_s + " **\n"
							seqend = refhash[chromo].seq.length - 1
						end
					end
					nuclen = (line_arr[9].to_i - line_arr[8].to_i).abs + 1 # Get length of matching nucleotide in reference
					if nuclen >= $options.blastlen && percid >= $options.percid
						unless $options.evalue_filter && evalue > $options.evalue
							if refhash[chromo].fasta
								seq = Fa_Seq.new(target + "_" + chromo + "_" + (seqst+1).to_s + "-" + (seqend+1).to_s, false, true) #Correct for 0-based counting
							else
								seq = Fa_Seq.new(target + "_" + chromo + "_" + (seqst+1).to_s + "-" + (seqend+1).to_s, false, false) #Correct for 0-based counting
								if seqcycles == -2
									seq.qual = refhash[chromo].qual[seqst..seqend]
								else
									seq.qual = refhash[chromo].qual[seqst..-1]
									(seqcycles+1).times { seq.qual << refhash[chromo].qual }
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
							seq.seq = reversecomp(seq.seq) if line_arr[9].to_i < line_arr[8].to_i	
							seq.bedheader = chromo
							seq.bedstart = seqst 
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
