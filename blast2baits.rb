#!/usr/bin/env ruby
#-----------------------------------------------------------------------------------------------
# aln2baits
BLAST2BAITSVER = "1.2.3"
# Michael G. Campana, 2017-2018
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
				if line_arr[9].to_i > line_arr[8].to_i
					seqst = line_arr[8].to_i - 1 - $options.pad # Correct for 0-based counting and padding
					seqend = line_arr[9].to_i - 1 + $options.pad # Correct for 0-based counting and padding
					strand = "+"
				else 
					seqend = line_arr[8].to_i - 1 + $options.pad # Correct for 0-based counting and padding
					seqst = line_arr[9].to_i - 1 - $options.pad # Correct for padding going off end
					strand = "-"
				end
				seqst = 0 if seqst < 0 # Correct for padding going off end
				seqend = refhash[chromo].seq.length - 1 if seqend > refhash[chromo].seq.length - 1 # Correct for padding going off end
				nuclen = (line_arr[9].to_i - line_arr[8].to_i).abs + 1 # Get length of matching nucleotide in reference
				if nuclen >= $options.blastlen && percid >= $options.percid
					unless $options.evalue_filter && evalue > $options.evalue
						if refhash[chromo].fasta
							seq = Fa_Seq.new(target + "_" + chromo + "_" + (seqst+1).to_s + "-" + (seqend+1).to_s, false, true) #Correct for 0-based counting
						else
							seq = Fa_Seq.new(target + "_" + chromo + "_" + (seqst+1).to_s + "-" + (seqend+1).to_s, false, false) #Correct for 0-based counting
							seq.qual = refhash[chromo].qual[seqst..seqend] 
							seq.calc_quality
						end
						if line_arr[9].to_i > line_arr[8].to_i	
							seq.seq = refhash[chromo].seq[seqst..seqend]
						else
							seq.seq = reversecomp(refhash[chromo].seq[seqst..seqend])
						end
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
	if regions.size == 0 # Exit if no BLAST hit found
		print "** No matching BLAST hit found. Exiting. **\n"
		exit
	else
		#Write fasta sequences from the files
		for reg in regions
			write_file("-regions.fa", ">" + reg.header + "\n" + reg.seq)
		end
		write_file(".log.txt", "\nTotalRegions\tTotalRegionLength\n" + regions.size.to_s + "\t" + totallength.to_s + "\n") if $options.log
		#Generate probes using methods from tilebaits
		tilebaits(regions)
	end
end
