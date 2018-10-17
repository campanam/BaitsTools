#!/usr/bin/env ruby
#-----------------------------------------------------------------------------------------------
# bed2baits
BED2BAITSVER = "1.2.3"
# Michael G. Campana, 2017-2018
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
	print "** Reading BED file **\n"
	regions = [] #Array to hold generated fasta sequences
	write_file(".log.txt", "ExtractedRegions\nRegion\tChromosome\tStart\tEnd\tLength") if $options.log
	totallength = 0
	File.open($options.infile, 'r') do |coord|
		while line = coord.gets
			line_arr = line.split("\t")
			chromosome = line_arr[0]
			seqst = line_arr[1].to_i - $options.pad
			seqend = line_arr[2].to_i + $options.pad
			if refhash.include?(chromosome)
				if seqst < 0
					print "** Chromosome " + chromosome + " starting coordinate set to 1. **\n"
					seqst = 0
				end
				if seqend > refhash[chromosome].seq.length 
					print "** Chromosome " + chromosome + " final coordinate set to " + refhash[chromosome].seq.length.to_s + " **\n"
					seqend = refhash[chromosome].seq.length
				end
				if refhash[chromosome].fasta
					seq = Fa_Seq.new(chromosome + "_" + (seqst+1).to_s + "-" + seqend.to_s, false, true)
				else
					seq = Fa_Seq.new(chromosome + "_" + (seqst+1).to_s + "-" + seqend.to_s, false, false)
					seq.qual = refhash[chromosome].qual[seqst..seqend-1] #Correct for 0/1-based counting
					seq.calc_quality
				end
				seq.seq = refhash[chromosome].seq[seqst..seqend-1] #Correct for 0-based counting
				seq.bedheader = chromosome
				seq.bedstart = seqst
				regions.push(seq)
				if $options.log
					write_file(".log.txt", seq.header + "\t" + chromosome + "\t" + (seqst+1).to_s + "\t" + seqend.to_s + "\t" + seq.seq.length.to_s)
					totallength += seq.seq.length
				end
			else
				print "** Chromosome " + chromosome + " not found in reference sequence file. **\n"
			end
		end
	end
	if regions.size == 0 # Break out if no regions found
		print "** No annotated regions found. Exiting.\n **"
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
