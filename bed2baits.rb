#!/usr/bin/ruby
#-----------------------------------------------------------------------------------------------
# bed2baits
BED2BAITSVER = "0.3"
# Michael G. Campana, 2017
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
	File.open($options.infile, 'r') do |coord|
		while line = coord.gets
			line_arr = line.split("\t")
			chromosome = line_arr[0]
			seqst = line_arr[1].to_i
			seqend = line_arr[2].to_i
			if refhash.include?(chromosome)
				if seqst < 1
					print "** Chromosome " + chromosome + " starting coordinate set to 1. **\n"
					seqst = 1
				end
				if seqend > refhash[chromosome].seq.length 
					print "** Chromosome " + chromosome + " final coordinate set to " + refhash[chromosome].seq.length.to_s + " **\n"
					seqend = refhash[chromosome].seq.length
				end
				if refhash[chromosome].fasta
					seq = Fa_Seq.new(chromosome + "_" + seqst.to_s + "-" + seqend.to_s, false, true)
				else
					seq = Fa_Seq.new(chromosome + "_" + seqst.to_s + "-" + seqend.to_s, false, false)
					seq.qual = refhash[chromosome].qual[seqst-1..seqend-1] #Correct for 0-based counting
					seq.calc_quality
				end
				seq.seq = refhash[chromosome].seq[seqst-1..seqend-1] #Correct for 0-based counting
				regions.push(seq)
			else
				print "** Chromosome " + chromosome + " not found in reference sequence file. **\n"
			end
		end
	end
	#Write fasta sequences from the files
	outfasta = ""
	for reg in regions
		outfasta += ">" + reg.header + "\n" + reg.seq + "\n"
	end
	File.open($options.infile+"-regions.fa", 'w') do |out|
		out.puts outfasta
	end
	#Generate probes using methods from tilebaits
	tilebaits(regions)
end
