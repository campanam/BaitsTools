#!/usr/bin/ruby
#-----------------------------------------------------------------------------------------------
# annot2baits
ANNOT2BAITSVER = "0.2"
# Michael G. Campana, 2017
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
	File.open($options.infile, 'r') do |annot|
		while line = annot.gets
			if line[0].chr != "#"
				line_arr = line.split("\t")
				chromo = line_arr[0]
				feature = line_arr[2]
				seqst = line_arr[3]
				seqend = line_arr[4]
				unless feature.nil? # Skip bad feature lines such as line breaks
					if $options.features.include?(feature.upcase) # test whether included feature
						if refhash[chromo].fasta
							seq = Fa_Seq.new(chromo + "_" + seqst + "-" + seqend, false, true)
						else
							seq = Fa_Seq.new(chromo + "_" + seqst + "-" + seqend, false, false)
							seq.qual = refhash[chromo].qual[seqst.to_i-1..seqend.to_i-1] #Correct for 0-based counting
							seq.calc_quality
						end
						seq.seq = refhash[chromo].seq[seqst.to_i-1..seqend.to_i-1]  #Correct for 0-based counting
						regions.push(seq)
					end
				end
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

