#!/usr/bin/ruby
#-----------------------------------------------------------------------------------------------
# coords2baits 0.2
# Michael G. Campana, 2017
# Smithsonian Conservation Biology Institute
#-----------------------------------------------------------------------------------------------

def coords2baits
	#Import reference sequence
	print "** Reading reference sequence **\n"
	refseq = read_fasta($options.refseq)
	refhash = {} # sort sequences by name
	for ref in refseq
		refhash[ref.header]=ref
	end
	#Read coordinates table
	print "** Reading coordinates table **\n"
	regions = [] #Array to hold generated fasta sequences
	File.open($options.infile, 'r') do |coord|
		while line = coord.gets
			chromosome = ""
			seqst = ""
			seqend = ""
			brk = 0
			for i in 0..line.length-2 # Last character is line break
				case brk
				when 0
					if line[i].chr != ":"
						chromosome += line[i].chr
					else
						brk += 1
					end
				when 1
					if line[i].chr != "-"
						seqst += line[i].chr
					else
						brk += 1
					end
				when 2
					seqend += line[i].chr
				end
			end
			if refhash.include?(chromosome)
				if seqst.to_i < 1
					print "** Chromosome " + chromosome + " starting coordinate set to 1. **\n"
					seqst = "1"
				end
				if seqend.to_i > refhash[chromosome].seq.length 
					print "** Chromosome " + chromosome + " final coordinate set to " + refhash[chromosome].seq.length.to_s + " **\n"
					seqend = refhash[chromosome].seq.length.to_s
				end
				if refhash[chromosome].fasta
					seq = Fa_Seq.new(chromosome + "_" + seqst + "-" + seqend, false, true)
				else
					seq = Fa_Seq.new(chromosome + "_" + seqst + "-" + seqend, false, false)
					seq.qual = refhash[chromosome].qual[seqst.to_i-1..seqend.to_i-1] #Correct for 0-based counting
					seq.calc_quality
				end
				seq.seq = refhash[chromosome].seq[seqst.to_i-1..seqend.to_i-1] #Correct for 0-based counting
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
