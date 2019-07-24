#!/usr/bin/env ruby
#-----------------------------------------------------------------------------------------------
# bed2baits
BED2BAITSVER = "1.3.0"
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
	if $options.list_format == "bed"
		print "** Reading BED file **\n" 
	else
		print "** Reading interval list **\n"
	end
	regions = [] #Array to hold generated fasta sequences
	write_file(".log.txt", "ExtractedRegions\nRegion\tChromosome\tStart\tEnd\tLength") if $options.log
	totallength = 0
	File.open($options.infile, 'r') do |coord|
		while line = coord.gets
			case $options.list_format
			when "bed"
				line_arr = line.split("\t")
				chromosome = line_arr[0]
				seqst = line_arr[1].to_i
				seqend = line_arr[2].to_i
			when "GATK"
				chromosome = line.split(":")[0]
				seqst = line.split(":")[1].split("-")[0].to_i
				seqend = line.split(":")[1].split("-")[1].to_i
			when "Picard"
				if line[0].chr == "@"
					next
				else
					line_arr = line.split("\t")
					chromosome = line_arr[0]
					seqst = line_arr[1].to_i - 1
					seqend = line_arr[2].to_i
				end
			end
			seqst -= $options.pad
			seqend += $options.pad
			seqcycles = -2 # Number of times around reference sequence
			if refhash.include?(chromosome)
				if refhash[chromosome].circular
					while seqend > refhash[chromosome].seq.length - 2 # Correct for 0/1-based counting
						seqcycles += 1
						seqend -= refhash[chromosome].seq.length # Correct for padding going off end
					end
					while seqst < 0
						seqcycles += 1
						seqst += refhash[chromosome].seq.length
					end	
				else
					if seqst < 0
						print "** Chromosome " + chromosome + " starting coordinate set to 1. **\n"
						seqst = 0
					end
					if seqend > refhash[chromosome].seq.length 
						print "** Chromosome " + chromosome + " final coordinate set to " + refhash[chromosome].seq.length.to_s + " **\n"
						seqend = refhash[chromosome].seq.length
					end
				end
				if refhash[chromosome].fasta
					seq = Fa_Seq.new(chromosome + "_" + (seqst+1).to_s + "-" + seqend.to_s, false, true)
				else
					seq = Fa_Seq.new(chromosome + "_" + (seqst+1).to_s + "-" + seqend.to_s, false, false)
					if seqcycles == -2
						seq.qual = refhash[chromosome].qual[seqst..seqend-1] #Correct for 0/1-based counting
					else
						seq.qual = refhash[chromosome].qual[seqst..-1]
						(seqcycles+1).times { seq.qual << refhash[chromosome].qual }
						seq.qual << refhash[chromosome].qual[0..seqend-1]
					end
					seq.calc_quality
				end
				if seqcycles == -2
					seq.seq = refhash[chromosome].seq[seqst..seqend-1] #Correct for 0-based counting
				else
					seq.seq = refhash[chromosome].seq[seqst..-1]
					(seqcycles+1).times { seq.seq << refhash[chromosome].seq }
					seq.seq << refhash[chromosome].seq[0..seqend-1]
				end 
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
	tile_regions(regions, totallength)
end
