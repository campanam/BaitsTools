#!/usr/bin/ruby
#-----------------------------------------------------------------------------------------------
# tilebaits
TILEBAITSVER = "0.6"
# Michael G. Campana, 2017
# Smithsonian Conservation Biology Institute
#-----------------------------------------------------------------------------------------------


def tilebaits(seq_array)
	# Process FASTA file
	baitsout = ""
	outfilter = ""
	paramline = "Chromosome:Coordinates\tBaitLength\t%GC\tTm\tMeanQuality\tMinQuality\tKept\n"
	coordline = ""
	filtercoordline = ""
	if seq_array.is_a?(String) #Read FASTA file if external file
		print "** Reading FASTA/FASTQ **\n"
		seq_array = read_fasta(seq_array)
	end
	print "** Generating and filtering baits **\n"
	for seq in seq_array
		qual = [$options.fasta_score]
		seqst = 1 # Beginning of bait coordinate
		while seqst < seq.seq.length+1 # Stop the loop once end of sequence reached
			seqend = seqst+$options.baitlength-1 #End of bait coordinate
			if seqend > seq.seq.length and !seq.circular #correct for running off end of linear sequence
				seqend = seq.seq.length
				prb = seq.seq[seqst-1..seqend-1] #Correct for 0-based counting
				qual = seq.numeric_quality[seqst-1..seqend-1] unless seq.fasta
			elsif seqend > seq.seq.length and seq.circular #add beginning of sequence to end
				seqend -= seq.seq.length
				prb = seq.seq[seqst-1..seq.seq.length-1] + seq.seq[0..seqend-1] #Correct for 0-based counting
				qual = seq.numeric_quality[seqst-1..seq.seq.length-1] + seq.numeric_quality[0..seqend-1] unless seq.fasta
			else
				prb = seq.seq[seqst-1..seqend-1] #Correct for 0-based counting
				qual = seq.numeric_quality[seqst-1..seqend-1] unless seq.fasta
			end
			prb.gsub!("T","U") if $options.rna # RNA output handling
			baitsout += ">" + seq.header + "_" + seqst.to_s + "-" + seqend.to_s + "\n" + prb + "\n"
			coordline += seq.header + "\t" + seqst.to_s + "\t" + seqend.to_s + "\n"
			if $options.filter
				flt = filter_baits(prb, qual) # U should not affection filtration
				if flt[0]
					outfilter += ">" + seq.header + "_" + seqst.to_s + "-" + seqend.to_s + "\n" + prb + "\n"
					filtercoordline += seq.header + "\t" + seqst.to_s + "\t" + seqend.to_s + "\n"
				end
				if $options.params
					paramline += seq.header + ":" + seqst.to_s + "-" + seqend.to_s + "\t" + flt[1]
				end
			end
			seqst += $options.tileoffset
		end
	end
	write_baits(baitsout, outfilter, paramline, coordline, filtercoordline)
end
