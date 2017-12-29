#!/usr/bin/ruby
#-----------------------------------------------------------------------------------------------
# checkbaits
CHECKBAITSVER = "0.2"
# Michael G. Campana, 2017
# Smithsonian Conservation Biology Institute
#-----------------------------------------------------------------------------------------------

def checkbaits
	# Process FASTA file
	print "** Reading FASTA/FASTQ **\n"
	outfilter = ""
	paramline = "Bait\tBaitLength\t%GC\tTm\tMeanQuality\tMinQuality\tKept\n"
	seq_array = read_fasta($options.infile)
	print "** Filtering baits **\n"
	for seq in seq_array
		if seq.fasta
			flt = filter_baits(seq.seq, [$options.fasta_score])
		else
			flt = filter_baits(seq.seq, seq.numeric_quality)
		end
		if flt[0]
			outfilter += ">" + seq.header + "\n" + seq.seq + "\n"
		end
		if $options.params
			paramline += seq.header + "\t" + flt[1]
		end
	end
	write_baits("", outfilter, paramline)
end
