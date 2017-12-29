#!/usr/bin/ruby
#-----------------------------------------------------------------------------------------------
# tilebaits 0.3
# Michael G. Campana, 2016
# Smithsonian Conservation Biology Institute
#-----------------------------------------------------------------------------------------------

def tilebaits(seq_array)
	# Process FASTA file
	baitsout = ""
	outfilter = ""
	paramline = "Chromosome:Coordinates\tBaitLength\t%GC\tTm\tKept\n"
	coordline = ""
	filtercoordline = ""
	seq_array = read_fasta(seq_array) if seq_array.is_a?(String) #Read FASTA file if external file
	for seq in seq_array
		sequence = seq.seq
		seqst = 1 # Beginning of bait coordinate
		while seqst < sequence.length+1 # Stop the loop once end of sequence reached
			seqend = seqst+$options.baitlength-1 #End of bait coordinate
			if seqend > sequence.length and !seq.circular #correct for running off end of linear sequence
				seqend = sequence.length
				prb = sequence[seqst-1..seqend-1] #Correct for 0-based counting
				rng = seqst.to_s + "-" + seqend.to_s
			elsif seqend > sequence.length and seq.circular #add beginning of sequence to end
				seqend -= sequence.length
				prb = sequence[seqst-1..sequence.length-1] + sequence[0..seqend-1] #Correct for 0-based counting
				rng = seqst.to_s + "-" + seqend.to_s
			else
				prb = sequence[seqst-1..seqend-1] #Correct for 0-based counting
				rng = seqst.to_s + "-" + seqend.to_s
			end
			baitsout += ">" + seq.header + "_" + rng + "\n" + prb + "\n"
			coordline += seq.header + ":" + rng + "\n"
			if $options.filter
				flt = filter_baits(prb)
				if flt[0]
					outfilter += ">" + seq.header + "_" + rng + "\n" + prb + "\n"
					filtercoordline += seq.header + ":" + rng + "\n"
				end
				if $options.params
					paramline += seq.header + ":" + rng + "\t" + flt[1]
				end
			end
			seqst += $options.tileoffset
		end
	end
	write_baits(baitsout, outfilter, paramline, coordline, filtercoordline)
end
