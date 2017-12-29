#!/usr/bin/env ruby
#-----------------------------------------------------------------------------------------------
# tilebaits
TILEBAITSVER = "0.7"
# Michael G. Campana, 2017
# Smithsonian Conservation Biology Institute
#-----------------------------------------------------------------------------------------------

def tilebaits(seq_array)
	# Process FASTA file
	baitsout = []
	outfilter = []
	paramline = ["Chromosome:Coordinates\tBaitLength\t%GC\tTm\tMeanQuality\tMinQuality\tKept\n"]
	coordline = []
	filtercoordline = []
	if seq_array.is_a?(String) #Read FASTA file if external file
		print "** Reading FASTA/FASTQ **\n"
		seq_array = read_fasta(seq_array)
	end
	seq_array.size.times do
		baitsout.push([])
		outfilter.push([])
		paramline.push([])
		coordline.push([])
		filtercoordline.push([])
	end
	print "** Generating and filtering baits **\n"
	threads = [] # Array to hold threading
	$options.threads.times do |i|
		threads[i] = Thread.new {
			for Thread.current[:j] in 0...seq_array.size
				if Thread.current[:j] % $options.threads == i
					Thread.current[:qual] = [$options.fasta_score]
					Thread.current[:seqst] = 1 # Beginning of bait coordinate
					while Thread.current[:seqst] < seq_array[Thread.current[:j]].seq.length+1 # Stop the loop once end of sequence reached
						Thread.current[:seqend] = Thread.current[:seqst]+$options.baitlength-1 #End of bait coordinate
						if Thread.current[:seqend] > seq_array[Thread.current[:j]].seq.length and !seq_array[Thread.current[:j]].circular #correct for running off end of linear sequence
							Thread.current[:seqend] = seq_array[Thread.current[:j]].seq.length
							Thread.current[:prb] = seq_array[Thread.current[:j]].seq[Thread.current[:seqst]-1..Thread.current[:seqend]-1] #Correct for 0-based counting
							Thread.current[:qual] = seq_array[Thread.current[:j]].numeric_quality[Thread.current[:seqst]-1..Thread.current[:seqend]-1] unless seq_array[Thread.current[:j]].fasta
						elsif Thread.current[:seqend] > seq_array[Thread.current[:j]].seq.length and seq_array[Thread.current[:j]].circular #add beginning of sequence to end
							Thread.current[:seqend] -= seq_array[Thread.current[:j]].seq.length
							Thread.current[:prb] = seq_array[Thread.current[:j]].seq[Thread.current[:seqst]-1..seq_array[Thread.current[:j]].seq.length-1] + seq_array[Thread.current[:j]].seq[0..Thread.current[:seqend]-1] #Correct for 0-based counting
							Thread.current[:qual] = seq_array[Thread.current[:j]].numeric_quality[Thread.current[:seqst]-1..seq.seq.length-1] + seq_array[Thread.current[:j]].numeric_quality[0..Thread.current[:seqend]-1] unless seq_array[Thread.current[:j]].fasta
						else
							Thread.current[:prb] = seq_array[Thread.current[:j]].seq[Thread.current[:seqst]-1..Thread.current[:seqend]-1] #Correct for 0-based counting
							Thread.current[:qual] = seq_array[Thread.current[:j]].numeric_quality[Thread.current[:seqst]-1..Thread.current[:seqend]-1] unless seq_array[Thread.current[:j]].fasta
						end
						Thread.current[:prb].gsub!("T","U") if $options.rna # RNA output handling
						Thread.current[:bait] = ">" + seq_array[Thread.current[:j]].header + "_" + Thread.current[:seqst].to_s + "-" + Thread.current[:seqend].to_s + "\n" + Thread.current[:prb] + "\n"
						Thread.current[:coord] = seq_array[Thread.current[:j]].header + "\t" + Thread.current[:seqst-1].to_s + "\t" + Thread.current[:seqend].to_s + "\n"
						baitsout[Thread.current[:j]].push(Thread.current[:bait])
						coordline[Thread.current[:j]].push(Thread.current[:coord])
						if $options.filter
							Thread.current[:flt] = filter_baits(Thread.current[:prb], Thread.current[:qual]) # U should not affection filtration
							if Thread.current[:flt][0]
								outfilter[Thread.current[:j]].push(Thread.current[:bait])
								filtercoordline[Thread.current[:j]].push(Thread.current[:coord])
							end
							if $options.params
								Thread.current[:param] = seq_array[Thread.current[:j]].header + ":" + Thread.current[:seqst].to_s + "-" + Thread.current[:seqend].to_s + "\t" + Thread.current[:flt][1]
								paramline[Thread.current[:j]+1].push(Thread.current[:param])
							end
						end
						Thread.current[:seqst] += $options.tileoffset
					end
				end
			end
		}
	end
	threads.each { |thr| thr.join }
	write_baits(baitsout, outfilter, paramline, coordline, filtercoordline)
end
