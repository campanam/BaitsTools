#!/usr/bin/env ruby
#-----------------------------------------------------------------------------------------------
# checkbaits
CHECKBAITSVER = "0.4"
# Michael G. Campana, 2017
# Smithsonian Conservation Biology Institute
#-----------------------------------------------------------------------------------------------

def checkbaits
	# Process FASTA file
	print "** Reading FASTA/FASTQ **\n"
	outfilter = []
	paramline = ["Bait\tBaitLength\t%GC\tTm\tMeanQuality\tMinQuality\tKept\n"]
	seq_array = read_fasta($options.infile)
	seq_array.size.times do
		outfilter.push([])
		paramline.push([])
	end
	print "** Filtering baits **\n"
	threads = []
	$options.threads.times do |i|
		threads[i] = Thread.new {
			for Thread.current[:j] in 0...seq_array.size
				if Thread.current[:j] % $options.threads == i
					if seq_array[Thread.current[:j]].fasta
						Thread.current[:flt] = filter_baits(seq_array[Thread.current[:j]].seq, [$options.fasta_score])
					else
						Thread.current[:flt] = filter_baits(seq_array[Thread.current[:j]].seq, seq.numeric_quality)
					end
					if Thread.current[:flt][0]
						seq_array[Thread.current[:j]].seq.gsub!("T","U") if $options.rna # RNA output handling
						Thread.current[:bait] = ">" + seq_array[Thread.current[:j]].header + "\n" + seq_array[Thread.current[:j]].seq + "\n"
						outfilter[Thread.current[:j]].push(Thread.current[:bait])
					end
					if $options.params
						Thread.current[:param] = seq_array[Thread.current[:j]].header + "\t" + Thread.current[:flt][1]
						paramline[Thread.current[:j]+1].push(Thread.current[:param])
					end
				end
			end
		}
	end
	threads.each { |thr| thr.join }
	write_baits([""], outfilter, paramline)
end
