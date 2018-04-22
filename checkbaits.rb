#!/usr/bin/env ruby
#-----------------------------------------------------------------------------------------------
# checkbaits
CHECKBAITSVER = "1.1.0"
# Michael G. Campana, 2017-2018
# Smithsonian Conservation Biology Institute
#-----------------------------------------------------------------------------------------------

def checkbaits
	# Process FASTA file
	print "** Reading FASTA/FASTQ **\n"
	if $options.params
		paramline = "Bait\tBaitLength\tGC%\tTm\tMasked%\tMaxHomopolymer\tSeqComplexity\tMeanQuality\tMinQuality\tNs\tGaps\tKept"
		write_file("-filtered-params.txt", paramline)
	end
	seq_array = read_fasta($options.infile)
	print "** Filtering baits **\n"
	@splits = setup_temp(seq_array.size)
	threads = []
	@kept = []
	$options.threads.times do |i|
		@kept.push(0)
		threads[i] = Thread.new {
			for Thread.current[:j] in @splits[i] ... @splits[i+1]
				if seq_array[Thread.current[:j]].fasta
					Thread.current[:flt] = filter_baits(seq_array[Thread.current[:j]].seq, [$options.fasta_score])
				else
					Thread.current[:flt] = filter_baits(seq_array[Thread.current[:j]].seq, seq.numeric_quality)
				end
				if Thread.current[:flt][0]
					seq_array[Thread.current[:j]].seq = reversecomp(seq_array[Thread.current[:j]].seq) if $options.rc # Output reverse complemented baits if requested
					seq_array[Thread.current[:j]].seq.gsub!("T","U") if $options.rna # RNA output handling
					seq_array[Thread.current[:j]].seq.gsub!("t","u") if $options.rna # RNA output handling
					Thread.current[:bait] = ">" + seq_array[Thread.current[:j]].header + "\n" + seq_array[Thread.current[:j]].seq
					write_file("-filtered-baits.fa", Thread.current[:bait], true, i)
					@kept[i] += 1
				end
				if $options.params
					Thread.current[:param] = seq_array[Thread.current[:j]].header + "\t" + Thread.current[:flt][1]
					write_file("-filtered-params.txt", Thread.current[:param], true, i)
				end
			end
		}
	end
	threads.each { |thr| thr.join }
	cat_files
	if $options.log
		kept = 0
		for good in @kept
			kept += good
		end
		write_file(".log.txt", "TotalBaits\tFilteredBaits\tExcludedBaits\n" + seq_array.size.to_s + "\t" + kept.to_s + "\t" + (seq_array.size-kept).to_s + "\n")
	end
end
