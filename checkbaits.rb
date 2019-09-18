#!/usr/bin/env ruby
#-----------------------------------------------------------------------------------------------
# checkbaits
CHECKBAITSVER = "1.6.0"
# Michael G. Campana, 2017-2019
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
	$options.used_threads.times do |i|
		@kept.push(0)
		threads[i] = Thread.new {
			for Thread.current[:j] in @splits[i] ... @splits[i+1]
				Thread.current[:prb] = seq_array[Thread.current[:j]].seq
				seq_array[Thread.current[:j]].fasta ? Thread.current[:qual] = [$options.fasta_score] : Thread.current[:qual] = seq_array[Thread.current[:j]].numeric_quality
				Thread.current[:prb], Thread.current[:qual], Thread.current[:completeprb] = revise_baits(Thread.current[:prb], Thread.current[:qual], seq_array[Thread.current[:j]],  Thread.current[:seqst], Thread.current[:seqend])
				Thread.current[:bait] = ">" + seq_array[Thread.current[:j]].header + "\n" + Thread.current[:completeprb]
				Thread.current[:flt] = filter_baits(Thread.current[:prb], Thread.current[:qual])
				if Thread.current[:flt][0]
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
