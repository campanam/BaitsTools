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
	$options.threads.times do |i|
		@kept.push(0)
		threads[i] = Thread.new {
			for Thread.current[:j] in @splits[i] ... @splits[i+1]
				Thread.current[:prb] = seq_array[Thread.current[:j]].seq
				seq_array[Thread.current[:j]].fasta ? Thread.current[:qual] = [$options.fasta_score] : Thread.current[:qual] = seq_array[Thread.current[:j]].numeric_quality
				Thread.current[:prb], Thread.current[:qual] = fill_in_baits(Thread.current[:prb], Thread.current[:qual], seq_array[Thread.current[:j]].fasta) if $options.fillin_switch
				Thread.current[:prb], Thread.current[:qual] = reversecomp(Thread.current[:prb], Thread.current[:qual]) if $options.rc # Output reverse complemented baits if requested
				Thread.current[:prb] = make_rna(seq_array[Thread.current[:j]].seq) if $options.rna # RNA output handling
				Thread.current[:completeprb] = $options.fiveprime + Thread.current[:prb] + $options.threeprime
				unless $options.noaddenda
					Thread.current[:prb] = Thread.current[:completeprb]
					Thread.current[:qual] = addend_qual(Thread.current[:qual]) unless seq_array[Thread.current[:j]].fasta
				end
				Thread.current[:bait] = ">" + seq_array[Thread.current[:j]].header + "\n" + Thread.current[:completeprb]
				Thread.current[:flt] = filter_baits(seq_array[Thread.current[:prb], Thread.current[:qual])
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
