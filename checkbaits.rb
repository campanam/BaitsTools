#!/usr/bin/env ruby
#-----------------------------------------------------------------------------------------------
# checkbaits
CHECKBAITSVER = "1.7.0"
# Michael G. Campana, 2017-2021
# Smithsonian Conservation Biology Institute
#-----------------------------------------------------------------------------------------------

def multi_checkbaits(seq_array, altbait = nil)
	logtext_infix = altbaits_infix(altbait)
	print "** Filtering " + $options.baitlength.to_s + " baits **\n"
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
				if altbait != nil # Adjust lengths for multi-length probes by cutting bases off 3' end
					Thread.current[:completeprb] = Thread.current[:completeprb][0...altbait]
					Thread.current[:prb] = Thread.current[:prb][0...altbait]
					Thread.current[:qual]= Thread.current[:qual][0...altbait]
				end
				Thread.current[:bait] = ">" + seq_array[Thread.current[:j]].header + "\n" + Thread.current[:completeprb]
				Thread.current[:flt] = filter_baits(Thread.current[:prb], Thread.current[:qual])
				if Thread.current[:flt][0]
					write_file("-filtered-baits.fa", Thread.current[:bait], true, i)
					@kept[i] += 1
					unless $options.inbed.nil?
						Thread.current[:bed] = seq_array[Thread.current[:j]].header + "\t" + $options.inbed[seq_array[Thread.current[:j]].header].join("\t")
						write_file("-filtered-baits.bed", Thread.current[:bed], true, i)
					end
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
		write_file(".log.txt", logtext_infix + "TotalBaits\tFilteredBaits\tExcludedBaits\n" + seq_array.size.to_s + "\t" + kept.to_s + "\t" + (seq_array.size-kept).to_s + "\n")
	end
end
#-----------------------------------------------------------------------------------------------
def checkbaits
	# Process FASTA file
	print "** Reading FASTA/FASTQ **\n"
	seq_array = read_fasta($options.infile)
	if $options.params
		paramline = "Bait\tBaitLength\tGC%\tTm\tMasked%\tMaxHomopolymer\tSeqComplexity\tMeanQuality\tMinQuality\tNs\tGaps\tKept"
		write_file("-filtered-params.txt", paramline)
	end
	unless $options.inbed.nil?
		bed_hash = {} # Hash of BED entries and their corresponding bait sequences
		print "** Reading BED file **\n"
		gz_file_open($options.inbed) do |inbed|
			while line = inbed.gets
				line_arr = line.split
				chromo = line_arr[0]
				seqst = line_arr[1].to_i
				seqend = line_arr[2].to_i
				bed_hash[chromo] = [seqst, seqend]
			end
		end
		$options.inbed = bed_hash # No longer need original file name, so limits having to pass variables
	end
	unless $options.altbaits.nil? # Tile baits under alternate lengths
		multi_checkbaits(seq_array, $options.baitlength)
		for altbait in $options.altbaits
			$options.baitlength = altbait
			unless $options.inbed.nil? # Update BED file for different sequence lengths
				$options.inbed = {} # Reset for later values
				bed_hash.each do | key, value |
					if value[0] + altbait > value[1]
						$options.inbed[key] = [value[0], value[1]]
					else
						$options.inbed[key] = [value[0], value[0] + altbait]
					end
				end
			end
			write_file(".log.txt", "") if $options.log # Add a linebreak between subsequent entries
			multi_checkbaits(seq_array, altbait)
		end
	else
		multi_checkbaits(seq_array) # Tile baits under original settings without infix
	end
end