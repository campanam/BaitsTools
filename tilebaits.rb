#!/usr/bin/env ruby
#-----------------------------------------------------------------------------------------------
# tilebaits
TILEBAITSVER = "1.4.0"
# Michael G. Campana, 2017-2018
# Smithsonian Conservation Biology Institute
#-----------------------------------------------------------------------------------------------

def tilebaits(seq_array)
	# Process FASTA file
	if $options.params
		paramline = "Chromosome:Coordinates\tBaitLength\tGC%\tTm\tMasked%\tMaxHomopolymer\tSeqComplexity\tMeanQuality\tMinQuality\tNs\tGaps\tKept"
		write_file("-filtered-params.txt", paramline)
	end
	if seq_array.is_a?(String) #Read FASTA file if external file
		print "** Reading FASTA/FASTQ **\n"
		seq_array = read_fasta(seq_array)
	end
	if $options.log
		logtext = "BaitCoverage\nSequence\tLength\tNumberBaits\tRetainedBaits\tExcludedBaits\tTotalBaitCoverage(x)\tFilteredBaitCoverage(x)"
		write_file(".log.txt", logtext)
	end
	@splits = setup_temp(seq_array.size)
	print "** Generating and filtering baits **\n"
	threads = [] # Array to hold threading
	$options.used_threads.times do |i|
		threads[i] = Thread.new {
			for Thread.current[:j] in @splits[i] ... @splits[i+1]
				Thread.current[:qual] = [$options.fasta_score]
				Thread.current[:seqst] = 1 # Beginning of bait coordinate
				Thread.current[:baitnum] = 0
				Thread.current[:filtnum] = 0
				while Thread.current[:seqst] < seq_array[Thread.current[:j]].seq.length+1 # Stop the loop once end of sequence reached
					Thread.current[:seqend] = Thread.current[:seqst]+$options.baitlength-1 #End of bait coordinate
					if Thread.current[:seqend] > seq_array[Thread.current[:j]].seq.length and !seq_array[Thread.current[:j]].circular #correct for running off end of linear sequence
						Thread.current[:seqend] = seq_array[Thread.current[:j]].seq.length
						Thread.current[:seqst] = Thread.current[:seqend] - $options.baitlength + 1 if $options.shuffle
						Thread.current[:prb] = seq_array[Thread.current[:j]].seq[Thread.current[:seqst]-1..Thread.current[:seqend]-1] #Correct for 0-based counting
						Thread.current[:qual] = seq_array[Thread.current[:j]].numeric_quality[Thread.current[:seqst]-1..Thread.current[:seqend]-1] unless seq_array[Thread.current[:j]].fasta
					elsif Thread.current[:seqend] > seq_array[Thread.current[:j]].seq.length and seq_array[Thread.current[:j]].circular #add beginning of sequence to end
						Thread.current[:prb] = ""
						while Thread.current[:seqend] > seq_array[Thread.current[:j]].seq.length
							Thread.current[:prb] << seq_array[Thread.current[:j]].seq[Thread.current[:seqst]-1..seq_array[Thread.current[:j]].seq.length-1] #
							Thread.current[:prb] << seq_array[Thread.current[:j]].seq[0..Thread.current[:seqst]-2] if Thread.current[:seqst] > 1 # Re-add beginning section but don't go backwards
							Thread.current[:qual] = seq_array[Thread.current[:j]].numeric_quality[Thread.current[:seqst]-1..seq_array[Thread.current[:j]].seq.length-1] + seq_array[Thread.current[:j]].numeric_quality[0..Thread.current[:seqst]-2] unless seq_array[Thread.current[:j]].fasta
							Thread.current[:seqend] -= seq_array[Thread.current[:j]].seq.length
						end
						Thread.current[:prb] << seq_array[Thread.current[:j]].seq[Thread.current[:seqst]-1..Thread.current[:seqend]-1] #Correct for 0-based counting
						Thread.current[:qual] << seq_array[Thread.current[:j]].numeric_quality[Thread.current[:seqst]-1..Thread.current[:seqend]-1] unless seq_array[Thread.current[:j]].fasta
					else
						Thread.current[:prb] = seq_array[Thread.current[:j]].seq[Thread.current[:seqst]-1..Thread.current[:seqend]-1] #Correct for 0-based counting
						Thread.current[:qual] = seq_array[Thread.current[:j]].numeric_quality[Thread.current[:seqst]-1..Thread.current[:seqend]-1] unless seq_array[Thread.current[:j]].fasta
					end
					Thread.current[:prb] = reversecomp(Thread.current[:prb]) if $options.rc # Output reverse complemented baits if requested
					Thread.current[:prb].gsub!("T","U") if $options.rna # RNA output handling
					Thread.current[:prb].gsub!("t","u") if $options.rna # RNA output handling
					Thread.current[:prb] = extend_baits(Thread.current[:prb],  seq_array[Thread.current[:j]].seq, Thread.current[:seqst]-1, Thread.current[:seqend]-1, seq_array[Thread.current[:j]].circular) if $options.gaps == "extend"
					Thread.current[:bait] = ">" + seq_array[Thread.current[:j]].header + "_" + Thread.current[:seqst].to_s + "-" + Thread.current[:seqend].to_s + "\n" + Thread.current[:prb]
					if $options.algorithm == "bed2baits" or $options.algorithm == "annot2baits" or $options.algorithm = "tilebaits"
						bedstart = (Thread.current[:seqst]-1 + seq_array[Thread.current[:j]].bedstart).to_s
						bedend = (Thread.current[:seqend] + seq_array[Thread.current[:j]].bedstart).to_s
						Thread.current[:coord] = seq_array[Thread.current[:j]].header + "\t" + bedstart + "\t" + bedend
						Thread.current[:rbed] = seq_array[Thread.current[:j]].header + "\t" + (Thread.current[:seqst]-1).to_s + "\t" + Thread.current[:seqend].to_s
					else
						Thread.current[:coord] = seq_array[Thread.current[:j]].header + "\t" + (Thread.current[:seqst]-1).to_s + "\t" + Thread.current[:seqend].to_s
					end
					Thread.current[:baitnum] += 1
					write_file("-baits.fa", Thread.current[:bait], true, i)
					write_file("-baits.bed", Thread.current[:coord], true, i) if $options.coords
					write_file("-baits-relative.bed", Thread.current[:rbed], true, i) if $options.rbed
					if $options.filter
						Thread.current[:flt] = filter_baits(Thread.current[:prb], Thread.current[:qual]) # U should not affection filtration
						if Thread.current[:flt][0]
							Thread.current[:filtnum] += 1
							write_file("-filtered-baits.fa", Thread.current[:bait], true, i)
							write_file("-filtered-baits.bed", Thread.current[:coord], true, i) if $options.coords
							write_file("-filtered-baits-relative.bed", Thread.current[:rbed], true, i) if $options.rbed
						end
						if $options.params
							Thread.current[:param] = seq_array[Thread.current[:j]].header + ":" + Thread.current[:seqst].to_s + "-" + Thread.current[:seqend].to_s + "\t" + Thread.current[:flt][1]
							write_file("-filtered-params.txt", Thread.current[:param], true, i)
						end
					end
					Thread.current[:seqst] += $options.tileoffset
				end
				if $options.log
					Thread.current[:log] = [seq_array[Thread.current[:j]].header, seq_array[Thread.current[:j]].seq.length, Thread.current[:baitnum]]
					if $options.filter
						Thread.current[:log].push(Thread.current[:filtnum], Thread.current[:baitnum] - Thread.current[:filtnum], (Thread.current[:baitnum] * $options.baitlength).to_f/seq_array[Thread.current[:j]].seq.length.to_f, (Thread.current[:filtnum] * $options.baitlength).to_f/seq_array[Thread.current[:j]].seq.length.to_f)
					else
						Thread.current[:log].push(Thread.current[:baitnum], "NA", (Thread.current[:baitnum] * $options.baitlength).to_f/seq_array[Thread.current[:j]].seq.length.to_f, "NA")
					end
					write_file(".log.txt", Thread.current[:log].join("\t"), true, i)
				end
			end
		}
	end
	threads.each { |thr| thr.join }
	cat_files
end
