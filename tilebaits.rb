#!/usr/bin/env ruby
#-----------------------------------------------------------------------------------------------
# tilebaits
TILEBAITSVER = "1.0.4"
# Michael G. Campana, 2017
# Smithsonian Conservation Biology Institute
#-----------------------------------------------------------------------------------------------

def tilebaits(seq_array)
	# Process FASTA file
	baitsout = []
	outfilter = []
	paramline = ["Chromosome:Coordinates\tBaitLength\tGC%\tTm\tMasked%\tMaxHomopolymer\tSeqComplexity\tMeanQuality\tMinQuality\tNs\tGaps\tKept\n"]
	coordline = []
	filtercoordline = []
	rbedline = [] # Relative bed coordinates
	rfilterline = [] # Filtered relative bed coordinates
	if seq_array.is_a?(String) #Read FASTA file if external file
		print "** Reading FASTA/FASTQ **\n"
		seq_array = read_fasta(seq_array)
	end
	if $options.log
		logs = []
		$options.logtext += "BaitCoverage\nSequence\tLength\tNumberBaits\tRetainedBaits\tExcludedBaits\tTotalBaitCoverage(×)\tFilteredBaitCoverage(×)\n"
	end
	seq_array.size.times do
		baitsout.push([])
		outfilter.push([])
		paramline.push([])
		coordline.push([])
		filtercoordline.push([])
		rbedline.push([])
		rfilterline.push([])
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
							Thread.current[:seqst] = Thread.current[:seqend] - $options.baitlength + 1 if $options.shuffle
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
						Thread.current[:prb] = reversecomp(Thread.current[:prb]) if $options.rc # Output reverse complemented baits if requested
						Thread.current[:prb].gsub!("T","U") if $options.rna # RNA output handling
						Thread.current[:prb].gsub!("t","u") if $options.rna # RNA output handling
						Thread.current[:prb] = extend_baits(Thread.current[:prb],  seq_array[Thread.current[:j]].seq, Thread.current[:seqst]-1, Thread.current[:seqend]-1) if $options.gaps == "extend"
						Thread.current[:bait] = ">" + seq_array[Thread.current[:j]].header + "_" + Thread.current[:seqst].to_s + "-" + Thread.current[:seqend].to_s + "\n" + Thread.current[:prb] + "\n"
						if $options.algorithm == "bed2baits" or $options.algorithm == "annot2baits" or $options.algorithm = "tilebaits"
							bedstart = (Thread.current[:seqst]-1 + seq_array[Thread.current[:j]].bedstart).to_s
							bedend = (Thread.current[:seqend] + seq_array[Thread.current[:j]].bedstart).to_s
							Thread.current[:coord] = seq_array[Thread.current[:j]].header + "\t" + bedstart + "\t" + bedend + "\n"
							Thread.current[:rbed] = seq_array[Thread.current[:j]].header + "\t" + (Thread.current[:seqst]-1).to_s + "\t" + Thread.current[:seqend].to_s + "\n"
						else
							Thread.current[:coord] = seq_array[Thread.current[:j]].header + "\t" + (Thread.current[:seqst]-1).to_s + "\t" + Thread.current[:seqend].to_s + "\n"
						end
						baitsout[Thread.current[:j]].push(Thread.current[:bait])
						coordline[Thread.current[:j]].push(Thread.current[:coord])
						rbedline[Thread.current[:j]].push(Thread.current[:rbed])
						if $options.filter
							Thread.current[:flt] = filter_baits(Thread.current[:prb], Thread.current[:qual]) # U should not affection filtration
							if Thread.current[:flt][0]
								outfilter[Thread.current[:j]].push(Thread.current[:bait])
								filtercoordline[Thread.current[:j]].push(Thread.current[:coord])
								rfilterline[Thread.current[:j]].push(Thread.current[:rbed])
							end
							if $options.params
								Thread.current[:param] = seq_array[Thread.current[:j]].header + ":" + Thread.current[:seqst].to_s + "-" + Thread.current[:seqend].to_s + "\t" + Thread.current[:flt][1]
								paramline[Thread.current[:j]+1].push(Thread.current[:param])
							end
						end
						Thread.current[:seqst] += $options.tileoffset
					end
					if $options.log
						logs[Thread.current[:j]] = [seq_array[Thread.current[:j]].header, seq_array[Thread.current[:j]].seq.length, baitsout[Thread.current[:j]].size]
						if $options.filter
							logs[Thread.current[:j]].push(outfilter[Thread.current[:j]].size, baitsout[Thread.current[:j]].size - outfilter[Thread.current[:j]].size, (baitsout[Thread.current[:j]].size * $options.baitlength).to_f/seq_array[Thread.current[:j]].seq.length.to_f, (outfilter[Thread.current[:j]].size * $options.baitlength).to_f/seq_array[Thread.current[:j]].seq.length.to_f)
						else
							logs[Thread.current[:j]].push(baitsout[Thread.current[:j]].size, "NA", (baitsout[Thread.current[:j]].size * $options.baitlength).to_f/seq_array[Thread.current[:j]].seq.length.to_f, "NA")
						end
					end
				end
			end
		}
	end
	threads.each { |thr| thr.join }
	if $options.log
		for i in 0 ... logs.size
			$options.logtext += logs[i].join("\t") + "\n"
		end
	end
	write_baits(baitsout, outfilter, paramline, coordline, filtercoordline, $options.outdir + "/" + $options.outprefix, rbedline, rfilterline)
end
