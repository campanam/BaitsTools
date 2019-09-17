#!/usr/bin/env ruby
#-----------------------------------------------------------------------------------------------
# aln2baits
ALN2BAITSVER = "1.6.0"
# Michael G. Campana, 2017-2019
# Smithsonian Conservation Biology Institute
#-----------------------------------------------------------------------------------------------

class Hap_Window # Object defining a haplotype window
	attr_accessor :header, :seqstart, :seqend, :haplotypes, :bedstarts, :locus
	def initialize(header, seqstart, seqend, haplotypes = [], bedstarts = [], locus = "")
		@header = header
		@seqstart = seqstart
		@seqend = seqend
		@haplotypes = haplotypes
		@bedstarts = bedstarts
		@locus = locus
	end
	def var_permutations(aln) # Get possible variant permutations
		variants = []
		for i in 0...@haplotypes[0].length
			variants.push([])
		end
		threads = [] # Array to hold threads
		$options.threads.times do |j|
			threads[j] = Thread.new {
				for Thread.current[:i] in 0...@haplotypes[0].length
					if Thread.current[:i] % $options.threads == j
						Thread.current[:vars] = []
						for Thread.current[:seq] in @haplotypes
							Thread.current[:vars].push(get_ambiguity(Thread.current[:seq][Thread.current[:i]]))
						end
						Thread.current[:vars] = Thread.current[:vars].flatten!.uniq # Remove duplicate variants
						variants[Thread.current[:i]] = Thread.current[:vars] # Minimize lock time
					end
				end
			}
		end
		threads.each { |thr| thr.join }
		varindex = 1 # Index for sequence
		for var in variants
			varindex *= var.size # Must be complete down here or interferes with multithreading
		end
		revised_haplos = []
		bedstarts = @bedstarts[0] # Reset bedstart array and assume coordinates of first array member
		@bedstarts = []
		for i in 1..varindex
			revised_haplos.push("")
			@bedstarts.push(bedstarts)
		end
		threads = []
		$options.threads.times do |j|
			threads[j] = Thread.new {
				for Thread.current[:k] in 0...varindex
					if Thread.current[:k] % $options.threads == j
						for Thread.current[:i] in 0...@haplotypes[0].length
							Thread.current[:var] = Thread.current[:k] % variants[Thread.current[:i]].size
							revised_haplos[Thread.current[:k]] << variants[Thread.current[:i]][Thread.current[:var]] # Minimize lock time
						end
						if $options.gaps == "extend"
							Thread.current[:refhap] = aln[rand(aln.size)] # Choose a random reference sequence
							revised_haplos[Thread.current[:k]] = extend_baits(revised_haplos[Thread.current[:k]], $options.fasta_score, Thread.current[:refhap], @seqstart, @seqend)[0]
						end
						revised_haplos[Thread.current[:k]] = fill_in_baits(revised_haplos[Thread.current[:k]], $options.fasta_score, true)[0] if $options.fillin_switch
					end
				end
			}
		end
		threads.each { |thr| thr.join }
		@haplotypes = revised_haplos.uniq # Remove new redundant haplotypes
	end
end
#-----------------------------------------------------------------------------------------------
def aln2baits(aln)
	if aln.is_a?(String) #Read FASTA file if external file
		print "** Reading alignment **\n"
		aln = read_fasta(aln)
	end
	print "** Sorting sequences by locus ** \n"
	aln_hash = {}
	for seq in aln
		if aln_hash[seq.locus].nil?
			aln_hash[seq.locus] = [seq]
		else
			aln_hash[seq.locus].push(seq)
		end
	end
	threads = []
	print "** Identifying haplotypes **\n"
	@windows = [] # Array to hold windows regardless of locus
	# Get haplotypes
	$options.threads.times do |i|
		threads[i] = Thread.new {
			for Thread.current[:j] in 0 ... aln_hash.keys.size
				if Thread.current[:j] % $options.threads == i
					Thread.current[:seqstart] = 0
					Thread.current[:windows] = []
					Thread.current[:aln] = aln_hash[aln_hash.keys[Thread.current[:j]]]
					while Thread.current[:seqstart] < Thread.current[:aln][0].seq.length
						Thread.current[:seqend] = Thread.current[:seqstart] + $options.baitlength - 1
						Thread.current[:seqcycles] = -1 # Number of complete cycles through short baits. First cycle through brings to 0
						if Thread.current[:seqend] > Thread.current[:aln][0].seq.length-1 and !Thread.current[:aln][0].circular
							Thread.current[:seqend] = Thread.current[:aln][0].seq.length - 1
							if $options.shuffle
								Thread.current[:seqstart] = Thread.current[:seqend] - $options.baitlength + 1
								Thread.current[:seqstart] = 0 if Thread.current[:seqstart] < 0
							end
						elsif Thread.current[:seqend] > Thread.current[:aln][0].seq.length-1 and Thread.current[:aln][0].circular
							while Thread.current[:seqend] > Thread.current[:aln][0].seq.length-1
								Thread.current[:seqcycles] += 1
								Thread.current[:seqend] -= Thread.current[:aln][0].seq.length
							end
						end
						Thread.current[:window] = Hap_Window.new([], Thread.current[:seqstart], Thread.current[:seqend])
						Thread.current[:window].locus = Thread.current[:aln][0].locus
						for Thread.current[:seq] in Thread.current[:aln]
							if Thread.current[:seqcycles] == -1
								Thread.current[:tmp_seq] = Thread.current[:seq].seq[Thread.current[:seqstart]..Thread.current[:seqend]]
							else
								Thread.current[:tmp_seq] = Thread.current[:seq].seq[Thread.current[:seqstart]..-1]
								Thread.current[:seqcycles].times { Thread.current[:tmp_seq] << Thread.current[:seq].seq }
								Thread.current[:tmp_seq] << Thread.current[:seq].seq[0..Thread.current[:seqend]]
							end
							Thread.current[:tmp_seq].upcase! unless $options.maxmask_filter # Treat all cases the same unless the masking filter is requested
							unless Thread.current[:window].haplotypes.include?(Thread.current[:tmp_seq]) # Add new sequences to haplotype list. Extra include test prevents some baits from being extended/filled-in multiple times.
								if $options.haplodef == "haplotype"
									Thread.current[:tmp_seq] = extend_baits(Thread.current[:tmp_seq], $options.fasta_score, Thread.current[:seq], Thread.current[:seqstart], Thread.current[:seqend])[0] if $options.gaps == "extend"
									Thread.current[:tmp_seq] = fill_in_baits(Thread.current[:tmp_seq], $options.fasta_score, true)[0] if $options.fillin_switch
								end
								unless Thread.current[:window].haplotypes.include?(Thread.current[:tmp_seq]) # Sequence revisions can make redundant baits
									Thread.current[:window].haplotypes.push(Thread.current[:tmp_seq])
									Thread.current[:window].bedstarts.push(Thread.current[:seq].bedstart)
									Thread.current[:window].header.push(Thread.current[:seq].header)
								end
							end
						end
						Thread.current[:window].var_permutations(aln) if $options.haplodef == "variant" # Call here for code efficiency
						Thread.current[:windows].push(Thread.current[:window])
						Thread.current[:seqstart] += $options.tileoffset # Control the window by tiling density
					end
					@windows[Thread.current[:j]] = Thread.current[:windows]
				end
			end
		}
	end
	threads.each { |thr| thr.join }
	@windows.flatten!
	# Get bait candidates
	print "** Generating and filtering baits **\n"
	if $options.params
		paramline = "Sequence:Locus:Coordinates:Haplotype\tBaitLength\tGC%\tTm\tMasked%\tMaxHomopolymer\tSeqComplexity\tMeanQuality\tMinQuality\tNs\tGaps\tKept"
		write_file("-filtered-params.txt", paramline)
	end
	if $options.log
		logs = []
		write_file(".log.txt", "Windows\nLocus\tWindowStart\tWindowEnd\tNumberHaplotypes\tRetainedBaits\tExcludedBaits") 
	end
	threads = []
	@splits = setup_temp(@windows.size)
	$options.used_threads.times do |i|
		threads[i] = Thread.new {
			for Thread.current[:j] in @splits[i] ... @splits[i+1]
				Thread.current[:filtnum] = 0
				for Thread.current[:hapno] in 1..@windows[Thread.current[:j]].haplotypes.size
					Thread.current[:prb] = @windows[Thread.current[:j]].haplotypes[Thread.current[:hapno]-1]
					Thread.current[:rng] = (@windows[Thread.current[:j]].seqstart+1).to_s+"-"+(@windows[Thread.current[:j]].seqend+1).to_s # Adjust for 1-based indexing
					Thread.current[:prb] = reversecomp(Thread.current[:prb])[0] if $options.rc  # Output reverse complemented baits if requested
					Thread.current[:prb] = make_rna(Thread.current[:prb]) if $options.rna # Correct RNA
					Thread.current[:completeprb] = $options.fiveprime + Thread.current[:prb] + $options.threeprime
					Thread.current[:prb] = Thread.current[:completeprb] unless $options.noaddenda
					$options.haplodef == "haplotype" ? Thread.current[:header] = @windows[Thread.current[:j]].header[Thread.current[:hapno]-1] : Thread.current[:header] = "Alignment"
					Thread.current[:bait] = ">" + Thread.current[:header] +"_locus" + @windows[Thread.current[:j]].locus + "_" + Thread.current[:rng] + "_haplotype" + Thread.current[:hapno].to_s + "\n" + Thread.current[:completeprb]
					Thread.current[:bedstart] = (@windows[Thread.current[:j]].seqstart + @windows[Thread.current[:j]].bedstarts[Thread.current[:hapno]-1]).to_s
					Thread.current[:bedend] =  (@windows[Thread.current[:j]].seqend+1 + @windows[Thread.current[:j]].bedstarts[Thread.current[:hapno]-1]).to_s
					Thread.current[:coord] = Thread.current[:header] + "\t" + Thread.current[:bedstart] + "\t" + Thread.current[:bedend]
					Thread.current[:rbed] = Thread.current[:header] + "\t" + @windows[Thread.current[:j]].seqstart.to_s + "\t" + (@windows[Thread.current[:j]].seqend+1).to_s
					write_file("-baits.fa", Thread.current[:bait], true, i)
					write_file("-baits.bed", Thread.current[:coord], true, i) if $options.coords
					write_file("-baits-relative.bed", Thread.current[:rbed], true, i) if $options.rbed
					if $options.filter
						Thread.current[:flt] = filter_baits(Thread.current[:prb]) # U won't affect filtration
						if Thread.current[:flt][0]
							Thread.current[:filtnum] += 1
							write_file("-filtered-baits.fa", Thread.current[:bait], true, i)
							write_file("-filtered-baits.bed", Thread.current[:coord], true, i) if $options.coords
							write_file("-filtered-baits-relative.bed", Thread.current[:rbed], true, i) if $options.rbed
						end
						if $options.params
							Thread.current[:param] = Thread.current[:header] + ":" + @windows[Thread.current[:j]].locus + Thread.current[:rng] + ":" + Thread.current[:hapno].to_s + "\t" + Thread.current[:flt][1]
							write_file("-filtered-params.txt", Thread.current[:param], true, i)
						end
					end
				end
				if $options.log
					logs[Thread.current[:j]] = [@windows[Thread.current[:j]].locus,@windows[Thread.current[:j]].seqstart+1, @windows[Thread.current[:j]].seqend+1, @windows[Thread.current[:j]].haplotypes.size]
					if $options.filter
						logs[Thread.current[:j]].push(Thread.current[:filtnum], @windows[Thread.current[:j]].haplotypes.size - Thread.current[:filtnum])
					else
						logs[Thread.current[:j]].push(@windows[Thread.current[:j]].haplotypes.size, "NA")
					end
				end
			end
		}
	end
	threads.each { |thr| thr.join }
	cat_files
	if $options.log
		vlogs = [[],[],[],[]]
		for i in 0 ... logs.size
			write_file(".log.txt", logs[i].join("\t"))
			vlogs[0].push(logs[i][1])
			vlogs[1].push(logs[i][2])
			vlogs[2].push(logs[i][3])
			vlogs[3].push(logs[i][4])
		end
		write_file(".log.txt", "\nAlignmentLength\tTotalBaitCoverage(x)\tFilteredBaitCoverage(x)")
		if $options.filter
			write_file(".log.txt", (vlogs[1].max-vlogs[0].min+1).to_s + "\t" + mean(vlogs[2]).to_s + "\t" + mean(vlogs[3]).to_s)
		else
			write_file(".log.txt", (vlogs[1].max-vlogs[0].min+1).to_s + "\t" + mean(vlogs[2]).to_s + "\tNA")
		end
	end
end
