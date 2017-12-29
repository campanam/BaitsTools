#!/usr/bin/ruby
#-----------------------------------------------------------------------------------------------
# selectsnps 0.6
# Michael G. Campana, 2016
# Smithsonian Conservation Biology Institute
#-----------------------------------------------------------------------------------------------

# To be implemented
# Output non-reference alleles
# Include/Exclude multi-allelic sites

#-----------------------------------------------------------------------------------------------
def selectsnps
	# Read VCF file
	@snps = {}
	chromosome = ""
	temp_snps = []
	File.open($options.infile, 'r') do |snpreg|
		while line = snpreg.gets
			if $options.scale and line[0..12] == "##contig=<ID="
				reg = ""
				len = ""
				count = 13
				brk = 0
				while brk == 0 
					if line[count].chr == ","
						count += 8
						brk += 1
					else
						reg += line[count].chr
						count += 1
					end
				end
				while brk == 1
					if line[count].chr == ">"
						brk += 1
					else
						len += line[count].chr
						count += 1
					end
				end
				scaled = (len.to_i/$options.distance).floor
				$options.scalehash[reg] = scaled
			elsif line[0].chr != "#"
				reg = ""
				snp = ""
				tab = 0
				count = 0
				while tab < 1
					if line[count].chr == "\t"
						@snps[chromosome]=temp_snps
						temp_snps = [] if chromosome != reg
						chromosome = reg
						tab += 1
					else
						reg += line[count].chr
					end
					count += 1
				end
				while tab < 2
					if line[count].chr == "\t"
						tab += 1
					else
						snp += line[count].chr
					end
					count += 1
				end
				temp_snps.push(snp.to_i)		
			end
		end
	end
	@snps[chromosome]=temp_snps 		#To add last state from last line since otherwise will not loop
	@snps.delete_if {|key, value| key == ""} # Delete dummy original value
	# Select SNPs
	temp_snps = @snps.dup #Avoid messing with the original hash just in case. Re-use the temporary variable.
	@selectsnps = {}
	if !$options.every
		for i in 1..$options.totalsnps
			break if temp_snps.size == 0 # Stop selecting SNPs when all contigs deleted from consideration
			selected_contig = temp_snps.keys[rand(temp_snps.size)] # Get name of contig
			selected_snp = temp_snps[selected_contig][rand(temp_snps[selected_contig].size)]
			# Delete SNPs that are too close
			tmp = temp_snps[selected_contig].dup # Duplicate this subsection so that during deletion there are not errors
			for i in 0 .. tmp.size - 1
				if tmp[i] < selected_snp && selected_snp - tmp[i] < $options.distance
					temp_snps[selected_contig].delete(tmp[i])
				elsif tmp[i] > selected_snp && tmp[i] - selected_snp < $options.distance
					temp_snps[selected_contig].delete(tmp[i])
				end
			end
			# Add SNP to selected pool and delete contigs if maximum number of SNPs reached or no remaining SNPs on contig
			@selectsnps[selected_contig] = [] if @selectsnps[selected_contig].nil?
			@selectsnps[selected_contig].push(selected_snp)
			temp_snps[selected_contig].delete(selected_snp) # So it cannot be reselected
			if $options.scale
				maxsize = $options.scalehash[selected_contig]
			else
				maxsize = $options.maxsnps
			end
			if @selectsnps[selected_contig].size == maxsize or temp_snps[selected_contig].size == 0
				temp_snps.delete_if {|key, value | key == selected_contig}
			end
		end
	else
		@selectsnps = @snps.dup
	end
	# Write VCF & baits
	search_keys = [] # Array to hold contig/SNP unique search values
	for i in 0..@selectsnps.size - 1
		for snp in @selectsnps[@selectsnps.keys[i]]
			search = @selectsnps.keys[i] + "\t" + snp.to_s + "\t"
			search_keys.push(search)
		end
	end
	filter_search_keys = search_keys.dup # Array for filtering
	vcfout = ""
	if !$options.every
		File.open($options.infile, 'r') do |snpsel|
			while line = snpsel.gets
				if line[0].chr == "#"
					vcfout += line
				else
					search = ""
					tab = 0
					count = 0
					while tab < 2
						tab += 1 if line[count].chr == "\t"
						search += line[count].chr
						count += 1
					end
					if search_keys.include?(search)
						vcfout += line
						search_keys.delete(search)
					end
				end
			end
		end
		File.open($options.infile + "-selected.vcf", 'w') do |write|
			write.puts vcfout
		end
	end
	# Output bait sequences if this option was selected
	if $options.baits
		baitsout = ""
		outfilter = ""
		paramline = "Chromosome:Coordinates\tSNP\tBaitLength\t%GC\tTm\tKept\n"
		coordline = ""
		filtercoordline = ""
		refseq = read_fasta($options.refseq)
		for rseq in refseq
			if @selectsnps.keys.include?(rseq.header)
				for snp in @selectsnps[rseq.header]
					if $options.tiling
						tiling = $options.tiledepth
					else
						tiling = 1
						$options.tileoffset = 1
					end
					for tile in 0...tiling
						be4 = snp - $options.lenbef + tile * $options.tileoffset
						after = snp + $options.lenaft + tile * $options.tileoffset
						be4 = 1 if (be4 < 1 and !rseq.circular)
						if after > rseq.seq.length and !rseq.circular
							after = rseq.seq.length
							if be4 < 1
								prb = rseq.seq[be4-1..-1] + rseq.seq[0..after-1]  #Correct for 0-based counting
							else
								prb = rseq.seq[be4-1..after-1]  #Correct for 0-based counting
							end 
						elsif after > rseq.seq.length and rseq.circular
							after -= rseq.seq.length
							if be4 < 1
								prb = rseq.seq[be4-1..-1] + rseq.seq[be4-1..rseq.seq.length-1]+rseq.seq[0..after-1]
							else
								prb = rseq.seq[be4-1..rseq.seq.length-1]+rseq.seq[0..after-1]  #Correct for 0-based counting and end of sequence
							end
						else
							if be4 < 1
								prb = rseq.seq[be4-1..-1] + rseq.seq[0..after-1]
							else
								prb = rseq.seq[be4-1..after-1]  #Correct for 0-based counting
							end
						end
						seq = ">" + rseq.header + "\t" + snp.to_s + "\n" + prb + "\n"
						baitsout += seq
						be4 = rseq.seq.length + be4 if be4 < 1
						coord = rseq.header + ":" + be4.to_s + "-" + after.to_s + "\n"
						coordline += coord
						if $options.filter
							parameters = filter_baits(prb) #Correct for 0-based counting and end of sequence
							if parameters[0]
								outfilter += seq
								filtercoordline += coord 
							else
								search = rseq.header + "\t" + snp.to_s + "\t"
								filter_search_keys.delete(search)
							end
							if $options.params
								paramline += rseq.header + ":" + be4.to_s + "-" + after.to_s + "\t" + snp.to_s + "\t" + parameters[1]
							end
						end
					end
				end
			end
		end
		write_baits(baitsout, outfilter, paramline, coordline, filtercoordline)
		if $options.filter
			vcffilt = ""
			File.open($options.infile, 'r') do |snpsel|
				while line = snpsel.gets
					if line[0].chr == "#"
						vcffilt += line
					else
						search = ""
						tab = 0
						count = 0
						while tab < 2
							tab += 1 if line[count].chr == "\t"
							search += line[count].chr
							count += 1
						end
						if filter_search_keys.include?(search)
							vcffilt += line
							filter_search_keys.delete(search)
						end
					end
				end
			end
			File.open($options.infile + "-filtered.vcf", 'w') do |write|
				write.puts vcffilt
			end
		end
	end
end
