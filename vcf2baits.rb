#!/usr/bin/ruby
#-----------------------------------------------------------------------------------------------
# vcf2baits 0.7
## Michael G. Campana, 2016
# Smithsonian Conservation Biology Institute
#-----------------------------------------------------------------------------------------------

# To be implemented
# Output non-reference alleles
# Include/Exclude multi-allelic sites

#-----------------------------------------------------------------------------------------------
def vcf2baits
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
	@selectsnps = selectsnps(@snps) 	# Select SNPs
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
		baits = snp_to_baits(@selectsnps, filter_search_keys)
		write_baits(baits[0], baits[1], baits[2], baits[3], baits[4])
		if $options.filter
			filter_search_keys = baits[5]
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
