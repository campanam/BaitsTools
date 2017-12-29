#!/usr/bin/env ruby
#-----------------------------------------------------------------------------------------------
# vcf2baits
VCF2BAITSVER = "0.10"
# Michael G. Campana, 2017
# Smithsonian Conservation Biology Institute
#-----------------------------------------------------------------------------------------------

def vcf2baits
	# Read VCF file
	@snps = {}
	chromosome = ""
	temp_snps = []
	vcfout = ""
	columns = ""
	print "** Reading VCF **\n"
	File.open($options.infile, 'r') do |snpreg|
		while line = snpreg.gets
			vcfout += line if line[0..1] == "##" # Add info lines to header
			columns = line if line[0..1] == "#C" # Get column header separately so I can insert new information
			if $options.scale and line[0..12] == "##contig=<ID="
				line_arr = line.split(",")
				reg = line_arr[0][13..-1] # Get region ID, but remove tag
				len = line_arr[1].split("=")[1].to_i # Remove length= tag, noting that spacing can be variable
				scaled = (len/$options.distance).floor
				$options.scalehash[reg] = scaled
			elsif line[0].chr != "#"
				line_arr = line.split("\t")
				reg = line_arr[0]
				snp = line_arr[1].to_i
				ref = line_arr[3]
				alt = line_arr[4].split(",")
				qual = line_arr[5].to_i
				@snps[chromosome]=temp_snps
				temp_snps = [] if chromosome != reg # Empty set if no longer on same chromosome
				chromosome = reg # Reset chromosome
				if !$options.varqual_filter or qual > $options.varqual
					temp_snps.push(Chromo_SNP.new(reg, snp, [], ref, alt, qual, line)) #Now using Chromo_SNPs to be consistent
				end
			end
		end
	end
	@snps[chromosome]=temp_snps 		#To add last state from last line since otherwise will not loop
	@snps.delete_if {|key, value| key == ""} # Delete dummy original value
	@snps.delete_if {|key, value| value == []} # Delete empty contigs due to QUAL filter
	print "** Selecting variants **\n" unless $options.every
	$options.logtext += "Variants\n" if $options.log
	@selectsnps = selectsnps(@snps) 	# Select SNPs
	filteredsnps = @selectsnps.dup # Hash for filtering
	# Write VCF & baits
	vcfout += "##baitstools_vcf2baitsVersion=" + VCF2BAITSVER + "+BaitsTools-" + BAITSTOOLSVER + "\n##baitstools_vcf2baitsCommand="
	cmdline = get_command_line
	vcffilt = vcfout + cmdline[0] + cmdline[1] + "\n" + columns
	vcfout += cmdline[0] + "\n" + columns
	if !$options.every
		for i in 0...@selectsnps.size
			for snp in @selectsnps[@selectsnps.keys[i]]
				vcfout += snp.line
			end
		end
		File.open($options.outdir + "/" + $options.outprefix + "-selected.vcf", 'w') do |write|
			write.puts vcfout
		end
	end
	# Output bait sequences unless -p
	if !$options.no_baits
		print "** Reading reference sequence **\n"
		refseq = read_fasta($options.refseq)
		print "** Generating and filtering baits **\n"
		$options.logtext += "VariantBaits\n" if $options.log
		baits = snp_to_baits(@selectsnps, refseq)
		write_baits(baits[0], baits[1], baits[2], baits[3], baits[4])
		if $options.filter
			for i in 0...baits[5].size # Filtered snps
				for snp in baits[5][baits[5].keys[i]]
					vcffilt += snp.line
				end
			end
			File.open($options.outdir + "/" + $options.outprefix + "-filtered.vcf", 'w') do |write|
				write.puts vcffilt
			end
		end
	end
end
