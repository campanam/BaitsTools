#!/usr/bin/env ruby
#-----------------------------------------------------------------------------------------------
# vcf2baits
VCF2BAITSVER = "1.1.0"
# Michael G. Campana, 2017-2018
# Smithsonian Conservation Biology Institute
#-----------------------------------------------------------------------------------------------

def vcf2baits
	# Read VCF file
	@snps = {}
	chromosome = ""
	vcfout = ""
	temp_snps = []
	columns = ""
	print "** Reading VCF **\n"
	File.open($options.infile, 'r') do |snpreg|
		while line = snpreg.gets
			vcfout << line if line[0..1] == "##" # Add info lines to header
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
	write_file(".log.txt", "Variants") if $options.log
	@selectsnps = selectsnps(@snps) 	# Select SNPs
	filteredsnps = @selectsnps.dup # Hash for filtering
	# Write VCF & baits
	vcfout = "##baitstools_vcf2baitsVersion=" + VCF2BAITSVER + "+BaitsTools-" + BAITSTOOLSVER + "\n##baitstools_vcf2baitsCommand="
	cmdline = get_command_line
	vcffilt = vcfout + cmdline[0] + cmdline[1] + "\n" + columns
	if !$options.every
		write_file("-selected.vcf", vcfout + cmdline[0] + "\n" + columns )
		for i in 0 ... @selectsnps.size
			for snp in @selectsnps[@selectsnps.keys[i]]
				write_file("-selected.vcf", snp.line)
			end
		end
	end
	# Output bait sequences unless -p
	if !$options.no_baits
		print "** Reading reference sequence **\n"
		refseq = read_fasta($options.refseq)
		print "** Generating and filtering baits **\n"
		write_file(".log.txt", "VariantBaits") if $options.log
		baits = snp_to_baits(@selectsnps, refseq)
		if $options.filter
			for i in 0 ... baits.size # Filtered snps
				write_file( "-filtered.vcf", vcffilt)
				for snp in baits[baits.keys[i]]
					write_file( "-filtered.vcf", snp.line)
				end
			end
		end
	end
end
