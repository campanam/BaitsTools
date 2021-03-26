#!/usr/bin/env ruby
#-----------------------------------------------------------------------------------------------
# vcf2baits
VCF2BAITSVER = "1.7.0"
# Michael G. Campana, 2017-2021
# Smithsonian Conservation Biology Institute
#-----------------------------------------------------------------------------------------------

def write_filtered_vcf(baits, vcffilt)
	if $options.filter
		write_file( "-filtered.vcf", vcffilt)
		for i in 0 ... baits.size # Filtered snps
			for snp in baits[baits.keys[i]]
				write_file( "-filtered.vcf", snp.line)
			end
		end
		if $options.gzip
			system("gzip #{resolve_unix_path($options.filestem + $options.infix + '-filtered.vcf')}")
		end
	end
end

#-----------------------------------------------------------------------------------------------
def vcf2baits
	# Read VCF file
	@snps = {}
	chromosome = ""
	vcfout = ""
	temp_snps = []
	columns = ""
	print "** Reading VCF **\n"
	gz_file_open($options.infile) do |snpreg|
		while line = snpreg.gets
			vcfout << line if line[0..1] == "##" # Add info lines to header
			if line[0..1] == "#C" # Get column header separately so I can insert new information
				columns = line[0..-2] # Remove line break or extraneous line is added by write_file
				@samples = line[0..-2].split("\t")[9..-1] # Get sample names
			end
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
					temp_snps.push(Chromo_SNP.new(reg, snp, [], ref, alt, qual, line[0..-2])) #Now using Chromo_SNPs to be consistent
				end
			end
		end
	end
	@snps[chromosome]=temp_snps 		#To add last state from last line since otherwise will not loop
	@snps.delete_if {|key, value| key == ""} # Delete dummy original value
	@snps.delete_if {|key, value| value == []} # Delete empty contigs due to QUAL filter
	unless $options.taxafile.nil?
		print "** Reading taxa file and sorting variants **\n"
		gz_file_open($options.taxafile) do |taxafile| # Get taxa assignments
			while line = taxafile.gets
				taxon = line[0..-2].split("\t")
				$options.taxa[taxon[0]] = taxon[1]
			end
		end
		tmptax = [] # Sort taxa by sample order
		for samp in @samples
			tmptax.push($options.taxa[samp])
		end
		$options.taxa = tmptax.dup
		unless $options.popcategories.nil? # Build popcategories hash for snp selection
			tmpcat = $options.popcategories
			$options.popcategories = {}
			tmptax = tmptax.uniq # Using uniq! on a uniqued array returns nil
			for taxon in 0 ... tmptax.size
				$options.popcategories[tmptax[taxon]] = tmpcat[taxon]
			end
		end
		write_file(".log.txt", "Variant Categories\nChromosome\tVariant\tCategory") if $options.log
		@snps.keys.each { |chromo| @snps[chromo].each { |snp| snp.categorize } }
		write_file(".log.txt", "") if $options.log # Add break line to make log easier to read
	end
	print "** Selecting variants **\n" unless $options.every
	write_file(".log.txt", "Variants") if $options.log
	@selectsnps = selectsnps(@snps) # Select SNPs
	# Write VCF & baits
	vcfout << "##baitstools_vcf2baitsVersion=" + VCF2BAITSVER + "+BaitsToolsVersion=" + BAITSTOOLSVER + "+baitslibVersion=" + BAITSLIBVER + "\n##baitstools_vcf2baitsCommand="
	cmdline = get_command_line
	vcffilt = vcfout + cmdline[0] + cmdline[1] + "\n" + columns
	if !$options.every
		write_file("-selected.vcf", vcfout + cmdline[0] + "\n" + columns )
		for i in 0 ... @selectsnps.size
			for snp in @selectsnps[@selectsnps.keys[i]]
				write_file("-selected.vcf", snp.line)
			end
		end
		if $options.gzip
			system("gzip #{resolve_unix_path($options.filestem + '-selected.vcf')}")
		end
	end
	# Output bait sequences unless -p
	if !$options.no_baits
		print "** Reading reference sequence **\n"
		refseq = read_fasta($options.refseq)
		write_file(".log.txt", "VariantBaits") if $options.log
		unless $options.altbaits.nil? # Tile baits under alternate lengths
			baits = snp_to_baits(@selectsnps, refseq, $options.baitlength)
			write_filtered_vcf(baits, vcffilt)
			for altbait in $options.altbaits
				$options.baitlength = altbait
				write_file(".log.txt", "") if $options.log # Add a linebreak between subsequent entries
				baits = snp_to_baits(@selectsnps, refseq, $options.baitlength)
				write_filtered_vcf(baits, vcffilt)
			end
		else
			baits = snp_to_baits(@selectsnps, refseq) # Baits without multibaits infix
			write_filtered_vcf(baits, vcffilt)
		end
		
	end
end
