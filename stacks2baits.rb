#!/usr/bin/env ruby
#-----------------------------------------------------------------------------------------------
# stacks2baits
STACKS2BAITSVER = "1.7.0"
# Michael G. Campana, 2017-2021
# Smithsonian Conservation Biology Institute
#-----------------------------------------------------------------------------------------------

class Popvar # Population-specific SNP data object
	attr_accessor :pop, :alleles, :no_ind, :pfreq, :hetobs, :line 
	def initialize(pop, alleles, no_ind, pfreq, hetobs, line)
		@pop = pop # Population
		@alleles = alleles # Array of alleles
		@no_ind = no_ind # Sample size,
		@pfreq = pfreq # Major frequency
		@hetobs = hetobs # Observed heterozygosity
		@line = line # Original SNP descriptor line
	end
	def monomorphic? # Determine if SNP is monomorphic within population
		if (@pfreq == 1.0 or @pfreq == 0.0)
			return true
		else
			return false
		end
	end
	def in_hwe? # Return whether variant is in HWE for a population
		qfreq = 1.0 - @pfreq # Calculate minor allele frequency
		p2exp = @pfreq ** 2.0 * @no_ind  # Calculate expected major allele homozygotes
		q2exp = qfreq ** 2.0 * @no_ind # Calculate expected minor allele homozygotes
		pqexp = (2.0 * @pfreq * qfreq) * @no_ind # Calculate expected heterozygotes
		pqobs = @hetobs * @no_ind
		p2obs = (@pfreq - @hetobs/2.0) * @no_ind # Calculate observed major alleles in homozygotes
		q2obs = (qfreq - @hetobs/2.0) * @no_ind # Calculate observed minor alleles in homozygotes
		hwe = ((p2obs - p2exp) ** 2.0)/p2exp + ((pqobs - pqexp) ** 2.0)/pqexp + ((q2obs - q2exp) ** 2.0)/q2exp # Calculate chi-square statistic
		if 1.0 - chi_cum_prob(hwe) < $options.alpha # Compare to alpha 0.05
			return true
		else
			return false
		end
	end
end
#-----------------------------------------------------------------------------------------------
def chi_cum_prob(test) # This calculates chi distribution cumulative probability distribution under special case of df = 1
	return Math.erf(Math.sqrt(test/2.0))
end
#-----------------------------------------------------------------------------------------------
def write_stacks(header, snps, tag) # Method to write stacks output since repeating over and over
	write_file(tag + ".tsv", header)
	for key in snps.keys
		for ssnp in snps[key]
			write_file(tag + ".tsv", ssnp.line)
		end
	end
end
#-----------------------------------------------------------------------------------------------
def stacks_altbaits(stacksheader, snpset, refseq, infix, logheader = "") # Reduce redundant programming for multiple baits
	write_file(".log.txt", logheader) if $options.log
	unless $options.altbaits.nil?
		baits = snp_to_baits(snpset, refseq, $options.baitlength, infix)
		write_stacks(stacksheader, baits, $options.baitlength.to_s + "-" + infix + "-filtered") if $options.filter
		for altbait in $options.altbaits
			$options.baitlength = altbait
			write_file(".log.txt", "") if $options.log # Add a linebreak between subsequent entries
			baits = snp_to_baits(snpset, refseq, altbaits, infix)
			write_stacks(stacksheader, baits, altbait.to_s + "-" + infix + "-filtered") if $options.filter
		end
	else
		baits = snp_to_baits(snpset, refseq, nil, infix)
		write_stacks(stacksheader, baits, infix + "-filtered") if $options.filter
	end
end
#-----------------------------------------------------------------------------------------------
def stacks2baits
	# Read stacks summary tsv file
	print "** Reading stacks tsv **\n"
	stacksvars = {} # Hash, keying by stacks Locus ID and SNP index
	stacksheader = "" # Stacks TSV header
	gz_file_open($options.infile) do |stacks|
		while line = stacks.gets
			if line[0].chr != "#"
				split_line = line.split("\t")
				locus = split_line[1]+split_line[4]
				chromo = split_line[2]
				len = split_line[3].to_i # Length
				snp = split_line[4].to_i + 1 #Stacks uses 0-based counting
				pop = split_line[5] # Population
				alleles = [split_line[6], split_line[7]] # Get major, minor alleles
				alleles.delete("-") # Remove non-alleles, separate command to avoid assigning alleles as "-"
				no_ind = split_line[8].to_f # No. of individuals
				pfreq = split_line[9].to_f # Major allele frequency
				hetobs = split_line[10].to_f # Observed heterozygosity
				if stacksvars.include?(locus)
					stacksvars[locus].popvar_data.push(Popvar.new(pop, alleles, no_ind, pfreq, hetobs, line))
				else
					stacksvars[locus]=Chromo_SNP.new(chromo, snp, [Popvar.new(pop, alleles, no_ind, pfreq, hetobs, line)])
					scaled = (len/$options.distance).floor
					$options.scalehash[chromo] = scaled
				end
			else
				stacksheader << line
			end
		end
	end
	# Sort SNPs and convert to usable form for selectsnps algorithm
	print "** Sorting SNPs **\n"
	between_pops = {} # Hash to hold SNPs that are only variable between populations (also all SNPs if not sorting)
	within_pops = {} # Hash to hold SNPs that are variable within populations (overrides between_pops)
	in_hwe = {} # Hash to hold SNPs that are in HWE within populations
	out_hwe = {} # Hash to hold SNPs that are not in HWE within populations
	for key in stacksvars.keys
		snp = stacksvars[key]
		if snp.within_pops? and $options.sort
			if within_pops.include?(snp.chromo)
				within_pops[snp.chromo].push(snp)
			else
				within_pops[snp.chromo]=[snp]
			end
			if $options.hwe
				hwe = true
				for pop in snp.popvar_data # If any population has the SNP out-of-HWE, exclude it
					if !pop.in_hwe?
						hwe = false
						break
					end
				end
				if hwe
					if in_hwe.include?(snp.chromo)
						in_hwe[snp.chromo].push(snp)
					else
						in_hwe[snp.chromo]=[snp]
					end
				else
					if out_hwe.include?(snp.chromo)
						out_hwe[snp.chromo].push(snp)
					else
						out_hwe[snp.chromo]=[snp]
					end
				end
			end
		else
			if between_pops.include?(snp.chromo)
				between_pops[snp.chromo].push(snp)
			else
				between_pops[snp.chromo]=[snp]
			end
		end
	end
	# Select SNPs -- Note that there is no cross-referencing between types
	print "** Selecting SNPs **\n"
	if $options.log
		if $options.sort
			write_file(".log.txt", "BetweenPopsVariants")
		else
			write_file(".log.txt", "AllVariants")
		end
	end
	selected_between = selectsnps(between_pops)
	if $options.sort
		write_stacks(stacksheader, selected_between, "-betweenpops")
	else
		write_stacks(stacksheader, selected_between, "-all")
	end
	if $options.sort and $options.hwe
		write_file(".log.txt", "InHWEVariants") if $options.log
		selected_inhwe = selectsnps(in_hwe)
		write_file(".log.txt", "OutHWEVariants") if $options.log
		selected_outhwe = selectsnps(out_hwe)
		write_stacks(stacksheader, selected_inhwe, "-inhwe")
		write_stacks(stacksheader, selected_outhwe, "-outhwe")	
	elsif $options.sort
		write_file(".log.txt", "WithinPopsVariants") if $options.log
		selected_within = selectsnps(within_pops)
		write_stacks(stacksheader, selected_within, "-withinpops")
	end
	# Output baits unless -p
	if !$options.no_baits
		print "** Reading reference sequence **\n"
		refseq = read_fasta($options.refseq)
		if $options.sort
			stacks_altbaits(stacksheader, selected_between, refseq, "-betweenpops", "BetweenPopsVariantBaits")
		else
			stacks_altbaits(stacksheader, selected_between, refseq, "-all", "AllVariantBaits")
		end
		if $options.sort and $options.hwe
			stacks_altbaits(stacksheader, selected_inhwe, refseq, "-inhwe", "InHWEVariantBaits")
			stacks_altbaits(stacksheader, selected_outhwe, refseq, "-outhwe", "OutHWEVariantBaits")
		elsif $options.sort
			stacks_altbaits(stacksheader, selected_within, refseq, "-withinpops", "WithinPopsVariantBaits")
		end
	end
end
