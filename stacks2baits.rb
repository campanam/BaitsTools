#!/usr/bin/ruby
#-----------------------------------------------------------------------------------------------
# stacks2baits 0.1
# Michael G. Campana, 2016
# Smithsonian Conservation Biology Institute
#-----------------------------------------------------------------------------------------------

class Popvar # Population-specific SNP data object
	attr_accessor :pop, :pnuc, :qnuc, :no_ind, :pfreq, :hetobs # Population, Major allele, minor allele, Major frequency, sample size, observed heterozygosity
	def initialize(pop, pnuc, qnuc, no_ind, pfreq, hetobs)
		@pop = pop
		@pnuc = pnuc
		@qnuc = qnuc
		@no_ind = no_ind
		@pfreq = pfreq
		@hetobs = hetobs
	end
	def monomorphic? # Determine if SNP is monomorphic within population
		if (@pfreq == 1 or @pfreq == 0)
			return true
		else
			return false
		end
	end
	def in_hwe? # Return whether variant is in HWE for a population
		qfreq = 1 - pfreq # Calculate minor allele frequency
		p2exp = pfreq ** 2 * @no_ind  # Calculate expected major allele homozygotes
		q2exp = qfreq ** 2 * @no_ind # Calculate expected minor allele homozygotes
		pqexp = (2 * pfreq * qfreq) * @no_ind # Calculate expected heterozygotes
		pqobs = (@hetobs * @no_ind).to_i
		p2obs = ((pfreq - @hetobs/2) * @no_ind).to_i # Calculate observed major alleles in homozygotes
		q2obs = ((qfreq - @hetobs/2) * @no_ind).to_i # Calculate observed minor alleles in homozygotes
		hwe = ((p2obs - p2exp) ** 2)/p2exp + ((pqobs - pqexp) ** 2)/pqexp + ((q2obs - q2exp) ** 2)/q2exp # Calculate chi-square statistic
		case $options.alpha
		when 0.1
			alpha = 2.706
		when 0.05
			alpha = 3.841
		when 0.025
			alpha = 5.024
		when 0.01
			alpha = 6.635
		end
		if hwe < alpha # Compare to alpha 0.05
			return true
		else
			return false
		end
	end
end
#-----------------------------------------------------------------------------------------------
def stacks2baits
	# Read stacks summary tsv file
	stacksvars = {} # Hash, keying by stacks Locus ID and SNP index
	File.open($options.infile, 'r') do |stacks|
		while line = stacks.gets
			if line[0].chr != "#"
				tabcount = 0
				locus = ""
				chromo = ""
				len = ""
				snp = ""
				pop = ""
				pnuc = ""
				qnuc = ""
				no_ind = ""
				pfreq = ""
				hetobs = ""
				for i in 0..line.length - 2 # Exclude final line break
					if line[i].chr == "\t"
						tabcount += 1
 					else
						case tabcount
						when 1
							locus += line[i].chr
						when 2
							chromo += line[i].chr
						when 3
							len += line[i].chr
						when 4
							snp += line[i].chr #Stacks uses 0-based counting
							locus += line[i].chr
						when 5
							pop += line[i].chr
						when 6
							pnuc += line[i].chr
						when 7
							qnuc += line[i].chr
						when 8
							no_ind += line[i].chr
						when 9
							pfreq += line[i].chr
						when 10
							hetobs += line[i].chr
						when 11
							break # Avoid reading rest of line
						end
					end
				end
				if stacksvars.include?(locus)
					stacksvars[locus].popvar_data.push(Popvar.new(pop, pnuc, qnuc, no_ind.to_i, pfreq.to_f, hetobs.to_f))
				else
					stacksvars[locus]=Chromo_SNP.new(chromo, snp.to_i + 1, [Popvar.new(pop, pnuc, qnuc, no_ind.to_i, pfreq.to_f, hetobs.to_f)])
					scaled = (len.to_i/$options.distance).floor
					$options.scalehash[chromo] = scaled
				end
					
			end
		end
	end
	# Sort SNPs and convert to usable form for selectsnps algorithm
	between_pops = {} # Hash to hold SNPs that are only variable between populations (also all SNPs if not sorting)
	within_pops = {} # Hash to hold SNPs that are variable within populations (overrides between_pops)
	in_hwe = {} # Hash to hold SNPs that are in HWE within populations
	out_hwe = {} # Hash to hold SNPs that are not in HWE within populations
	for key in stacksvars.keys
		snp = stacksvars[key]
		if snp.within_pops? and $options.sort
			if within_pops.include?(snp.chromo)
				within_pops[snp.chromo].push(snp.snp)
			else
				within_pops[snp.chromo]=[snp.snp]
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
						in_hwe[snp.chromo].push(snp.snp)
					else
						in_hwe[snp.chromo]=[snp.snp]
					end
				else
					if out_hwe.include?(snp.chromo)
						out_hwe[snp.chromo].push(snp.snp)
					else
						out_hwe[snp.chromo]=[snp.snp]
					end
				end
			end
		else
			if between_pops.include?(snp.chromo)
				between_pops[snp.chromo].push(snp.snp)
			else
				between_pops[snp.chromo]=[snp.snp]
			end
		end
	end
	# Select SNPs -- Note that there is no cross-referencing between types
	selected_between = selectsnps(between_pops)
	if $options.sort and $options.hwe
		selected_inhwe = selectsnps(in_hwe)
		selected_outhwe = selectsnps(out_hwe)
	elsif $options.sort
		selected_within = selectsnps(within_pops)
	end
	if $options.baits
		bbaits = snp_to_baits(selected_between)
		write_baits(bbaits[0], bbaits[1], bbaits[2], bbaits[3], bbaits[4], $options.infile+"-betweenpops")
		if $options.sort and $options.hwe
			ihbaits = snp_to_baits(selected_inhwe)
			ohbaits = snp_to_baits(selected_outhwe)
			write_baits(ihbaits[0], ihbaits[1], ihbaits[2], ihbaits[3], ihbaits[4], $options.infile+"-inhwe")
			write_baits(ohbaits[0], ohbaits[1], ohbaits[2], ohbaits[3], ohbaits[4], $options.infile+"-outhwe")
		elsif $options.sort
			wbaits = snp_to_baits(selected_within_pops)
			write_baits(wbaits[0], wbaits[1], wbaits[2], wbaits[3], wbaits[4], $options.infile+"-withinpops")
		end
	end
end
