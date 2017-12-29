#!/usr/bin/env ruby
#-----------------------------------------------------------------------------------------------
# pyrad2baits
PYRAD2BAITSVER = "0.1"
# Michael G. Campana, 2017
# Smithsonian Conservation Biology Institute
#-----------------------------------------------------------------------------------------------

class Pyrad_Locus
	attr_accessor :locus_id, :diplotypes, :snp_sites, :phylo_sites, :refseq
	def initialize(locus_id = 0, diplotypes = [], snp_sites = [], phylo_sites = [], refseq = "")
		@locus_id = locus_id
		@diplotypes = diplotypes # Diplotype sequences
		@snp_sites = snp_sites # All SNP locations, 1-indexed
		@phylo_sites = phylo_sites # Phylogenetically informative SNP locations,1-indexed
		@refseq = refseq
	end
	def get_variants(snp)
		vars = []
		for dip in @diplotypes
			vars.push(get_ambiguity(dip[snp - 1]))
		end
		vars.flatten!.uniq!
	end
	def make_reference # Create a random reference for bait generation
		seq = collapse_ambiguity(@diplotypes[rand(@diplotypes.size)].seq)
		@refseq = Fa_Seq.new(locus_id.to_s, false, true)
		@refseq.seq = seq
	end
end
#-----------------------------------------------------------------------------------------------
def pyrad_to_selectsnps(pyradloci)
	selectsnps = {}
	for locus in pyradloci
		selectsnps[locus.locus_id.to_s] = []
		$options.strategy == "informative" ? snps = locus.phylo_sites : snps = locus.snp_sites
		for snp in snps
			vars = []
			for dip in locus.diplotypes
				vars.push(get_ambiguity(dip.seq[snp - 1]))
			end
			vars.flatten!.uniq!
			ref = locus.refseq.seq[snp - 1]
			vars.delete(ref)
			chrsnp = Chromo_SNP.new(locus.locus_id.to_s, snp, [], ref, vars)
			selectsnps[locus.locus_id.to_s].push(chrsnp)
		end
	end
	selectsnps.delete_if {|key, value| value == []} # Delete empty contigs
	return selectsnps
end
#-----------------------------------------------------------------------------------------------
def pyrad2baits
	print "** Reading LOCI file **\n"
	@pyrad_loci = []
	@pyrad = false # file is in pyrad format
	loc = Pyrad_Locus.new
	locid = 0
	$options.logtext += "Locus Coverage:\nLocus\tNumberIndividuals\tKept\n" if $options.log
	File.open($options.infile, 'r') do |pyr|
		while line = pyr.gets
			if line[0..1] != "//"
				diplo = line.split
				if line[0].chr == ">" # > indicates pyrad format
					@pyrad = true
					diplo[0] = diplo[0][1..-1] # Remove > character
				end
				seq = Fa_Seq.new(diplo[0],false,true)
				seq.seq = diplo[1].gsub("x","n") # Replace x PE separators with n separators
				seq.make_dna
				seq.calc_quality
				loc.diplotypes.push(seq)
			else
				site_info = line.delete("\n").split("|") # Remove line break in case doesn't split
				if @pyrad
					loc.locus_id = locid
					locid += 1
				else
					loc.locus_id = site_info[1].to_i
				end
				for seq in loc.diplotypes
					seq.locus = loc.locus_id.to_s
				end
				sites = site_info[0][diplo[1].length * -1 .. -1]
				for i in 1..sites.length # Get lists of snps and phylogenetically-informative snps (1-indexed)
					if sites[i - 1].chr == "*"
						loc.snp_sites.push(i)
						loc.phylo_sites.push(i)
					elsif sites[i - 1].chr == "-"
						loc.snp_sites.push(i)
					end
				end
				@pyrad_loci.push(loc) unless loc.diplotypes.size < $options.minind # Keep locus unless too few individuals
				if $options.log
					loc.diplotypes.size < $options.minind ? kept = "false\n" : kept = "true\n"
					$options.logtext += loc.locus_id.to_s + "\t" + loc.diplotypes.size.to_s + "\t" + kept
				end
				loc = Pyrad_Locus.new
			end
		end
	end
	$options.logtext += "\n" if $options.log
	if $options.strategy == "alignment"
		# Treat like an alignment
		aln = []
		for loc in @pyrad_loci
			aln.push(loc.diplotypes)
		end
		aln.flatten!
		aln2baits(aln)
	else
		# Treat like a vcf/stacks file.
		print "** Generating reference sequences **\n"
		refseq = []
		refs = ""
		for loc in @pyrad_loci
			loc.make_reference
			refseq.push(loc.refseq)
			refs += ">" + loc.refseq.header + "\n" + loc.refseq.seq + "\n"
		end
		File.open($options.outdir + "/" + $options.outprefix + "-refseq.fa", 'w') do |write|
			write.puts refs
		end
		print "** Selecting variants **\n"
		snp_hash = pyrad_to_selectsnps(@pyrad_loci)
		@selectsnps = selectsnps(snp_hash)
		print "** Generating and filtering baits **\n"
		baits = snp_to_baits(@selectsnps, refseq)
		write_baits(baits[0], baits[1], baits[2], baits[3], baits[4])
	end
end