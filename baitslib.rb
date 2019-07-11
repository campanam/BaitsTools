#!/usr/bin/env ruby
#-----------------------------------------------------------------------------------------------
# baitslib
BAITSLIBVER = "1.3.3"
# Michael G. Campana, 2017-2019
# Smithsonian Conservation Biology Institute
#-----------------------------------------------------------------------------------------------

class Fa_Seq #container for fasta/fastq sexquences
	attr_accessor :header, :circular, :fasta, :seq, :qual, :qual_array, :bedstart, :locus, :bedheader
	def initialize(header = "", circular = false, fasta = false, seq = "", qual = "", bedstart = 0, locus = "")
		@header = header # Sequence header
		@circular = circular # Circular sequence flag
		@seq = seq # DNA sequence
		@qual = qual # Original quality scores
		@qual_array = [] # Array of numeric quality scores
		@fasta = fasta # FASTA format flag
		@bedstart = bedstart # Sequence start in absolute BED coordinates
		@locus = locus # locus id for alignment
		@bedheader = header.dup # Default is header, change for region extraction
	end
	def make_dna # Replace uracils with thymines for internal consistency
		@seq.gsub!("u", "t")
		@seq.gsub!("U", "T")
	end
	def calc_quality # Convert quality scores to numeric values so only needed once
		for i in 0...@qual.length
			@qual_array.push($fq_hash[@qual[i]])
		end
	end
	def numeric_quality
		return @qual_array
	end
end
#-----------------------------------------------------------------------------------------------
class Chromo_SNP # Container for SNPs
	attr_accessor :chromo, :snp, :popvar_data, :ref, :alt, :qual, :category, :popcategory
	def initialize(chromo, snp, popvar_data =[], ref = nil, alt = nil, qual = nil, line = nil)
		@chromo = chromo # Chromosome
		@snp = snp # SNP index in 1-based indexing
		@popvar_data = popvar_data # Array holding stacks population-specific allele data
		@ref = ref # Reference allele (vcf)
		@alt = alt # Alternate alleles (vcf)
		@qual = qual # Allele quality for vcf
		@line = line # Original datalines for recall
		@category = nil # Sorting category for vcf variants
		@popcategory = {} # Sorting category for vcf variants within populations
	end
	def categorize # Categorize VCF variants
		alleles = line.split("\t")[9..-1]
		taxa = $options.taxa.uniq
		taxa_hash = {}
		for taxon in taxa
			taxa_hash[taxon] = []
		end
		for i in 0 ... alleles.size
			taxa_hash[$options.taxa[i]].push(alleles[i][0]) unless alleles[i][0] == "."
			taxa_hash[$options.taxa[i]].push(alleles[i][2]) unless alleles[i][2] == "."
			unless alleles[i][0] == "." or alleles[i][2] == "." # Get population-level variation
				@popcategory[$options.taxa[i]] = true if alleles[i][0] != alleles[i][2]
			end
			taxa_hash[$options.taxa[i]] = taxa_hash[$options.taxa[i]].uniq 
			taxa_hash[$options.taxa[i]].sort! unless taxa_hash[$options.taxa[i]].nil?
		end
		observed_combinations = {}
		for taxon in taxa_hash
			if observed_combinations[taxon[1]].nil?
				observed_combinations[taxon[1]] = 1
			else
				observed_combinations[taxon[1]] += 1
			end
		end
		observed_combinations = observed_combinations.keys
		if observed_combinations.include?([])
			@category = "MissingData"
		elsif observed_combinations.size == 1
			if observed_combinations[0].size == 1
				@category = "Invariant"
			else
				@category = "AllPopulations"
			end
		else
			homozygous = true
			for i in 0 ... observed_combinations.size
				if observed_combinations[i].size > 1
					homozygous = false
				end
			end
			if homozygous
				@category = "BetweenPopulations"
			else
				@category = "WithinPopulations"
			end
		end
		write_file(".log.txt", @chromo + "\t" + @snp.to_s + "\t" + @category) if $options.log
	end
	def within_pops? # Determine whether stacks SNP variable within populations
		within_pops = false
		for pop in @popvar_data
			if !pop.monomorphic?
				within_pops = true
				break
			end
		end
		return within_pops
	end
	def line
		if @line.nil? # Combine data from stacks2baits Popvars if line not previously defined
			@line = ""
			for pop in @popvar_data
				@line << pop.line
			end
		end
		return @line[0..-2] # Remove final line break
	end
	def alt_alleles
		if @alt.nil? # Get alternate alleles for stacks2baits (once)
			@alt = []
			for pop in @popvar_data
				@alt += pop.alleles
			end
			@alt = @alt.uniq
		end
		return @alt
	end
end
#-----------------------------------------------------------------------------------------------
def checkpop # Method to determine whether popcategories input is reasonable
	popbad = false
	for pop in $options.popcategories
		if pop < 0 or pop > $options.taxacount[2]
			popbad = true
			break
		end
	end
	return popbad
end
#-----------------------------------------------------------------------------------------------
def mean(val = [])
	mean = 0.0
	for i in 0.. val.size-1
		mean += val[i]
	end
	mean /= val.size
	return mean
end
#-----------------------------------------------------------------------------------------------
def get_ambiguity(base)
	vars = []
	case base
	when "A","a","T","t","G","g","C","c","-"
		vars.push(base)
	when "R"
		vars.push("A","G")
	when "Y"
		vars.push("C","T")
	when "M"
		vars.push("A","C")
	when "K"
		vars.push("G","T")
	when "S"
		vars.push("C","G")
	when "W"
		vars.push("A","T")
	when "H"
		vars.push("A","C","T")
	when "B"
		vars.push("C","G","T")
	when "V"
		vars.push("A","C","G")
	when "D"
		vars.push("A","G","T")
	when "N"
		vars.push("A","C","G","T")
	when "r"
		vars.push("a","g")
	when "y"
		vars.push("c","t")
	when "m"
		vars.push("a","c")
	when "k"
		vars.push("g","t")
	when "s"
		vars.push("c","g")
	when "w"
		vars.push("a","t")
	when "h"
		vars.push("a","c","t")
	when "b"
		vars.push("c","g","t")
	when "v"
		vars.push("a","c","g")
	when "d"
		vars.push("a","g","t")
	when "n"
		vars.push("a","c","g","t")
	end
	return vars
end
#-----------------------------------------------------------------------------------------------
def collapse_ambiguity(bait, force = false) # Ambiguity handling, force turns off no_Ns
	r = ["A","G"]
	rl = ["a","g"]
	y = ["C","T"]
	yl = ["c","t"]
	m = ["A","C"]
	ml = ["a","c"]
	k = ["G","T"]
	kl = ["g","t"]
	s = ["C","G"]
	sl = ["c","g"]
	w = ["A","T"]
	wl = ["a","t"]
	h = ["A","C","T"]
	hl = ["a","c","t"]
	b = ["C","G","T"]
	bl = ["c","g","t"]
	v = ["A","C","G"]
	vl = ["a","c","g"]
	d = ["A","G","T"]
	dl = ["a","g","t"]
	n = ["A","C","G","T"]
	nl = ["a","c","g","t"]
	bait.gsub!("R",r[rand(2)])
	bait.gsub!("Y",y[rand(2)])
	bait.gsub!("M",m[rand(2)])
	bait.gsub!("K",k[rand(2)])
	bait.gsub!("S",s[rand(2)])
	bait.gsub!("W",w[rand(2)])
	bait.gsub!("H",h[rand(3)])
	bait.gsub!("B",b[rand(3)])
	bait.gsub!("V",v[rand(3)])
	bait.gsub!("D",d[rand(3)])
	bait.gsub!("N",n[rand(4)]) unless ($options.no_Ns && !force)
	bait.gsub!("T","U") if $options.rna # Remove any residual Ts after ambiguity collapsing
	bait.gsub!("r",rl[rand(2)])
	bait.gsub!("y",yl[rand(2)])
	bait.gsub!("m",ml[rand(2)])
	bait.gsub!("k",kl[rand(2)])
	bait.gsub!("s",sl[rand(2)])
	bait.gsub!("w",wl[rand(2)])
	bait.gsub!("h",hl[rand(3)])
	bait.gsub!("b",bl[rand(3)])
	bait.gsub!("v",vl[rand(3)])
	bait.gsub!("d",dl[rand(3)])
	bait.gsub!("n",nl[rand(4)]) unless ($options.no_Ns && !force)
	bait.gsub!("t","u") if $options.rna # Remove any residual Ts after ambiguity collapsing
	return bait
end
#-----------------------------------------------------------------------------------------------
def build_fq_hash # method to build fastq quality hash
	fq_val = ["!","\"","#","$","%","&","\'","(",")","*","+",",","-",".","/","0","1","2","3","4","5","6","7","8","9",
		":",";","<","=",">","?","@","A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T",
		"U","V","W","X","Y","Z","[","\\","]","^","_","`","a","b","c","d","e","f","g","h","i","j","k","l","m","n","o",
		"p","q","r","s","t","u","v","w","x","y","z","{","|","}","~"]
	$fq_hash = {}
	if $options.phred64
		for i in 0 .. 40
			$fq_hash[fq_val[i]] = i + 31
		end
	else
		for i in 0 .. 93
			$fq_hash[fq_val[i]] = i
		end
	end
end
#-----------------------------------------------------------------------------------------------
def extend_baits(bait, reference, seqst, seqend) # Extend baits with gap characters
	bait.delete!("-")
	while bait.length < $options.baitlength
		seqend += 1
		break if seqend >= reference.length
		bait += reference[seqend].chr
		bait.delete!("-")
	end
	while bait.length < $options.baitlength
		seqst -=1
		break if seqst < 0
		bait.insert(0,reference[seqst].chr)
		bait.delete!("-")
	end
	return bait
end
#-----------------------------------------------------------------------------------------------
def reversecomp(bait)
	revcomp = ""
	rchash = {"A" => "T","G" => "C", "Y" => "R", "S" => "S", "W" => "W", "K" => "M", "N" => "N", "B" => "V", "D" => "H", "-" => "-"}
	rchash.merge!(rchash.invert)
	rchash.merge!({ "U" => "A"})
	rclc = {}
	for key in rchash.keys
		rclc[key.downcase] = rchash[key].downcase
	end
	rchash.merge!(rclc)
	for base in 0 ... bait.length
		revcomp += rchash[bait[base]]
	end
	revcomp.gsub!("T", "U") if $options.rna
	revcomp.gsub!("t", "u") if $options.rna
	revcomp.reverse!
	return revcomp
end
#-----------------------------------------------------------------------------------------------
def linguistic_complexity(bait)
	testbait = collapse_ambiguity(bait.dup, true) # Don't mess up original bait sequence, force N collapsing
	testbait.delete("-").upcase!
	ngram = 0
	maxngram = 0
	for i in 1 .. testbait.length
		testbait.length - i + 1 < 4 ** i ? maxigram = testbait.length - i + 1 : maxigram = 4 ** i
		ngrams = []
		for j in 0 .. testbait.length - i
			ngr = testbait[j..j+i -1]
			ngrams.push(ngr) unless ngrams.include?(ngr)
			break if ngrams.size == maxigram
		end
		ngram += ngrams.size
		maxngram += maxigram
	end
	return (ngram.to_f/maxngram.to_f)
end
#-----------------------------------------------------------------------------------------------
def max_homopolymer(testbait)
	bases = {"A" => ["A","R","W","M","D","H","V","N"], "T" => ["T","U","Y","W","K","B","D","H","N"], "U" => ["T","U","Y","W","K","B","D","H","N"], "G" => ["G","R","S","K","B","D","V","N"], "C" => ["C","Y","S","M","B","H","V","N"], "" => []}
	homopoly = 1
	max_homopoly = 1
	bait = testbait.dup.upcase # Don't mess up original bait sequence
	bait.gsub!("U","T")
	bait.delete!("-")
	base = ""
	for i in 0 ... bait.length
		if bases[base].include?(bait[i])
			homopoly += 1
		else
			posbases = get_ambiguity(bait[i])
			for j in i + 1 ... bait.length # Check if ambiguities allow a possible forward homopolymer
				nextposbases = get_ambiguity(bait[j])
				for posbase in posbases
					posbases.delete(posbase) unless nextposbases.include?(posbase)
				end
				if posbases.size == 0
					base = collapse_ambiguity(bait[i])
					break
				elsif posbases.size == 1
					base = posbases[0]
					break
				end
			end
			max_homopoly = homopoly if homopoly > max_homopoly
			homopoly = 1
		end
	end
	max_homopoly = homopoly if homopoly > max_homopoly # Allow update if last string max homopolymer
	return max_homopoly
end
#-----------------------------------------------------------------------------------------------
def filter_baits(bait, qual = [0])
	# To be implemented: Self-complementarity filter
	bait = collapse_ambiguity(bait) if $options.collapse_ambiguities # Ambiguity handling
	keep = true
	keep = false if (bait.length < $options.baitlength && $options.completebait)
	keep = false if (bait.include?("N") and $options.no_Ns)
	keep = false if (bait.include?("-") and $options.gaps == "exclude")
	gc = 0.0
	gc += bait.count("GCSgcs").to_f
	gc += 0.5 * bait.count("RYKMNrykmn").to_f
	gc += 0.667 * bait.count("BVbv").to_f
	gc += 0.333 * bait.count("DHdh").to_f
	mask = bait.count("atgcurymkswhbvdn").to_f
	gccont = gc/bait.length.to_f
	maskcont = mask/bait.length.to_f
	keep = false if (gccont * 100.0 < $options.mingc && $options.mingc_filter)
	keep = false if (gccont * 100.0 > $options.maxgc && $options.maxgc_filter)
	keep = false if (maskcont * 100.0 > $options.maxmask && $options.maxmask_filter)
	meanqual = melt = minqual = maxhomopoly = seqcomp = "" # Default values if not calculated
	if $options.mint_filter or $options.maxt_filter or $options.params
		case $options.bait_type
		when "RNA-RNA"
			melt = 79.8 + (58.4 * gccont) + (11.8 * (gccont**2.0)) - 820.0/bait.length.to_f + 18.5 * Math.log($options.na) - 0.35 * $options.formamide
		when "RNA-DNA"
			melt = 79.8 + (58.4 * gccont) + (11.8 * (gccont**2.0)) - 820.0/bait.length.to_f + 18.5 * Math.log($options.na) - 0.50 * $options.formamide
		when "DNA-DNA"
			melt = 81.5 + (41.0 * gccont) - 500.0/bait.length.to_f + 16.6 * Math.log($options.na) - 0.62 * $options.formamide
		end
		keep = false if (melt < $options.mint && $options.mint_filter)
		keep = false if (melt > $options.maxt && $options.maxt_filter)
	end
	if $options.meanqual_filter or $options.params
		meanqual = mean(qual)
		keep = false if (meanqual < $options.meanqual && $options.meanqual_filter)
	end
	if $options.minqual_filter or $options.params
		minqual = qual.min
		keep = false if (minqual < $options.minqual && $options.minqual_filter)
	end
	if $options.maxhomopoly_filter or $options.params
		maxhomopoly = max_homopolymer(bait)
		keep = false if (maxhomopoly > $options.maxhomopoly && $options.maxhomopoly_filter)
	end
	if $options.lc_filter or $options.params
		if $options.no_lc && !$options.lc_filter
			seqcomp = "NA"
		else
			seqcomp = linguistic_complexity(bait)
			keep = false if (seqcomp < $options.lc && $options.lc_filter)
		end
	end
	return [keep, bait.length.to_s + "\t" + (gccont * 100.0).to_s + "\t" + melt.to_s + "\t" + maskcont.to_s + "\t" +  maxhomopoly.to_s + "\t" + seqcomp.to_s + "\t" + meanqual.to_s + "\t" + minqual.to_s + "\t" + bait.include?("N").to_s + "\t" + bait.include?("-").to_s + "\t" + keep.to_s]
end
#-----------------------------------------------------------------------------------------------
def setup_output
	$options.filestem = $options.outdir + "/" + $options.outprefix
	$options.default_files.push(".log.txt") if $options.log
	unless $options.no_baits
		$options.default_files.push("-baits.fa") unless $options.algorithm == "checkbaits"
		$options.default_files.push("-baits.bed") if $options.coords
		$options.default_files.push("-baits-relative.bed") if $options.rbed
	end
	if $options.filter
		$options.default_files.push("-filtered-baits.fa")
		$options.default_files.push("-filtered-params.txt") if $options.params
		$options.default_files.push("-filtered-baits.bed") if $options.coords
		$options.default_files.push("-filtered-baits-relative.bed") if $options.rbed
	end
	if $options.algorithm == "stacks2baits"
		$options.default_files.push("-betweenpops.tsv")
		if $options.sort and $options.hwe
			$options.default_files.push("-inhwe.tsv")
			$options.default_files.push("-outhwe.tsv")	
		elsif $options.sort
			$options.default_files.push("-withinpops.tsv")
		end
	end
end
#-----------------------------------------------------------------------------------------------
def setup_temp(datasize)
	$options.threads > datasize ? $options.used_threads = datasize : $options.used_threads = $options.threads
	splits = []
	val = 0
	crit = (datasize.to_f/$options.used_threads.to_f).ceil
	while val < datasize
		splits.push(val)
		val += crit
	end
	splits.push(datasize)
	return splits
end
#-----------------------------------------------------------------------------------------------
def write_file(fileid, message, tmp = false, tmpfile = 0)
	(tmp && tmpfile != 0) ? file = $options.filestem + ".tmp" + tmpfile.to_s + fileid : file = $options.filestem + fileid 
	File.open(file, 'a') do |write|
		write << message + "\n"
	end
end
#-----------------------------------------------------------------------------------------------
def resolve_unix_path(path)
	reserved = ["\\", ";", "&", "(", ")","*","?","[","]","~",">","<","!","\"","\'", "$", " "] # \ first to prevent repeated gsub issues
	outpath = path.dup # Avoid rewriting original string
	for reschar in reserved
		outpath.gsub!(reschar) {"\\" + reschar} # Use odd syntax because escape character causes issues with backslash in sub
	end
	return outpath
end
#-----------------------------------------------------------------------------------------------
def cat_files(files = $options.default_files)
	for file in files
		for i in 1 ... $options.used_threads
			tmpfile = $options.filestem + '.tmp' + i.to_s + file
			if File.exist?(tmpfile)
				system("cat #{resolve_unix_path(tmpfile)} >> #{resolve_unix_path($options.filestem + file)}")
				system("rm", tmpfile)
			end
		end
	end
end
#-----------------------------------------------------------------------------------------------
def read_fasta(file) # Read fasta and fastq files
	seq_array = []
	faseq = nil # Dummy value
	qual = false # Denotes whether sequence or quality string
	File.open(file, 'r') do |seq|
		while line = seq.gets
			unless line == "\n" # Remove extraneous line breaks
				if line[0].chr == ">" and qual == false
					if !faseq.nil? # push previously completed sequence into array
						seq_array.push(faseq)
					end
					header = line[1...-1] # Remove final line break and beginning >
					header_array = header.split("#")
					if header_array.include?("circ")  # Determine circularity
						circular = true
						header_array.delete("circ") # Delete array value so bed value always at same index
					else
						circular = false
					end
					# Get locus information
					locus = header_array.find { |gene| /^loc/ =~ gene }
					header_array.delete(locus) unless locus.nil?
					if header_array[1].nil?
						bedstart = 0 # Set default bedstart to 0
					else
						bedstart = header_array[1][3..-1].to_i
					end
					header_array = header_array[0].split(" ") if $options.ncbi # Remove sequence description from sequence name for NCBI-formatted files
					header = header_array[0]
					faseq = Fa_Seq.new(header, circular, true)
					faseq.bedstart = bedstart
					faseq.locus = locus[3..-1] unless locus.nil?
				elsif qual and faseq.qual.length < faseq.seq.length # Test if quality line is complete
					faseq.qual << line[0...-1] #Exclude final line break, allow multi-line fastq
			 	elsif line[0].chr == "@" # Start new sequence
					qual = false
					if !faseq.nil? # push previously completed sequence into array
						seq_array.push(faseq)
					end
					header = line[1...-1] # Remove final line break and beginning >
					header_array = header.split("#")
					if header_array.include?("circ")  # Determine circularity
						circular = true
						header_array.delete!("circ") # Delete array value so bed value always at same index
					else
						circular = false
					end
					# Get locus information
					locus = header_array.find { |gene| /^loc/ =~ gene }
					header_array.delete(locus) unless locus.nil?
					if header_array[1].nil?
						bedstart = 0 # Set default bedstart to 0
					else
						bedstart = header_array[1][3..-1].to_i
					end
					header_array = header_array[0].split(" ") if $options.ncbi # Remove sequence description from sequence name for NCBI-formatted files
					header = header_array[0]
					faseq = Fa_Seq.new(header, circular, false)
					faseq.bedstart = bedstart
					faseq.locus = locus[3..-1] unless locus.nil?
				elsif line[0].chr == "+"
					qual = true
				else
					faseq.seq << line[0...-1] #Exclude final line break, allow multi-line fasta/fastq
				end
			end
		end
	end
	seq_array.push(faseq) # Put last sequence into fasta array
	threads = [] # Array to hold threads
	$options.threads.times do |i|
		threads[i] = Thread.new {
			for Thread.current[:j] in 0 ... seq_array.size
				if Thread.current[:j] % $options.threads == i
					seq_array[Thread.current[:j]].calc_quality
					seq_array[Thread.current[:j]].make_dna
				end
			end
		}
	end
	threads.each { |thr| thr.join }
	if $options.log
		logtext = "SequencesRead\nType\tNumberCircular\tNumberLinear\tTotalNumber\tCircularBp\tLinearBp\tTotalBp\n"
		facirc = 0
		falin = 0
		facircbp = 0
		falinbp = 0
		fqcirc = 0
		fqlin = 0
		fqcircbp = 0
		fqlinbp = 0
		for seq in seq_array
			if seq.fasta && seq.circular
				facirc += 1
				facircbp += seq.seq.length
			elsif seq.fasta && !seq.circular
				falin += 1
				falinbp += seq.seq.length
			elsif !seq.fasta && seq.circular
				fqcirc += 1
				fqcircbp += seq.seq.length
			else
				fqlin += 1
				fqlinbp += seq.seq.length
			end
		end
		logtext << "FASTA\t" + facirc.to_s + "\t" + falin.to_s + "\t" + (facirc + falin).to_s + "\t" + facircbp.to_s + "\t" + falinbp.to_s + "\t" + (facircbp + falinbp).to_s + "\n"
		logtext << "FASTQ\t" + fqcirc.to_s + "\t" + fqlin.to_s + "\t" + (fqcirc + fqlin).to_s + "\t" + fqcircbp.to_s + "\t" + fqlinbp.to_s + "\t" + (fqcircbp + fqlinbp).to_s + "\n"
		logtext << "Total\t" + (facirc + fqcirc).to_s + "\t" + (falin + fqlin).to_s + "\t" + (facirc + falin + fqcirc + fqlin).to_s + "\t" + (facircbp + fqcircbp).to_s + "\t" + (falinbp + fqlinbp).to_s + "\t" + (facircbp + falinbp + fqcircbp + fqlinbp).to_s + "\n"
		write_file(".log.txt", logtext)
	end
	return seq_array
end
#-----------------------------------------------------------------------------------------------
def selectsnps(snp_hash) # Choose SNPs based on input group of SNPSs
	# Sort chromosomal SNPs in case unsorted
	totalvar = 0
	selectvar = 0
	for chromo in snp_hash.keys
		totalvar += snp_hash[chromo].size
		snp_hash[chromo].sort_by! { |snp| snp.snp }
		if $options.taxafile != nil
			snp_hash[chromo].delete_if { |snp| snp.category == "MissingData" || snp.category == "Invariant" }
			snp_hash[chromo].delete_if { |snp| snp.category == "AllPopulations" } if $options.taxacount[0] == 0
			snp_hash[chromo].delete_if { |snp| snp.category == "BetweenPopulations" } if $options.taxacount[1] == 0
			snp_hash[chromo].delete_if { |snp| snp.category == "WithinPopulations"} if $options.taxacount[2] == 0
			unless $options.popcategories.nil? # Remove unwanted population-specific SNPs
				for popcat in $options.popcategories.keys
					snp_hash[chromo].delete_if { |snp| snp.popcategory.include?(popcat) } if $options.popcategories[popcat] == 0
				end
			end
			snp_hash.delete_if {|key, value | key == chromo} if snp_hash[chromo].size == 0
		end
	end
	selectsnps = {}
	all_populations = 0
	between_populations = 0
	within_populations = 0
	unless $options.popcategories.nil? # Set count of selected popcategories to 0
		popcategories = $options.popcategories.dup
		for key in popcategories.keys
			popcategories[key] = 0
		end
	end
	if !$options.every
		for i in 1..$options.totalsnps
			selected_contig = snp_hash.keys[rand(snp_hash.size)] # Get name of contig
			snpindex = rand(snp_hash[selected_contig].size)
			selected_snp = snp_hash[selected_contig][snpindex]
			# Delete SNPs that are too close
			if snpindex + 1 < snp_hash[selected_contig].size
				while (snp_hash[selected_contig][snpindex + 1].snp - selected_snp.snp).abs < $options.distance
					snp_hash[selected_contig].delete(snp_hash[selected_contig][snpindex + 1])
					break if snpindex + 1 >= snp_hash[selected_contig].size
				end
			end
			if snpindex > 0
				while (snp_hash[selected_contig][snpindex - 1].snp - selected_snp.snp).abs < $options.distance
					snp_hash[selected_contig].delete(snp_hash[selected_contig][snpindex - 1])
					snpindex -= 1
					break if snpindex < 1
				end
			end
			# Add SNP to selected pool and delete contigs if maximum number of SNPs reached or no remaining SNPs on contig
			selectsnps[selected_contig] = [] if selectsnps[selected_contig].nil?
			selectsnps[selected_contig].push(selected_snp)
			snp_hash[selected_contig].delete(selected_snp) # So it cannot be reselected
			# Delete remaining SNPs of category if maximum amount reached
			if $options.taxafile != nil
				case selected_snp.category
				when "AllPopulations"
					all_populations += 1
					if all_populations == $options.taxacount[0]
						for chromo in snp_hash.keys
							snp_hash[chromo].delete_if { |snp| snp.category == "AllPopulations" }
							snp_hash.delete_if {|key, value | key == chromo} if snp_hash[chromo].size == 0
						end
					end
				when "BetweenPopulations"
					between_populations += 1
					if between_populations == $options.taxacount[1]
						for chromo in snp_hash.keys
							snp_hash[chromo].delete_if { |snp| snp.category == "BetweenPopulations" }
							snp_hash.delete_if {|key, value | key == chromo} if snp_hash[chromo].size == 0
						end
					end
				when "WithinPopulations"
					within_populations += 1
					unless $options.popcategories.nil? # Remove unwanted population-specific SNPs
						for popcat in $options.popcategories.keys
							popcategories[popcat] += 1 if selected_snp.popcategory.include?(popcat)
							if popcategories[popcat] == $options.popcategories[popcat]
								for chromo in snp_hash.keys
									snp_hash[chromo].delete_if { |snp| snp.popcategory.include?(popcat) }
									snp_hash.delete_if {|key, value | key == chromo} if snp_hash[chromo].size == 0
								end
							end
						end
					end
					if within_populations == $options.taxacount[2]
						for chromo in snp_hash.keys
							snp_hash[chromo].delete_if { |snp| snp.category == "WithinPopulations" }
							snp_hash.delete_if {|key, value | key == chromo} if snp_hash[chromo].size == 0
						end
					end
				end
			end
			if $options.scale
				maxsize = $options.scalehash[selected_contig]
			else
				maxsize = $options.maxsnps
			end
			unless snp_hash[selected_contig].nil?
				if selectsnps[selected_contig].size == maxsize or snp_hash[selected_contig].size == 0
					snp_hash.delete_if {|key, value | key == selected_contig}
				end
			end
			break if snp_hash.size == 0 # Stop selecting SNPs when all contigs deleted from consideration
		end
	else
		selectsnps = snp_hash
	end
	write_file(".log.txt", "Chromosome\tSelectedVariants") if $options.log
	for chromo in selectsnps.keys
		selectsnps[chromo].sort_by! { |snp| snp.snp }
		if $options.log
			selectvar += selectsnps[chromo].size
			logtext = chromo + "\t"
			for snp in selectsnps[chromo]
				logtext << snp.snp.to_s + ","
			end
			write_file(".log.txt",  logtext[0...-1]) # Remove final comma
		end
	end
	if $options.log
		write_file(".log.txt", "\nNumberTotalVariants\tNumberSelectedVariants\n" + totalvar.to_s + "\t" + selectvar.to_s + "\n")
		if $options.taxafile != nil
			write_file(".log.txt", "NumberAllPopulations\tNumberBetweenPopulations\tNumberWithinPopulations\n" + all_populations.to_s + "\t" + between_populations.to_s + "\t" + within_populations.to_s + "\n")
			unless $options.popcategories.nil?
				popcatline = "Population-Specific Variants\n" + $options.taxa.uniq.join("\t") + "\n" + popcategories.values.join("\t") + "\n"
				write_file(".log.txt",popcatline)
			end
		end
	end
	return selectsnps
end
#-----------------------------------------------------------------------------------------------
def snp_to_baits(selectedsnps, refseq, filext = "")
	if $options.params
		paramline = "Chromosome:Coordinates\tSNP\tBaitLength\tGC%\tTm\tMasked%\tMaxHomopolymer\tSeqComplexity\tMeanQuality\tMinQuality\tNs\tGaps\tKept"
		write_file("-filtered-params.txt", paramline)
	end
	filteredsnps = {}
	threads = []
	if $options.log
		logs = []
		refseq.size.times { logs.push([]) }
		logtext = "Chromosome\tVariant\tNumberBaits\tRetainedBaits\tExcludedBaits"
		write_file(".log.txt", logtext)
	end
	@splits = setup_temp(refseq.size)
	$options.threads.times do |i|
		threads[i] = Thread.new {
			for Thread.current[:j] in @splits[i] ... @splits[i+1]
				if selectedsnps.keys.include?(refseq[Thread.current[:j]].header)
					for Thread.current[:snp] in selectedsnps[refseq[Thread.current[:j]].header]
						Thread.current[:totalbaits] = 0 if $options.log
						Thread.current[:retainedbaits] = 0 if $options.log
						Thread.current[:tiling] = $options.tiledepth
						Thread.current[:rseq_vars] = [refseq[Thread.current[:j]]]
						if $options.alt_alleles
							for Thread.current[:altvar] in Thread.current[:snp].alt_alleles
								Thread.current[:tmpseq] = refseq[Thread.current[:j]].seq.dup # Need to duplicate to prevent changes to original sequence
								if Thread.current[:snp].ref.nil?
									Thread.current[:tmpseq][Thread.current[:snp].snp-1]=Thread.current[:altvar]
								else
									Thread.current[:tmpseq][Thread.current[:snp].snp-1..Thread.current[:snp].snp+Thread.current[:snp].ref.length-2]=Thread.current[:altvar]
								end
								Thread.current[:tmpheader] = refseq[Thread.current[:j]].header.dup + "_" + Thread.current[:altvar] + "-variant"
								Thread.current[:tmp] = Fa_Seq.new(Thread.current[:tmpheader], refseq[Thread.current[:j]].circular, refseq[Thread.current[:j]].fasta, Thread.current[:tmpseq], refseq[Thread.current[:j]].qual)
								unless refseq[Thread.current[:j]].fasta
									Thread.current[:qual_arr] = []
									Thread.current[:tmpqual] = refseq[Thread.current[:j]].qual_array.dup
									if $options.algorithm == "vcf2baits"
										Thread.current[:qual] = Thread.current[:snp].qual
										Thread.current[:qual] = 40 if Thread.current[:qual] > 40 # Adjust for quality scores beyond typical sequence capability
									else
										Thread.current[:qual] = refseq[Thread.current[:j]].numeric_quality[Thread.current[:snp].snp-1] # Assume quality is equal across the variants since no quality information
									end
									for Thread.current[:i] in 1..Thread.current[:altvar].length # Apply quality score to all bases in alternate allele
										Thread.current[:qual_arr].push(Thread.current[:qual])
									end
									Thread.current[:tmpqual][Thread.current[:snp].snp-1]=Thread.current[:qual_arr]
									Thread.current[:tmp].qual_array = Thread.current[:tmpqual].flatten
								end
								Thread.current[:rseq_vars].push(Thread.current[:tmp]) # Add alternate variants to the sequence
							end
						end
						for Thread.current[:rseq_var] in Thread.current[:rseq_vars]
							for Thread.current[:tile] in 0...Thread.current[:tiling]
								Thread.current[:be4] = Thread.current[:snp].snp - $options.lenbef + Thread.current[:tile] * $options.tileoffset
								Thread.current[:after] = Thread.current[:snp].snp + $options.baitlength - 1 - $options.lenbef + Thread.current[:tile] * $options.tileoffset
								Thread.current[:be4] = 1 if (Thread.current[:be4] < 1 and !Thread.current[:rseq_var].circular)
								Thread.current[:qual] = [$options.fasta_score]
								if Thread.current[:after] > Thread.current[:rseq_var].seq.length and !Thread.current[:rseq_var].circular
									Thread.current[:after] = Thread.current[:rseq_var].seq.length
									if $options.shuffle
										Thread.current[:be4] = $options.baitlength - Thread.current[:after] - 1
										Thread.current[:be4] = 1 if Thread.current[:be4] < 1
									end
									Thread.current[:prb] = Thread.current[:rseq_var].seq[Thread.current[:be4]-1..Thread.current[:after]-1]  #Correct for 0-based counting
									Thread.current[:qual] = Thread.current[:rseq_var].numeric_quality[Thread.current[:be4]-1..Thread.current[:after]-1] unless Thread.current[:rseq_var].fasta
								elsif Thread.current[:after] > Thread.current[:rseq_var].seq.length and Thread.current[:rseq_var].circular
									Thread.current[:after] -= Thread.current[:rseq_var].seq.length
									if Thread.current[:be4] < 1
										Thread.current[:prb] = Thread.current[:rseq_var].seq[Thread.current[:be4]-1..-1] + Thread.current[:rseq_var].seq[Thread.current[:be4]-1..Thread.current[:rseq_var].seq.length-1]+Thread.current[:rseq_var].seq[0..Thread.current[:after]-1] #Correct for 0-based counting and end of sequence
										Thread.current[:qual] = Thread.current[:rseq_var].numeric_quality[Thread.current[:be4]-1..-1] + Thread.current[:rseq_var].numeric_quality[Thread.current[:be4]-1..Thread.current[:rseq_var].seq.length-1]+Thread.current[:rseq_var].numeric_quality[0..Thread.current[:after]-1] unless Thread.current[:rseq_var].fasta
									else
										Thread.current[:prb] = Thread.current[:rseq_var].seq[Thread.current[:be4]-1..Thread.current[:rseq_var].seq.length-1]+Thread.current[:rseq_var].seq[0..Thread.current[:after]-1]  #Correct for 0-based counting and end of sequence
										Thread.current[:qual] = Thread.current[:rseq_var].numeric_quality[Thread.current[:be4]-1..Thread.current[:rseq_var].seq.length-1]+rseq.numeric_quality[0..Thread.current[:after]-1] unless Thread.current[:rseq_var].fasta
									end
								else
									if Thread.current[:be4] < 1
										Thread.current[:prb] = Thread.current[:rseq_var].seq[Thread.current[:be4]-1..-1] + Thread.current[:rseq_var].seq[0..Thread.current[:after]-1]
										Thread.current[:qual] = Thread.current[:rseq_var].numeric_quality[Thread.current[:be4]-1..-1] + Thread.current[:rseq_var].numeric_quality[0..Thread.current[:after]-1] unless Thread.current[:rseq_var].fasta
									else
										Thread.current[:prb] = Thread.current[:rseq_var].seq[Thread.current[:be4]-1..Thread.current[:after]-1]  #Correct for 0-based counting
										Thread.current[:qual] = Thread.current[:rseq_var].numeric_quality[Thread.current[:be4]-1..Thread.current[:after]-1] unless Thread.current[:rseq_var].fasta
									end
								end
								Thread.current[:prb] = reversecomp(Thread.current[:prb]) if $options.rc # Output reverse complemented baits if requested
								Thread.current[:prb].gsub!("T","U") if $options.rna # RNA output handling
								Thread.current[:prb].gsub!("t","u") if $options.rna # RNA output handling
								Thread.current[:prb] = extend_baits(Thread.current[:prb], Thread.current[:rseq_var].seq, Thread.current[:be4]-1, Thread.current[:after]-1) if $options.gaps == "extend" # Basic gap extension
								Thread.current[:seq] = ">" + Thread.current[:rseq_var].header + "_site" + Thread.current[:snp].snp.to_s + "\n" + Thread.current[:prb]
								Thread.current[:be4] = Thread.current[:rseq_var].seq.length + Thread.current[:be4] if Thread.current[:be4] < 1
								Thread.current[:coord] = Thread.current[:rseq_var].header + "\t" + (Thread.current[:be4]-1).to_s + "\t" + Thread.current[:after].to_s
								write_file(filext + "-baits.fa", Thread.current[:seq], true, i)
								Thread.current[:totalbaits] += 1 if $options.log
								write_file(filext + "-baits.bed", Thread.current[:coord], true, i) if $options.coords
								if $options.filter
									Thread.current[:parameters] = filter_baits(Thread.current[:prb], Thread.current[:qual]) # U should not affect filtration
									if Thread.current[:parameters][0]
										write_file(filext + "-filtered-baits.fa", Thread.current[:seq], true, i)
										write_file(filext + "-filtered-baits.bed", Thread.current[:coord], true, i) if $options.coords
										Thread.current[:retainedbaits] += 1 if $options.log
										if filteredsnps[refseq[Thread.current[:j]].header].nil?
											filteredsnps[refseq[Thread.current[:j]].header] = [Thread.current[:snp]]
										else
											filteredsnps[refseq[Thread.current[:j]].header].push(Thread.current[:snp])
											filteredsnps[refseq[Thread.current[:j]].header] = filteredsnps[refseq[Thread.current[:j]].header].uniq
										end
									end
									if $options.params
										Thread.current[:param] = Thread.current[:rseq_var].header + ":" + Thread.current[:be4].to_s + "-" + Thread.current[:after].to_s + "\t" + Thread.current[:snp].snp.to_s + "\t" + Thread.current[:parameters][1]
										write_file(filext + "-filtered-params.txt", Thread.current[:param], true, i)
									end
								end
							end
						end
						if $options.log
							Thread.current[:vals] = [refseq[Thread.current[:j]].header, Thread.current[:snp].snp, Thread.current[:totalbaits]]
							if $options.filter
								Thread.current[:vals].push(Thread.current[:retainedbaits], Thread.current[:totalbaits] - Thread.current[:retainedbaits])
							else
								Thread.current[:vals].push(Thread.current[:totalbaits], "NA")
							end
							logs[Thread.current[:j]].push(Thread.current[:vals])
						end
					end
				end
			end
		}
	end
	threads.each { |thr| thr.join }
	cat_files([filext + "-baits.fa", filext + "-baits.bed", filext + "-filtered-baits.fa", filext + "-filtered-baits.bed", filext + "-filtered-params.txt"])
	if $options.log
		vlogs = [[],[]]
		for i in 0 ... logs.size
			for l in 0 ... logs[i].size
				write_file(".log.txt", logs[i][l].join("\t") )
				vlogs[0].push(logs[i][l][2])
				vlogs[1].push(logs[i][l][3])
			end
		end
		write_file(".log.txt", "\nTotalBaitCoverage(x)\tFilteredBaitCoverage(x)")
		if $options.filter
			write_file(".log.txt", mean(vlogs[0]).to_s + "\t" + mean(vlogs[1]).to_s)
		else
			write_file(".log.txt", mean(vlogs[0]).to_s + "\tNA")
		end
	end
	return filteredsnps
end
#-----------------------------------------------------------------------------------------------
def tile_regions(regions, totallength)
	if regions.size == 0 # Break out if no regions found
		print "** No matching regions found. Exiting. **\n"
		exit
	else
		#Write fasta sequences from the files
		for reg in regions
			write_file("-regions.fa", ">" + reg.header + "\n" + reg.seq)
		end
		write_file(".log.txt", "\nTotalRegions\tTotalRegionLength\n" + regions.size.to_s + "\t" + totallength.to_s + "\n") if $options.log
		#Generate probes using methods from tilebaits
		tilebaits(regions)
	end
end
#-----------------------------------------------------------------------------------------------
def get_command_line # Get command line for summary output
	# Generate basic command line
	cmdline = $options.algorithm + " -i " + resolve_unix_path($options.infile)
	case $options.algorithm
	when "vcf2baits", "stacks2baits"
		if $options.every
			cmdline << " -e"
			cmdline << " -L" + $options.baitlength.to_s + " -O" + $options.tileoffset.to_s + " -b" + $options.lenbef.to_s + " -k" + $options.tiledepth.to_s
		else
			cmdline << " -t" + $options.totalsnps.to_s if $options.taxafile.nil?
			cmdline << " -d" + $options.distance.to_s
			if $options.scale
				cmdline << " -j"
			else
				cmdline << " -m" + $options.maxsnps.to_s
			end
			if $options.no_baits
				cmdline << " -p"
			else
				cmdline << " -L" + $options.baitlength.to_s + " -O" + $options.tileoffset.to_s + " -b" + $options.lenbef.to_s + " -k" + $options.tiledepth.to_s
			end
		end
		cmdline << " -r " + resolve_unix_path($options.refseq) unless $options.no_baits
		cmdline << " -a" if $options.alt_alleles
		cmdline << " -V" + $options.varqual.to_s if $options.varqual_filter
		if $options.taxafile != nil
			cmdline << " --taxafile " + resolve_unix_path($options.taxafile)
			cmdline << " --taxacount " + $options.taxacount.join(",")
			if $options.popcategories != nil
				if $options.popcategories.is_a?(Hash) # Not converted to hash until later in vcf2baits
					cmdline << " --popcategories " + $options.popcategories.values.join(",")
				else
					cmdline << " --popcategories " + $options.popcategories.join(",")
				end
			end
		end
		cmdline << " -S" if $options.sort
		cmdline << " -H -A" + $options.alpha.to_s if $options.hwe
	else
		if $options.algorithm == "annot2baits" or $options.algorithm == "bed2baits" or $options.algorithm == "blast2baits"
			cmdline << " --list " + $options.list_format if $options.algorithm == "bed2baits"
			cmdline << " -r " + resolve_unix_path($options.refseq) + " -P" + $options.pad.to_s
		end
		cmdline << " -L" + $options.baitlength.to_s
		cmdline << " -O" + $options.tileoffset.to_s unless $options.algorithm == "checkbaits"
		if $options.algorithm == "pyrad2baits"
			cmdline << " -I" + $options.minind.to_s
			cmdline << " -W " + $options.strategy
		end
		if $options.algorithm == "aln2baits" or ($options.algorithm == "pyrad2baits" && $options.strategy == "alignment")
			cmdline << " -H " + $options.haplodef
		elsif $options.algorithm == "annot2baits"
			cmdline << " -U "
			for feature in $options.features
				cmdline << feature + ","
			end
			cmdline = cmdline[0...-1] # Remove final , from feature list
		elsif $options.algorithm == "blast2baits"
			cmdline << " --percid " + $options.percid.to_s + " --blastlen " + $options.blastlen.to_s
			cmdline << " --evalue " + $options.evalue.to_s if $options.evalue_filter
		end
		if $options.algorithm == "pyrad2baits" && $options.strategy != "alignment"
			cmdline << " --uncollapsedref" if $options.uncollapsed_ref
			cmdline << " -t" + $options.totalsnps.to_s + " -m" + $options.maxsnps.to_s + " -d" + $options.distance.to_s + " -k" + $options.tiledepth.to_s
			cmdline << " -a" if $options.alt_alleles
		end
	end
	cmdline << " -o " + resolve_unix_path($options.outprefix)
	cmdline << " -Z " + resolve_unix_path($options.outdir)
	cmdline << " -l" if $options.log
	cmdline << " -B" if $options.coords
	cmdline << " -E" if $options.rbed
	cmdline << " --shuffle" if $options.shuffle
	cmdline << " -D" if $options.ncbi
	cmdline << " -Y" if $options.rna
	cmdline << " -R" if $options.rc
	cmdline << " --phred64" if $options.phred64
	cmdline << " -G " + $options.gaps + " -X" + $options.threads.to_s + " --rng " + $options.rng.to_s
	# Generate filtration options
	fltline = ""
	fltline << " -w" if $options.params
	fltline << " --disable-lc" if $options.no_lc
	fltline << " -c" if $options.completebait
	fltline << " -N" if $options.no_Ns
	fltline << " -C" if $options.collapse_ambiguities
	fltline << " -n" + $options.mingc.to_s if $options.mingc_filter
	fltline << " -x" + $options.maxgc.to_s if $options.maxgc_filter
	fltline << " -q" + $options.mint.to_s if $options.mint_filter
	fltline << " -z" + $options.maxt.to_s if $options.maxt_filter
	if $options.mint_filter or $options.maxt_filter
		fltline << " -T " + $options.bait_type + " -s" + $options.na.to_s + " -f" + $options.formamide.to_s
	end
	fltline << " -K" + $options.maxmask.to_s if $options.maxmask_filter
	fltline << " -J" + $options.maxhomopoly.to_s if $options.maxhomopoly_filter
	fltline << " -y" + $options.lc.to_s if $options.lc_filter
	fltline << " -Q" + $options.meanqual.to_s if $options.meanqual_filter
	fltline << " -M" + $options.minqual.to_s if $options.minqual_filter
	fltline << " -F" + $options.fasta_score.to_s if ($options.meanqual_filter or $options.minqual_filter)
	return [cmdline,fltline]
end
