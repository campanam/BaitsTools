#!/usr/bin/ruby
#-----------------------------------------------------------------------------------------------
# select_snps 0.3
# Michael G. Campana, 2015
# Smithsonian Conservation Biology Institute
#-----------------------------------------------------------------------------------------------

::Version = "select_snps 0.3"
require 'optparse'
require 'ostruct'
class Parser
	def self.parse(options)
		# Set defaults
		args = OpenStruct.new
		args.vcfin = ""
		args.totalsnps = 30000
		args.maxsnps = 2
		args.distance = 10000
		args.probes = false
		args.refseq = ""
		args.lenbef = 60
		args.lenaft = 59
		args.filter = false
		args.completeprobe = false
		args.probelength = 120
		args.mingc = 0.0
		args.maxgc = 100.0
		opt_parser = OptionParser.new do |opts|
			opts.banner = "Usage: ruby select_snps-0.3.rb [options]"
			opts.separator ""
			opts.separator "Specific options:"
			opts.separator ""
			opts.on("-i","--input [FILE]", String, "Input VCF file") do |vcf|
				args.vcfin = vcf
			end
			opts.on("-t", "--totalsnps [VALUE]", Integer, "Total requested SNPs (Default = 30,000)") do |tsnps|
				args.totalsnps = tsnps
			end
			opts.on("-m","--maxsnps [VALUE]", Integer, "Maximum number of SNPs per contig (Default = 2)") do |msnps|
				args.maxsnps = msnps
			end
			opts.on("-d","--distance [VALUE]", Integer, "Minimum distance between SNPs within a contig (Default = 10,000)") do |dist|
				args.distance = dist
			end
			opts.on("-p","--probes", "Output probe sequences") do
				args.probes = true
			end
			opts.on("-r","--refseq [FILE]", String, "If probe ouput selected, the reference sequence file") do |ref|
				args.refseq= ref
			end
			opts.on("-b", "--lenbef [VALUE]", Integer, "If probe output selected, the number of bases before the SNP (Default = 60)") do |b4|
				args.lenbef = b4
			end
			opts.on("-a", "--lenaft [VALUE]", Integer, "If probe output selected, the number of bases before the SNP (Default = 59)") do |af|
				args.lenaft = af
			end
			opts.on("-f","--filter", "If probes selected, filter probe sequences") do
				args.filter = true
			end
			opts.on("-c","--complete", "If filtering, require probes be full length") do
				args.completeprobe = true
			end
			opts.on("-n","--mingc [VALUE]", Float, "If filtering, minimum GC content (Default = 0.0)") do |mgc|
				args.mingc = mgc
			end
			opts.on("-g","--maxgc [VALUE]", Float, "If filtering, maximum GC content (Default = 100.0)") do |xgc|
				args.maxgc = xgc
			end
			opts.separator ""
			opts.separator "Common options:"
			opts.on_tail("-h","--help", "Show help") do
				print "Welcome to select_snps.\n\n"
				print "To use the interactive interface, enter <ruby select_snps-0.3.rb>.\n\n"
				puts opts
				exit
			end
			opts.on_tail("-v","--version","Show version") do
				puts ::Version
				exit
			end
		end
		opt_parser.parse!(options)
		return args
	end
end
def filter_probes(probe)
	# Sequence complexity
	# Include/Exclude multi-allelic sites
	# Minimum PIC/Calculate PIC
	# Melting Temperature
	# Self-complementarity
	keep = true
	keep = false if probe.length < $options.probelength
	gc = 0
	for i in 0 ... probe.length
		if probe[i].chr.upcase == "G" or probe[i].chr.upcase == "C"
			gc += 1
		end
	end
	keep = false if (gc * 100/probe.length) < $options.mingc
	keep = false if (gc * 100/probe.length) > $options.maxgc
	return keep
end
begin
	# Data options block
	interact = true if ARGV.size == 0
	$options = Parser.parse(ARGV)
	if interact
		print "Enter input file.\n"
		$options.vcfin = gets.chomp
	end
	while !FileTest.exist?($options.vcfin)
		print "VCF file not found. Please re-enter.\n"
		$options.vcfin = gets.chomp
	end
	if interact
		print "Enter total number of requested SNPs.\n"
		$options.totalsnps = gets.chomp.to_i
	end
	while $options.totalsnps < 1
		print "The total number of SNPs must be greater than 0. Please re-enter.\n"
		$options.totalsnps = gets.chomp.to_i
	end
	if interact
		print "Enter maximum number of SNPs per contig.\n"
		$options.maxsnps = gets.chomp.to_i
	end
	while $options.maxsnps < 1
		print "The maximum number of SNPs per contig must be greater than 0. Please re-enter.\n"
		$options.maxsnps = gets.chomp.to_i
	end
	if interact
		print "Enter minimum distance between SNPs within a contig.\n"
		$options.distance = gets.chomp.to_i
	end
	while $options.distance < 1
		print "The minimum distance between SNPs must be greater than 0. Please re-enter.\n"
		$options.distance = gets.chomp.to_i
	end
	if interact
		print "Output probes? (y/n)\n"
		t = gets.chomp.upcase
		if t == "Y" or t == "YES"
			$options.probes = true
		end
	end
	if $options.probes
		if interact
			print "Enter reference sequence.\n"
			$options.refseq = gets.chomp
		end
		while !FileTest.exist?($options.refseq)
			print "Reference sequence not found. Please re-enter.\n"
			$options.refseq = gets.chomp
		end
		if interact
			print "Enter number of bases before SNP in probe.\n"
			$options.lenbef = gets.chomp.to_i
		end
		while $options.lenbef < 0
			print "The number of probe bases before the SNP must be at least 0. Please re-enter.\n"
			$options.lenbef = gets.chomp.to_i
		end
		if interact
			print "Enter number of bases after SNP in probe.\n"
			$options.lenaft = gets.chomp.to_i
		end
		while $options.lenaft < 0
			print "The number of probe bases after the SNP must be at least 0. Please re-enter.\n"
			$options.lenaft = gets.chomp.to_i
		end
		$options.probelength = $options.lenbef + $options.lenaft + 1
		if interact
			print "Filter probes?\n"
			t = gets.chomp.upcase
			if t == "Y" or t == "YES"
				$options.filter = true
			end
		end
		if $options.filter
			if interact
				print "Require complete length probe?\n"
				t = gets.chomp.upcase
				if t == "Y" or t == "YES"
					$options.completeprobe = true
				end
				print "Enter minimum GC content.\n"
				$options.mingc = gets.chomp.to_f
			end
			while $options.mingc < 0.0 or $options.mingc > 100.0
				print "Minimum GC content must be between 0 and 100. Please re-enter.\n"
				$options.mingc = gets.chomp.to_f
			end
			if interact
				print "Enter maximum GC content.\n"
				$options.maxgc = gets.chomp.to_f
			end
			while $options.maxgc < $options.mingc or $options.maxgc > 100.0
				print "Maximum GC content must be between minimum GC content and 100. Re-enter.\n"
				$options.maxgc = gets.chomp.to_f
			end
		end
	end
	# Read VCF file
	@snps = {}
	chromosome = ""
	temp_snps = []
	File.open($options.vcfin, 'r') do |snpreg|
		while line = snpreg.gets
			if line[0].chr != "#"
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
	@select_snps = {}
	for i in 1..$options.totalsnps
		break if temp_snps.size == 0 # Stop selecting SNPs when all contigs deleted from consideration
		selected_contig = temp_snps.keys[rand(temp_snps.size)] # Get name of contig
		selected_snp = temp_snps[selected_contig][rand(temp_snps[selected_contig].size)]
		#Check SNP/Probe for suitability
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
		@select_snps[selected_contig] = [] if @select_snps[selected_contig].nil?
		@select_snps[selected_contig].push(selected_snp)
		temp_snps[selected_contig].delete(selected_snp) # So it cannot be reselected
		if @select_snps[selected_contig].size == $options.maxsnps or temp_snps[selected_contig].size == 0
			temp_snps.delete_if {|key, value | key == selected_contig}
		end
	end
	# Write VCF & Probes
	search_keys = [] # Array to hold contig/SNP unique search values
	for i in 0..@select_snps.size - 1
		for snp in @select_snps[@select_snps.keys[i]]
			search = @select_snps.keys[i] + "\t" + snp.to_s + "\t"
			search_keys.push(search)
		end
	end
	vcfout = ""
	File.open($options.vcfin, 'r') do |snpsel|
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
	File.open($options.vcfin + "-selected.vcf", 'w') do |write|
		write.puts vcfout
	end
	# Output probe sequences if this option was selected
	if $options.probes
		probesout = ""
		outfilter = ""
		File.open($options.refseq, 'r') do |seqget|
		while line = seqget.gets
			if line[0].chr == ">"
				count = 1 #start at 1 to avoid the > symbol
				chromo = ""
				while line[count].chr != "\n"
					chromo += line[count].chr
					count += 1
				end
				if @select_snps.keys.include?(chromo)
					getsnps = true
				else
					getsnps = false
				end
			elsif getsnps
				for snp in @select_snps[chromo]
					be4 = snp - $options.lenbef
					after = snp + $options.lenaft
					be4 = 0 if be4 < 0
					after = line.length - 2 if after >= line.length - 1 # Final char is line break
					seq = ">" + chromo + "\t" + snp.to_s + "\n" + line[be4 .. after] + "\n"
					probesout += seq
						if $options.filter
							outfilter += seq if filter_probes(line[be4 .. after])
						end
					end
				end
			end
		end
		File.open($options.vcfin + "-selected-probes.fa", 'w') do |write|
			write.puts probesout
		end
		if $options.filter
			File.open($options.vcfin + "-filtered-probes.fa", 'w') do |write|
				write.puts outfilter
			end
		end
	end
end