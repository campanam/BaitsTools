#!/usr/bin/ruby
#-----------------------------------------------------------------------------------------------
# select_snps 0.4
# Michael G. Campana, 2015
# Smithsonian Conservation Biology Institute
#-----------------------------------------------------------------------------------------------

# To be implemented
# Output non-reference alleles
# Sequence complexity filter
# Include/Exclude multi-allelic sites
# Minimum PIC/Calculate PIC
# Self-complementarity filter

require 'optparse'
require 'ostruct'
#-----------------------------------------------------------------------------------------------
class Parser
	def self.parse(options)
		# Set defaults
		args = OpenStruct.new
		args.vcfin = ""
		args.every = false
		args.totalsnps = 30000
		args.scale = false
		args.scalehash = {}
		args.maxsnps = 2
		args.distance = 10000
		args.probes = false
		args.refseq = ""
		args.lenbef = 60
		args.lenaft = 59
		args.tiling = false
		args.tileoffset = 60
		args.tiledepth = 2
		args.filter = false
		args.completeprobe = false
		args.probelength = 120
		args.gc = false
		args.mingc = 30.0
		args.maxgc = 50.0
		args.melt = false
		args.na = 0.9
		args.mint = 0.0
		args.maxt = 120.0
		args.params = false
		args.coords = false
		opt_parser = OptionParser.new do |opts|
			opts.banner = "Usage: ruby select_snps-0.4.rb [options]"
			opts.separator ""
			opts.separator "Specific options:"
			opts.separator ""
			opts.on("-i","--input [FILE]", String, "Input VCF file") do |vcf|
				args.vcfin = vcf
			end
			opts.on("-e","--every", "Output probe sequences for every SNP in the input VCF file") do
				args.every = true
				args.probes = true
			end
			opts.on("-t", "--totalsnps [VALUE]", Integer, "Total requested SNPs (Default = 30,000)") do |tsnps|
				args.totalsnps = tsnps
			end
			opts.on("-j", "--scale", "Scale maximum number of SNPs per contig by contig length") do
				args.scale = true
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
			opts.on("-r","--refseq [FILE]", String, "If probe output selected, the reference sequence file") do |ref|
				args.refseq= ref
			end
			opts.on("-b", "--lenbef [VALUE]", Integer, "If probe output selected, the number of bases before the SNP (Default = 60)") do |b4|
				args.lenbef = b4
			end
			opts.on("-a", "--lenaft [VALUE]", Integer, "If probe output selected, the number of bases after the SNP (Default = 59)") do |af|
				args.lenaft = af
			end
			opts.on("-u", "--tiling", "If probe output selected, tile probe sequences.") do
				args.tiling = true
			end
			opts.on("-v", "--offset [VALUE]", Integer, "If tiling selected, base pair offset between tiled probes.") do |toff|
				args.tileoffset = toff
			end
			opts.on("-k", "--depth [VALUE]", Integer, "If tiling selected, requested probes per SNP.") do |tdep|
				args.tiledepth = tdep
			end
			opts.on("-f","--filter", "If probes selected, filter probe sequences") do
				args.filter = true
			end
			opts.on("-w", "--params", "If filtering, output probe statistics table") do
				args.params = true
			end
			opts.on("-c","--complete", "If filtering, require probes be full length") do
				args.completeprobe = true
			end
			opts.on("-g", "--gc", "If filtering, filter by GC content") do
				args.gc = true
			end				
			opts.on("-n","--mingc [VALUE]", Float, "If filtering, minimum GC content (Default = 30.0)") do |mgc|
				args.mingc = mgc
			end
			opts.on("-x","--maxgc [VALUE]", Float, "If filtering, maximum GC content (Default = 50.0)") do |xgc|
				args.maxgc = xgc
			end
			opts.on("-l", "--melt", "If filtering, filter by melting temperature") do
				args.melt = true
			end				
			opts.on("-q","--mint [VALUE]", Float, "If filtering, minimum melting temperature (Default = 0.0)") do |mt|
				args.mint = mt
			end
			opts.on("-z","--maxt [VALUE]", Float, "If filtering, maximum melting temperature (Default = 120.0)") do |xt|
				args.maxt = xt
			end
			opts.on("-s","--na [VALUE]", Float, "If filtering by melting temperature, sodium concentration (Default = 0.9)") do |na|
				args.na = na
			end
			opts.on("-o", "--coords", "Output coordinate tables for the candidate baits.") do
				args.coords = true
			end
			opts.separator ""
			opts.separator "Common options:"
			opts.on_tail("-h","--help", "Show help") do
				print "Welcome to select_snps.\n\n"
				print "To use the interactive interface, enter <ruby select_snps-0.4.rb>.\n\n"
				puts opts
				exit
			end
			opts.on_tail("-v","--version","Show version") do
				print "select_snps 0.4\n"
				exit
			end
		end
		opt_parser.parse!(options)
		return args
	end
end
#-----------------------------------------------------------------------------------------------
def filter_probes(probe)
	keep = true
	keep = false if probe.length < $options.probelength
	gc = 0.0
	for i in 0 ... probe.length
		if probe[i].chr.upcase == "G" or probe[i].chr.upcase == "C"
			gc += 1.0
		end
	end
	gccont = gc/probe.length.to_f
	melt = 79.8 + (58.4 * gccont) + (11.8 * (gccont**2.0)) - 820.0/probe.length.to_f + 18.5 * Math::log($options.na)
	if $options.gc
		keep = false if gccont * 100.0 < $options.mingc
		keep = false if gccont * 100.0 > $options.maxgc
	end
	if $options.melt
		keep = false if melt < $options.mint
		keep = false if melt > $options.maxt
	end
	return [keep, probe.length.to_s + "\t" + gc.to_s + "\t" + melt.to_s + "\t" + keep.to_s + "\n"]
end
#-----------------------------------------------------------------------------------------------
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
		print "Output probes for all SNPs in VCF? (y/n)\n"
		t = gets.chomp.upcase
		if t == "Y" or t == "YES"
			$options.every = true
			$options.probes = true
		end
	end
	if interact and !$options.every
		print "Enter total number of requested SNPs.\n"
		$options.totalsnps = gets.chomp.to_i
	end
	while $options.totalsnps < 1
		print "The total number of SNPs must be greater than 0. Please re-enter.\n"
		$options.totalsnps = gets.chomp.to_i
	end
	if interact and !$options.every
		print "Scale maximum number of SNPs per contig by contig length? (y/n)\n"
		t = gets.chomp.upcase
		if t == "Y" or t == "YES"
			$options.scale = true
		end
	end
	if interact and !$options.every and !$options.scale
		print "Enter maximum number of SNPs per contig.\n"
		$options.maxsnps = gets.chomp.to_i
	end
	while $options.maxsnps < 1 and !$options.scale
		print "The maximum number of SNPs per contig must be greater than 0. Please re-enter.\n"
		$options.maxsnps = gets.chomp.to_i
	end
	if interact and !$options.every
		print "Enter minimum distance between SNPs within a contig.\n"
		$options.distance = gets.chomp.to_i
	end
	while $options.distance < 1
		print "The minimum distance between SNPs must be greater than 0. Please re-enter.\n"
		$options.distance = gets.chomp.to_i
	end
	if interact and !$options.every
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
			print "Tile probes? (y/n)\n"
			t = gets.chomp.upcase
			if t == "Y" or t == "YES"
				$options.tiling = true
				print "Enter probe length.\n"
				$options.probelength = gets.chomp.to_i
				while $options.probelength < 1
					print "Probes must be at least 1 bp long. Re-enter.\n"
					$options.probelength = gets.chomp.to_i
				end
				$options.lenbef = $options.probelength - 1
				$options.lenaft = 0
				print "Enter tiling bp offset.\n"
				$options.tileoffset = gets.chomp.to_i
				print "Enter number of probes per SNP.\n"
				$options.tiledepth = gets.chomp.to_i
			end
		end
		if interact and !$options.tiling
			print "Enter number of bases before SNP in probe.\n"
			$options.lenbef = gets.chomp.to_i
		end
		while $options.lenbef < 0
			print "The number of probe bases before the SNP must be at least 0. Please re-enter.\n"
			$options.lenbef = gets.chomp.to_i
		end
		if interact and !$options.tiling
			print "Enter number of bases after SNP in probe.\n"
			$options.lenaft = gets.chomp.to_i
		end
		while $options.lenaft < 0
			print "The number of probe bases after the SNP must be at least 0. Please re-enter.\n"
			$options.lenaft = gets.chomp.to_i
		end
		$options.probelength = $options.lenbef + $options.lenaft + 1
		if $options.tiling
			$options.lenbef = $options.probelength - 1
			$options.lenaft = 0
			while $options.tileoffset > $options.probelength or $options.tileoffset < 1
				print "Tiling offset cannot be less than 1 or greater than probe length. Please re-enter.\n"
				$options.tileoffset = gets.chomp.to_i
			end
			while $options.tiledepth > $options.probelength/$options.tileoffset or $options.tiledepth < 1
				print "Tiling depth cannot be less than 1 or greater than probe length/tiling offset ratio. Re-enter.\n"
				$options.tiledepth = gets.chomp.to_i
			end
		end
		if interact
			print "Filter probes? (y/n)\n"
			t = gets.chomp.upcase
			if t == "Y" or t == "YES"
				$options.filter = true
			end
			print "Output coordinates table(s) for candidate probes? (y/n)\n"
			t = gets.chomp.upcase
			if t == "Y" or t == "YES"
				$options.coords = true
			end
		end
		if $options.filter
			if interact
				print "Output probe statistics table? (y/n)\n"
				t = gets.chomp.upcase
				if t == "Y" or t == "YES"
					$options.params = true
				end
				print "Require complete length probe? (y/n)\n"
				t = gets.chomp.upcase
				if t == "Y" or t == "YES"
					$options.completeprobe = true
				end
				print "Filter by GC content? (y/n)\n"
				t = gets.chomp.upcase
				if t == "Y" or t == "YES"
					$options.gc = true
				end
				print "Filter by melting temperature? (y/n)\n"
				t = gets.chomp.upcase
				if t == "Y" or t == "YES"
					$options.melt = true
				end
			end
			if $options.gc
				if interact
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
			if $options.melt
				if interact
					print "Enter minimum melting temperature.\n"
					$options.mint = gets.chomp.to_f
					print "Enter maximum melting temperature.\n"
					$options.maxt = gets.chomp.to_f
				end
				while $options.maxt < $options.mint
					print "Maximum melting temperature cannot be less than minimum melting temperature. Re-enter.\n"
					$options.maxt = gets.chomp.to_f
				end
				if interact
					print "Enter sodium concentration.\n"
					$options.na = gets.chomp.to_f
				end
				while $options.na < 0.0
					print "Sodium concentrations cannot be negative. Re-enter.\n"
					$options.na = gets.chomp.to_f
				end
			end
		end
	end
	# Read VCF file
	@snps = {}
	chromosome = ""
	temp_snps = []
	File.open($options.vcfin, 'r') do |snpreg|
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
	@select_snps = {}
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
			@select_snps[selected_contig] = [] if @select_snps[selected_contig].nil?
			@select_snps[selected_contig].push(selected_snp)
			temp_snps[selected_contig].delete(selected_snp) # So it cannot be reselected
			if $options.scale
				maxsize = $options.scalehash[selected_contig]
			else
				maxsize = $options.maxsnps
			end
			if @select_snps[selected_contig].size == maxsize or temp_snps[selected_contig].size == 0
				temp_snps.delete_if {|key, value | key == selected_contig}
			end
		end
	else
		@select_snps = @snps.dup
	end
	# Write VCF & Probes
	search_keys = [] # Array to hold contig/SNP unique search values
	for i in 0..@select_snps.size - 1
		for snp in @select_snps[@select_snps.keys[i]]
			search = @select_snps.keys[i] + "\t" + snp.to_s + "\t"
			search_keys.push(search)
		end
	end
	filter_search_keys = search_keys.dup # Array for filtering
	vcfout = ""
	if !$options.every
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
	end
	# Output probe sequences if this option was selected
	if $options.probes
		probesout = ""
		outfilter = ""
		paramline = "Chromosome:Coordinates\tSNP\tProbeLength\t%GC\tTm\tKept\n"
		coordline = ""
		filtercoordline = ""
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
						if $options.tiling
							tiling = $options.tiledepth
						else
							tiling = 1
							$options.tileoffset = 1
						end
						for tile in 0...tiling
							be4 = snp - $options.lenbef + tile * $options.tileoffset
							after = snp + $options.lenaft + tile * $options.tileoffset
							be4 = 0 if be4 < 0
							after = line.length - 2 if after >= line.length - 1 # Final char is line break
							seq = ">" + chromo + "\t" + snp.to_s + "\n" + line[be4 .. after] + "\n"
							probesout += seq
							coord = chromo + ":" + be4.to_s + "-" + after.to_s + "\n"
							coordline += coord
							if $options.filter
								parameters = filter_probes(line[be4 .. after])
								if parameters[0]
									outfilter += seq
									filtercoordline += coord 
								else
									search = chromo + "\t" + snp.to_s + "\t"
									filter_search_keys.delete(search)
								end
								if $options.params
									paramline += chromo + ":" + be4.to_s + "-" + after.to_s + "\t" + snp.to_s + "\t" + parameters[1]
								end
							end
						end
					end
				end
			end
		end
		File.open($options.vcfin + "-selected-probes.fa", 'w') do |write|
			write.puts probesout
		end
		if $options.coords
			File.open($options.vcfin + "-selected-probes-coords.txt", 'w') do |write|
				write.puts coordline
			end
		end
		if $options.filter
			File.open($options.vcfin + "-filtered-probes.fa", 'w') do |write|
				write.puts outfilter
			end
			vcffilt = ""
			File.open($options.vcfin, 'r') do |snpsel|
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
			File.open($options.vcfin + "-filtered.vcf", 'w') do |write|
				write.puts vcffilt
			end
			if $options.params
				File.open($options.vcfin + "-filtered-params.txt", 'w') do |write|
					write.puts paramline
				end
			end
			if $options.coords
				File.open($options.vcfin + "-filtered-probes-coords.txt", 'w') do |write|
					write.puts filtercoordline
				end
			end
		end
	end
end
