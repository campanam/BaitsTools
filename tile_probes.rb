#!/usr/bin/ruby
#-----------------------------------------------------------------------------------------------
# tile_probes 0.1
# Michael G. Campana, 2016
# Smithsonian Conservation Biology Institute
#-----------------------------------------------------------------------------------------------

require 'optparse'
require 'ostruct'
require './filter_probes.rb'

#-----------------------------------------------------------------------------------------------
class Parser
	def self.parse(options)
		# Set defaults
		args = OpenStruct.new
		args.fa = ""
		args.probelength = 120
		args.tileoffset = 25
		args.filter = false
		args.completeprobe = false
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
			opts.banner = "Usage: ruby tile_probes.rb [options]"
			opts.separator ""
			opts.separator "Specific options:"
			opts.separator ""
			opts.on("-i","--input [FILE]", String, "Input FASTA file") do |fa|
				args.fa = fa
			end
			opts.on("-p", "--length [VALUE]", Integer, "Requested probe length.") do |prblength|
				args.probelength = prblength
			end
			opts.on("-v", "--offset [VALUE]", Integer, "Base pair offset between tiled probes.") do |toff|
				args.tileoffset = toff
			end
			opts.on("-f","--filter", "Filter probe sequences") do
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
				print "Welcome to tile_probes.\n\n"
				print "To use the interactive interface, enter <ruby tile_probes.rb>.\n\n"
				puts opts
				exit
			end
			opts.on_tail("-v","--version","Show version") do
				print "tile_probes 0.1\n"
				exit
			end
		end
		opt_parser.parse!(options)
		return args
	end
end
#-----------------------------------------------------------------------------------------------
begin
	# Data options block
	interact = true if ARGV.size == 0
	$options = Parser.parse(ARGV)
	if interact
		print "Enter input FASTA.\n"
		$options.fa = gets.chomp
	end
	while !FileTest.exist?($options.fa)
		print "FASTA file not found. Please re-enter.\n"
		$options.fa = gets.chomp
	end
	if interact
		print "Enter probe length.\n"
		$options.probelength = gets.chomp.to_i
		while $options.probelength < 1
			print "Probes must be at least 1 bp long. Re-enter.\n"
			$options.probelength = gets.chomp.to_i
		end
		print "Enter tiling bp offset.\n"
		$options.tileoffset = gets.chomp.to_i
		while $options.tileoffset > $options.probelength or $options.tileoffset < 1
			print "Tiling offset cannot be less than 1 or greater than probe length. Please re-enter.\n"
			$options.tileoffset = gets.chomp.to_i
		end
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
	# Read FASTA file
	probesout = ""
	outfilter = ""
	paramline = "Chromosome:Coordinates\tProbeLength\t%GC\tTm\tKept\n"
	coordline = ""
	filtercoordline = ""
	chrtp = "" #dummy variable for chromosome name
	circular = false #variable to determine whether a sequence is circular
	File.open($options.fa, 'r') do |seq|
		while line = seq.gets
			if line[0].chr == ">"
				chrtp = line[0...-1]
				if chrtp[-5..-1] == "#circ" #Adding this to the end of a sequence denotes as circular
					circular = true
				else
					circular = false
				end
			else
				sequence = line
				seqst = 1 # Beginning of probe coordinate
				while seqst < line.length # Stop the loop once end of sequence reached -- keeping in mind line break
					seqend = seqst+$options.probelength-1 #End of probe coordinate
					if seqend > (line.length-1) and !circular #correct for running off end of linear sequence
						seqend = line.length-1 # Last character is line break
						prb = line[seqst-1..seqend-1] #Correct for 0-based counting
						rng = seqst.to_s + "-" + seqend.to_s
					elsif seqend > line.length-1 and circular #add beginning of sequence to end
						seqend -= (line.length-1)
						prb = line[seqst-1..line.length-2] + line[0..seqend-1] #Correct for 0-based counting and final line break
						rng = seqst.to_s + "-" + seqend.to_s
					else
						prb = line[seqst-1..seqend-1] #Correct for 0-based counting
						rng = seqst.to_s + "-" + seqend.to_s
					end
					probesout += chrtp + "_" + rng + "\n" + prb + "\n"
					coordline += chrtp + ":" + rng + "\n"
					if $options.filter
						flt = filter_probes(prb)
						if flt[0]
							outfilter += chrtp + "_" + rng + "\n" + prb + "\n"
							filtercoordline += chrtp + ":" + rng + "\n"
						end
						if $options.params
							paramline += chrtp + ":" + rng + "\t" + "\t" + flt[1]
						end
					end
					seqst += $options.tileoffset
				end
					
			end
		end
	end
	File.open($options.fa + "-probes.fa", 'w') do |write|
		write.puts probesout
	end
	if $options.coords
		File.open($options.fa + "-probes-coords.txt", 'w') do |write|
			write.puts coordline
		end
	end
	if $options.filter
		File.open($options.fa + "-filtered-probes.fa", 'w') do |write|
			write.puts outfilter
		end
		if $options.params
			File.open($options.fa + "-filtered-params.txt", 'w') do |write|
				write.puts paramline
			end
		end
		if $options.coords
			File.open($options.fa + "-filtered-probes-coords.txt", 'w') do |write|
				write.puts filtercoordline
			end
		end
	end
end
