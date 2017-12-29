#!/usr/bin/ruby
#-----------------------------------------------------------------------------------------------
# baitslib 0.2
# Michael G. Campana, 2016
# Smithsonian Conservation Biology Institute
#-----------------------------------------------------------------------------------------------

# To be implemented
# Sequence complexity filter
# Minimum PIC/Calculate PIC
# Self-complementarity filter

#-----------------------------------------------------------------------------------------------
class Fa_Seq #container for fasta sequences
	attr_accessor :header, :circular, :seq
	def initialize(header = "", circular = false, seq = "")
		@header = header
		@circular = circular
		@seq = seq
	end 
end
#-----------------------------------------------------------------------------------------------
def filter_baits(bait)
	keep = true
	keep = false if bait.length < $options.baitlength
	gc = 0.0
	for i in 0 ... bait.length
		if bait[i].chr.upcase == "G" or bait[i].chr.upcase == "C"
			gc += 1.0
		end
	end
	gccont = gc/bait.length.to_f
	melt = 79.8 + (58.4 * gccont) + (11.8 * (gccont**2.0)) - 820.0/bait.length.to_f + 18.5 * Math::log($options.na)
	if $options.gc
		keep = false if gccont * 100.0 < $options.mingc
		keep = false if gccont * 100.0 > $options.maxgc
	end
	if $options.melt
		keep = false if melt < $options.mint
		keep = false if melt > $options.maxt
	end
	return [keep, bait.length.to_s + "\t" + gc.to_s + "\t" + melt.to_s + "\t" + keep.to_s + "\n"]
end
#-----------------------------------------------------------------------------------------------
def write_baits(baitsout = "", outfilter = "", paramline = "", coordline = "", filtercoordline = "")
	File.open($options.infile + "-baits.fa", 'w') do |write|
		write.puts baitsout
	end
	if $options.coords
		File.open($options.infile + "-baits-coords.txt", 'w') do |write|
			write.puts coordline
		end
	end
	if $options.filter
		File.open($options.infile + "-filtered-baits.fa", 'w') do |write|
			write.puts outfilter
		end
		if $options.params
			File.open($options.infile + "-filtered-params.txt", 'w') do |write|
				write.puts paramline
			end
		end
		if $options.coords
			File.open($options.infile + "-filtered-baits-coords.txt", 'w') do |write|
				write.puts filtercoordline
			end
		end
	end
end
#-----------------------------------------------------------------------------------------------
def read_fasta(file)
	seq_array = []
	faseq = nil # Dummy value
	File.open(file, 'r') do |seq|
		while line = seq.gets
			if line[0].chr == ">"
				seq_array.push(faseq) if !faseq.nil? # push previously completed sequence into array
				header = line[1...-1] # Remove final line break and beginning >
				if header[-5..-1] == "#circ" #Adding this to the end of a sequence denotes as circular
					circular = true
					header = header[0...-5] # Remove circ hashtag
				else
					circular = false
				end
				faseq = Fa_Seq.new(header, circular)
			else
				faseq.seq += line[0...-1] #Exclude final line break, allow multi-line fasta
			end
		end
	end
	seq_array.push(faseq) # Put last sequence into fasta array
	return seq_array
end
	