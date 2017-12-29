#!/usr/bin/ruby
#-----------------------------------------------------------------------------------------------
# checkbaits
CHECKBAITSVER = "0.3"
# Copyright (C) 2017 Michael G. Campana
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#-----------------------------------------------------------------------------------------------

def checkbaits
	# Process FASTA file
	print "** Reading FASTA/FASTQ **\n"
	outfilter = ""
	paramline = "Bait\tBaitLength\t%GC\tTm\tMeanQuality\tMinQuality\tKept\n"
	seq_array = read_fasta($options.infile)
	print "** Filtering baits **\n"
	for seq in seq_array
		if seq.fasta
			flt = filter_baits(seq.seq, [$options.fasta_score])
		else
			flt = filter_baits(seq.seq, seq.numeric_quality)
		end
		if flt[0]
			seq.seq.gsub!("T","U") if $options.rna # RNA output handling
			outfilter += ">" + seq.header + "\n" + seq.seq + "\n"
		end
		if $options.params
			paramline += seq.header + "\t" + flt[1]
		end
	end
	write_baits("", outfilter, paramline)
end
