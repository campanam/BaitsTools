#!/usr/bin/ruby
#-----------------------------------------------------------------------------------------------
# filter_probes 0.1
# Michael G. Campana, 2016
# Smithsonian Conservation Biology Institute
#-----------------------------------------------------------------------------------------------

# To be implemented
# Sequence complexity filter
# Minimum PIC/Calculate PIC
# Self-complementarity filter

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
