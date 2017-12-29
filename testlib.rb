#!/usr/bin/env ruby

#-----------------------------------------------------------------------------------------------
# testlib
# Michael G. Campana, 2017
# Smithsonian Conservation Biology Institute
#-----------------------------------------------------------------------------------------------

require_relative 'baitslib'

#====DUMMY OPTIONS BLOCK========================================================================
class Option
	attr_accessor :rna, :no_Ns
	def initialize(rna = false, noNs = false)
		@rna = rna
		@no_Ns = noNs
	end
end
#====TESTED METHODS=============================================================================

#====TEST CONTROL CODE==========================================================================
$options = Option.new
@exit = ""
while @exit == ""
	print "Exit? Return for no. Any string for yes.\n"
	@exit = gets.chomp
	break if @exit != ""
end
