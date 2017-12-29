#!/usr/bin/ruby
#-----------------------------------------------------------------------------------------------
# select_snps 0.1
# Michael G. Campana, 2015
# Smithsonian Conservation Biology Institute
#-----------------------------------------------------------------------------------------------

begin
	# Data options block
	if ARGV[0] == "--help" or ARGV[0] == "-h"
		print "Welcome to select_snps.\n\n"
		print "To use the interactive interface, enter no options after the command, e.g. <ruby select_snps.rb>\n"
		print "Otherwise, the parameters should be entered in the following order:\n\n"
		print "* Input vcf file\n"
		print "* Total requested SNPs\n"
		print "* Maximum number of SNPs per contig\n"
		print "* Minimum distance between SNPs within a contig\n"
        # GC content range
		# Sequence complexity
		# Include/Exclude multi-allelic sites
		# Minimum PIC/Calculate PIC
        print "* -y to output probes (otherwise leave blank)\n"
        print "* If probe ouput selected, the reference sequence file (otherwise blank)\n"
        print "* If probe output selected, the number of bases before the SNP (otherwise blank)\n"
        print "* If probe output selected, the number of bases after the SNP (otherwise blank)\n\n"
		print "For example, the command <ruby select_snps.rb input.vcf 5000 2 1000> would select 5000 SNPs, with\n"
		print "a maximum of 2 SNPs per contig and a minimum distance of 1000 bp between SNPs from file input.vcf.\n"
        print "Appending, <-y refseq.fa 60 59> to the previous command would ouput probes based on refseq.fa, with\n"
        print "the SNP appearing at the 61st postion of a 120 bp probe.\n\n"
		exit # Kill the program rather than proceed without data
	elsif ARGV[0].nil?
		print "Enter input file.\n"
		@vcfin = gets.chomp
		while !FileTest.exist?(@vcfin)
			print "VCF file not found. Please re-enter.\n"
			@vcfin = gets.chomp
		end
		print "Total requested SNPs.\n"
		@total_snps = gets.chomp.to_i
		while @total_snps < 1
			print "The total number of SNPs must be greater than 0. Please re-enter.\n"
			@total_snps = gets.chomp.to_i
		end
		print "Maximum number of SNPs per contig\n"
		@max_snps = gets.chomp.to_i
		while @max_snps < 1
			print "The maximum number of SNPs per contig must be greater than 0. Please re-enter.\n"
			@max_snps = gets.chomp.to_i
		end
		print "Minimum distance between SNPs within a contig\n"
		@min_distance = gets.chomp.to_i
		while @min_distance < 1
			print "The minimum distance between SNPs must be greater than 0. Please re-enter.\n"
			@min_distance = gets.chomp.to_i
		end
		print "Output probes? (y/n)\n"
        @probes = false
        @probes = true if gets.chomp == "y"
        if @probes
            print "Enter reference sequence.\n"
            @refseq = gets.chomp
            while !FileTest.exist?(@refseq)
                print "Reference sequence not found. Please re-enter.\n"
                @refseq = gets.chomp
            end
            print "Number of bases before SNP in probe\n"
            @lenbef = gets.chomp.to_i
            while @lenbef < 0
                print "The number of probe bases before the SNP must be at least 0. Please re-enter.\n"
                @lenbef = gets.chomp.to_i
            end
            print "Number of bases after SNP in probe\n"
            @lenaft = gets.chomp.to_i
            while @lenaft < 0
                print "The number of probe bases after the SNP must be at least 0. Please re-enter.\n"
                @lenaft = gets.chomp.to_i
            end
        end
	else
		@vcfin = ARGV[0]
		while !FileTest.exist?(@vcfin)
			print "VCF file not found. Please re-enter.\n"
			@vcfin = $stdin.gets.chomp
		end
		@total_snps = ARGV[1].to_i
		while @total_snps < 1
			print "The total number of SNPs must be greater than 0. Please re-enter.\n"
			@total_snps = $stdin.gets.chomp.to_i
		end
		@max_snps = ARGV[2].to_i
		while @max_snps < 1
			print "The maximum number of SNPs per contig must be greater than 0. Please re-enter.\n"
			@max_snps = $stdin.gets.chomp.to_i
		end
		@min_distance = ARGV[3].to_i
		while @min_distance < 1
			print "The minimum distance between SNPs must be greater than 0. Please re-enter.\n"
			@min_distance = $stdin.gets.chomp.to_i
		end
        @probes = false
        @probes = true if ARGV[4] == "-y"
        if @probes
            @refseq = ARGV[5]
            while !FileTest.exist?(@refseq)
                print "Reference sequence not found. Please re-enter.\n"
                @refseq = $stdin.gets.chomp
            end
            @lenbef = ARGV[6].to_i
            while @lenbef < 0
                print "The number of probe bases before the SNP must be at least 0. Please re-enter.\n"
                @lenbef = $stdin.gets.chomp.to_i
            end
            @lenaft = ARGV[7].to_i
            while @lenaft < 0
                print "The number of probe bases after the SNP must be at least 0. Please re-enter.\n"
                @lenaft = $stdin.gets.chomp.to_i
            end
        end
	end
	# Read VCF file
	@snps = {}
	chromosome = ""
	temp_snps = []
	File.open(@vcfin, 'r') do |snpreg|
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
	for i in 1..@total_snps
		break if temp_snps.size == 0 # Stop selecting SNPs when all contigs deleted from consideration
		selected_contig = temp_snps.keys[rand(temp_snps.size)] # Get name of contig
		selected_snp = temp_snps[selected_contig][rand(temp_snps[selected_contig].size)]
		#Check SNP/Probe for suitability
		# Delete SNPs that are too close
        tmp = temp_snps[selected_contig].dup # Duplicate this subsection so that during deletion there are not errors
		for i in 0 .. tmp.size - 1
			if tmp[i] < selected_snp && selected_snp - tmp[i] < @min_distance
				temp_snps[selected_contig].delete(tmp[i])
			elsif tmp[i] > selected_snp && tmp[i] - selected_snp < @min_distance
				temp_snps[selected_contig].delete(tmp[i])
			end
		end
		# Add SNP to selected pool and delete contigs if maximum number of SNPs reached or no remaining SNPs on contig
        @select_snps[selected_contig] = [] if @select_snps[selected_contig].nil?
		@select_snps[selected_contig].push(selected_snp)
        temp_snps[selected_contig].delete(selected_snp) # So it cannot be reselected
		if @select_snps[selected_contig].size == @max_snps or temp_snps[selected_contig].size == 0
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
    File.open(@vcfin, 'r') do |snpsel|
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
	File.open(@vcfin + "-selected.vcf", 'w') do |write|
		write.puts vcfout
	end
    # Output probe sequences if this option was selected
    if @probes
        probesout = ""
        File.open(@refseq, 'r') do |seqget|
            while line = seqget.gets
                if line[0].chr == ">"
                    count = 1 #start at 1 to avoid the > symbol
                    chromo = ""
                    while line[count].chr != " "
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
                        be4 = snp - @lenbef
                        after = snp + @lenaft
                        be4 = 0 if be4 < 0
                        after = line.length - 2 if after >= line.length - 1 # Final char is line break
                        seq = ">" + chromo + "\t" + snp.to_s + "\n" + line[be4 .. after] + "\n"
                        probesout += seq
                    end
                end
            end
        end
        File.open(@vcfin + "-selected-probes.fa", 'w') do |write|
            write.puts probesout
        end
    end
end
