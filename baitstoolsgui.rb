#!/usr/bin/env ruby
#-----------------------------------------------------------------------------------------------
# baitstoolsgui
BAITSTOOLSGUI = "1.7.0"
# Michael G. Campana, 2017-2021
# Smithsonian Conservation Biology Institute
#-----------------------------------------------------------------------------------------------

require 'tk'
require 'tkextlib/tile'
require 'ostruct'
require 'shell'
require_relative 'baitslib.rb'

#-----------------------------------------------------------------------------------------------
def start_baitstools
	# Generate basic command line
	cmdline = "baitstools.rb " + $options.algorithm + " -i " + resolve_unix_path($options.infile.to_s)
	case $options.algorithm
	when "vcf2baits", "stacks2baits"
		if $options.every == 1
			cmdline << " -e"
			cmdline << " -L" + $options.baitlength + " -O" + $options.tileoffset + " -b" + $options.lenbef + " -k" + $options.tiledepth
		else
			cmdline << " -t" + $options.totalsnps + " -d" + $options.distance if $options.taxafile == ""
			if $options.scale == 1
				cmdline << " -j"
			else
				cmdline << " -m" + $options.maxsnps
			end
			if $options.no_baits == 1
				cmdline << " -p"
			else
				cmdline << " -L" + $options.baitlength + " -O" + $options.tileoffset + " -b" + $options.lenbef + " -k" + $options.tiledepth
			end
			if $options.taxafile != ""
				cmdline << " --taxafile " + resolve_unix_path($options.taxafile.to_s)
				cmdline << " --taxacount " + $options.taxacount
				cmdline << " --popcategories " + $options.popcategories if $options.popcategories != ""
			end
			cmdline << " --previousbaits " + resolve_unix_path($options.previousbaits.to_s) if $options.previousbaits != ""
		end
		cmdline << " -r " + resolve_unix_path($options.refseq.to_s) unless $options.no_baits == 1
		cmdline << " -a" if $options.alt_alleles == 1
		cmdline << " -V " + $options.varqual if $options.varqual_filter == 1
		cmdline << " -S" if $options.sort == 1
		cmdline << " -H -A" + $options.alpha.to_s if $options.hwe == 1
	else
		if $options.algorithm == "annot2baits" or $options.algorithm == "bed2baits" or $options.algorithm == "blast2baits"
			cmdline << " --list " + $options.list_format if $options.algorithm == "bed2baits"
			cmdline << " -r " + resolve_unix_path($options.refseq.to_s) + " -P" + $options.pad			
		end
		cmdline << " -L" + $options.baitlength
		cmdline << " -O" + $options.tileoffset unless $options.algorithm == "checkbaits"
		if $options.algorithm == "pyrad2baits"
			cmdline << " -I" + $options.minind
			cmdline << " -W " + $options.strategy
		end
		if $options.algorithm == "aln2baits" or ($options.algorithm == "pyrad2baits" && $options.strategy == "alignment")
			cmdline << " -H " + $options.haplodef
		elsif $options.algorithm == "annot2baits"
			cmdline << " -U " + $options.features.value.upcase
		elsif $options.algorithm == "blast2baits"
			cmdline << " --percid " + $options.percid + " --blastlen " + $options.blastlen
			cmdline << " --evalue " + $options.evalue if $options.evalue_filter == 1
		end
		if $options.algorithm == "pyrad2baits" && $options.strategy != "alignment"
			cmdline << " --uncollapsedref" if $options.uncollapsed_ref == 1
			cmdline << " -t" + $options.totalsnps + " -m" + $options.maxsnps + " -d" + $options.distance + " -k" + $options.tiledepth
			cmdline << " -a" if $options.alt_alleles
		end
	end
	cmdline << " --inbed " + resolve_unix_path($options.inbed.to_s) if $options.inbed != ""
	cmdline << " -o " + resolve_unix_path($options.outprefix.to_s)
	cmdline << " -Z " + resolve_unix_path($options.outdir.to_s)
	cmdline << " -l" if $options.log == 1
	cmdline << " -B" if $options.coords == 1
	cmdline << " -E" if $options.rbed == 1
	cmdline << " --shuffle" if $options.shuffle == 1
	cmdline << " -D" if $options.ncbi == 1
	cmdline << " --phred64" if $options.phred64 == 1
	cmdline << " -C" if $options.collapse_ambiguities == 1
	cmdline << " -Y" if $options.rna == 1
	cmdline << " -G " + $options.gaps
	cmdline << " -5 " + $options.fiveprime if $options.fiveprime != ""
	cmdline << " -3 " + $options.threeprime if $options.threeprime != ""
	cmdline << " --fillin " + $options.fillin if $options.fillin != ""
	cmdline << " -X" + $options.threads + " --rng " + $options.rng
	cmdline << " --gzip" if $options.gzip == 1
	# Generate filtration options
	cmdline << " -w" if $options.params == 1
	cmdline << " --disable-lc" if $options.no_lc == 1
	cmdline << " -c" if $options.completebait == 1
	cmdline << " -N" if $options.no_Ns == 1
	cmdline << " --noaddenda" if $options.noaddenda == 1
	cmdline << " -n" + $options.mingc if $options.mingc_filter == 1
	cmdline << " -x" + $options.maxgc if $options.maxgc_filter == 1
	cmdline << " -q" + $options.mint if $options.mint_filter == 1
	cmdline << " -z" + $options.maxt if $options.maxt_filter == 1
	if $options.mint_filter == 1 or $options.maxt_filter == 1
		cmdline << " -T " + $options.bait_type + " -s" + $options.na + " -f" + $options.formamide
	end
	cmdline << " -K" + $options.maxmask if $options.maxmask_filter == 1
	cmdline << " -J" + $options.maxhomopoly if $options.maxhomopoly_filter == 1
	cmdline << " -y" + $options.lc if $options.lc_filter == 1
	cmdline << " -Q" + $options.meanqual if $options.meanqual_filter == 1
	cmdline << " -M" + $options.minqual if $options.minqual_filter == 1
	cmdline << " -F" + $options.fasta_score if ($options.meanqual_filter == 1 or $options.minqual_filter == 1)
	Tk::messageBox :message => 'Starting BaitsTools with the command: ' + cmdline
	system(cmdline)
	exit
end
#-----------------------------------------------------------------------------------------------
def pyrad_windows
	minind = TkLabel.new($root) do
		text 'Min individuals'
		font TkFont.new('times 20')
		place('x' => 540, 'y' => 150)
		pady 10
	end
	minindentry = TkEntry.new($root) do
		textvariable $options.minind
		borderwidth 5
		font TkFont.new('times 12')
		place('x' => 690, 'y' => 160)
		width 10
	end
	strategy = TkLabel.new($root) do
		text 'Strategy'
		font TkFont.new('times 20')
		place('x' => 50, 'y' => 200)
		pady 10
	end
	strategyselect = Tk::Tile::Combobox.new($root) do
		textvariable $options.strategy
		values ["alignment", "SNPs", "informative"]
		state "readonly"
		width 10
		height 2
		place('x' => 240, 'y' => 210)
	end
	strategyselect.bind("<ComboboxSelected>") do
		update_strategy
	end
	$widgets.push($haplo, $haploselect)
	$totalsnps = TkLabel.new($root) do
		text 'Total SNPs'
		font TkFont.new('times 20')
		place('x' => 50, 'y' => 300)
		pady 10
	end
	$totalsnpentry = TkEntry.new($root) do
		textvariable $options.totalsnps
		borderwidth 5
		font TkFont.new('times 12')
		place('x' => 180, 'y' => 310)
		width 10
	end
	$alts = TkCheckButton.new($root) do
		variable $options.alt_alleles
		text "Alternate alleles"
		place('x' => 550, 'y' => 350)
	end
	$uncollapsedref = TkCheckButton.new($root) do
		variable $options.uncollapsed_ref
		text "Uncollapsed reference"
		place('x' => 550, 'y' => 200)
	end
	$maxsnps = TkLabel.new($root) do
		text 'SNPs per locus'
		font TkFont.new('times 20')
		place('x' => 280, 'y' => 300)
		pady 10
	end
	$maxsnpentry = TkEntry.new($root) do
		textvariable $options.maxsnps
		borderwidth 5
		font TkFont.new('times 12')
		place('x' => 420, 'y' => 310)
		width 10
	end
	$distance = TkLabel.new($root) do
		text 'Variant distance'
		font TkFont.new('times 20')
		place('x' => 540, 'y' => 300)
		pady 10
	end
	$distanceentry = TkEntry.new($root) do
		textvariable $options.distance
		borderwidth 5
		font TkFont.new('times 12')
		place('x' => 690, 'y' => 310)
		width 10
	end
	$lenbef = TkLabel.new($root) do
		text 'Length before'
		font TkFont.new('times 20')
		place('x' => 50, 'y' => 350)
		pady 10
	end
	$lenbefentry = TkEntry.new($root) do
		textvariable $options.lenbef
		borderwidth 5
		font TkFont.new('times 12')
		place('x' => 180, 'y' => 360)
		width 10
	end
	$tiledepth = TkLabel.new($root) do
		text 'Tiles per SNP'
		font TkFont.new('times 20')
		place('x' => 280, 'y' => 350)
		pady 10
	end
	$tiledepthentry = TkEntry.new($root) do
		textvariable $options.tiledepth
		borderwidth 5
		font TkFont.new('times 12')
		place('x' => 421, 'y' => 360)
		width 10
	end
	$widgets.push(minind, minindentry, strategy, strategyselect, $alts, $uncollapsedref, $maxsnps, $maxsnpentry, $distance, $distanceentry, $totalsnps, $totalsnpentry, $lenbef, $lenbefentry, $tiledepth, $tiledepthentry)
	configure_buttons([$alts, $uncollapsedref]) # Do not configure widgets since will configure everything
	update_strategy
	$alts.width = $uncollapsedref.width = 20
end
#-----------------------------------------------------------------------------------------------
def update_strategy
	if $options.strategy == "alignment"
		$maxsnpentry.state = $maxsnps.state = "disabled"
		$distanceentry.state = $distance.state = "disabled"
		$totalsnpentry.state = $totalsnps.state = "disabled"
		$lenbef.state = $lenbefentry.state = "disabled"
		$tiledepth.state = $tiledepthentry.state = "disabled"
		$alts.state = "disabled"
		$uncollapsedref.state = "disabled"
		$haplo.state = $haploselect.state = "normal"
	else
		$maxsnpentry.state = $maxsnps.state = "normal"
		$distanceentry.state = $distance.state = "normal"
		$totalsnpentry.state = $totalsnps.state = "normal"
		$lenbef.state = $lenbefentry.state = "normal"
		$tiledepth.state = $tiledepthentry.state = "normal"
		$alts.state = "normal"
		$uncollapsedref.state = "normal"
		$haplo.state = $haploselect.state = "disabled"
	end
end
#-----------------------------------------------------------------------------------------------
def snp_windows
	$totalsnps = TkLabel.new($root) do
		text 'Total SNPs'
		font TkFont.new('times 20')
		place('x' => 50, 'y' => 200)
		pady 10
	end
	$totalsnpentry = TkEntry.new($root) do
		textvariable $options.totalsnps
		borderwidth 5
		font TkFont.new('times 12')
		place('x' => 180, 'y' => 210)
		width 10
	end
	$nobaits = TkCheckButton.new($root) do
		variable $options.no_baits
		text "No baits"
		place('x' => 50, 'y' => 250)
		command '$reffile.state == "disabled" ? $reffile.state = $reflabel.state = $every.state = $alts.state = "normal" : $reffile.state = $reflabel.state = $every.state = $alts.state = "disabled"'
	end
	$every = TkCheckButton.new ($root) do
		variable $options.every
		text "Every variant"
		place('x' => 300, 'y' => 250)
		command 'update_every'
	end
	$alts = TkCheckButton.new ($root) do
		variable $options.alt_alleles
		text "Alternate alleles"
		place('x' => 550, 'y' => 250)
		command 'update_alt_alleles'
	end
	$scale = TkCheckButton.new($root) do
		variable $options.scale
		text "Scale variants"
		place('x' => 50, 'y' => 300)
		command '$maxsnpentry.state == "disabled" ? $maxsnpentry.state = $maxsnps.state = "normal" : $maxsnpentry.state = $maxsnps.state = "disabled"'
	end
	$maxsnps = TkLabel.new($root) do
		text 'Variants per contig'
		font TkFont.new('times 20')
		place('x' => 250, 'y' => 300)
		pady 10
	end
	$maxsnpentry = TkEntry.new($root) do
		textvariable $options.maxsnps
		borderwidth 5
		font TkFont.new('times 12')
		place('x' => 420, 'y' => 310)
		width 10
	end
	$distance = TkLabel.new($root) do
		text 'Variant distance'
		font TkFont.new('times 20')
		place('x' => 540, 'y' => 300)
		pady 10
	end
	$distanceentry = TkEntry.new($root) do
		textvariable $options.distance
		borderwidth 5
		font TkFont.new('times 12')
		place('x' => 690, 'y' => 310)
		width 10
	end
	$lenbef = TkLabel.new($root) do
		text 'Length before'
		font TkFont.new('times 20')
		place('x' => 540, 'y' => 350)
		pady 10
	end
	$lenbefentry = TkEntry.new($root) do
		textvariable $options.lenbef
		borderwidth 5
		font TkFont.new('times 12')
		place('x' => 690, 'y' => 360)
		width 10
	end
	length_window(350)
	offset_window(350, 290)
	$tiledepth = TkLabel.new($root) do
		text 'Tiles per SNP'
		font TkFont.new('times 20')
		place('x' => 50, 'y' => 400)
		pady 10
	end
	$tiledepthentry = TkEntry.new($root) do
		textvariable $options.tiledepth
		borderwidth 5
		font TkFont.new('times 12')
		place('x' => 180, 'y' => 410)
		width 10
	end
	$widgets.push($nobaits, $every, $alts, $scale, $maxsnps, $maxsnpentry, $distance, $distanceentry, $totalsnps, $totalsnpentry, $lenbef, $lenbefentry, $tiledepth, $tiledepthentry)
	configure_buttons([$nobaits, $every, $alts, $scale]) # Do not configure widgets since will configure everything
	$reffile.state = $reflabel.state = $every.state = $alts.state = "disabled" if $options.no_baits == 1
	$nobaits.state = "disabled" if $options.alt_alleles == 1
	$maxsnpentry.state = $maxsnps.state = "disabled" if $options.scale == 1
	if $options.every == 1
		$maxsnpentry.state = $maxsnps.state = "disabled"
		$distanceentry.state = $distance.state = "disabled"
		$totalsnpentry.state = $totalsnps.state = "disabled"
		$nobaits.state = "disabled"
		$scale.state = "disabled"
	end
	$alts.width = 20
	$scale.width = 15
end
#-----------------------------------------------------------------------------------------------
def update_every
	update_alt_alleles
	if $scale.state == "normal"
		$scale.state = "disabled"
		$maxsnpentry.state = $maxsnps.state = "disabled"
		$distanceentry.state = $distance.state = "disabled"
		$totalsnpentry.state = $totalsnps.state = "disabled"
	else
		$scale.state = "normal"
		$maxsnpentry.state = $maxsnps.state = "normal"
		$distanceentry.state = $distance.state = "normal"
		$totalsnpentry.state = $totalsnps.state = "normal"
	end
end
#-----------------------------------------------------------------------------------------------
def update_alt_alleles
	if $options.every == 1
		$reffile.state = "normal"
		$reflabel.state = "normal"
		$nobaits.state = "disabled"
	elsif $options.every == 0 && $options.alt_alleles == 1
		$options.no_baits.value = 0
		$reffile.state = "normal"
		$reflabel.state = "normal"
		$nobaits.state = "disabled"
	else
		$options.no_baits.value = 0
		$nobaits.state = "normal"
		$reffile.state = "normal"
		$reflabel.state = "normal"
	end
end
#-----------------------------------------------------------------------------------------------
def list_format_window
	$listformat = TkLabel.new($root) do
		text 'File format'
		font TkFont.new('times 20')
		place('x' => 50, 'y' => 150)
		pady 10
	end
	$listselect = Tk::Tile::Combobox.new($root) do
		textvariable $options.list_format
		values ["BED", "GATK", "Picard"]
		state "readonly"
		width 10
		height 2
		place('x' => 160, 'y' => 160)
	end
	$widgets.push($listformat, $listselect)
end
#-----------------------------------------------------------------------------------------------
def haplodef_window
	$haplo = TkLabel.new($root) do
		text 'Haplotype definition'
		font TkFont.new('times 20')
		place('x' => 50, 'y' => 250)
		pady 10
	end
	$haploselect = Tk::Tile::Combobox.new($root) do
		textvariable $options.haplodef
		values ["haplotype", "variant"]
		state "readonly"
		width 10
		height 2
		place('x' => 240, 'y' => 260)
	end
	$widgets.push($haplo, $haploselect)
end
#-----------------------------------------------------------------------------------------------
def reference_window(winy = 150)
	$reffile = TkButton.new($root) do
		text 'Reference sequence'
		command '$options.refseq.value = Tk.getOpenFile'
		place('x' => 50, 'y' => winy)
	end
	configure_buttons([$reffile])
	$reffile.width = 20
	$reflabel = TkLabel.new($root) do
		textvariable $options.refseq
		font TkFont.new('times 12')
		place('x' => 300, 'y' => winy)
		pady 10
	end
	$widgets.push($reffile, $reflabel)
end
#-----------------------------------------------------------------------------------------------
def inbed_window(winy = 150)
	$inbedfile = TkButton.new($root) do
		text 'BED file (optional)'
		command '$options.inbed.value = Tk.getOpenFile'
		place('x' => 50, 'y' => winy)
	end
	configure_buttons([$inbedfile])
	$inbedfile.width = 20
	$inbedlabel = TkLabel.new($root) do
		textvariable $options.inbed
		font TkFont.new('times 12')
		place('x' => 300, 'y' => winy)
		pady 10
	end
	$widgets.push($inbedfile, $inbedlabel)
end
#-----------------------------------------------------------------------------------------------
def feature_window
	feat = TkLabel.new($root) do
		text 'Features (comma-separated list)'
		font TkFont.new('times 20')
		place('x' => 50, 'y' => 350)
		pady 10
	end
	featentry = TkEntry.new($root) do
		textvariable $options.features
		borderwidth 5
		font TkFont.new('times 12')
		place('x' => 340, 'y' => 360)
		width 50
	end
	$widgets.push(feat, featentry)
end
#-----------------------------------------------------------------------------------------------
def blast_window
	percid = TkLabel.new($root) do
		text 'Min. % identity'
		font TkFont.new('times 20')
		place('x' => 50, 'y' => 350)
		pady 10
	end
	percidentry = TkEntry.new($root) do
		textvariable $options.percid
		borderwidth 5
		font TkFont.new('times 12')
		place('x' => 180, 'y' => 360)
		width 10
	end
	blastlen = TkLabel.new($root) do
		text 'Min. hit length'
		font TkFont.new('times 20')
		place('x' => 330, 'y' => 350)
		pady 10
	end
	blastlenentry = TkEntry.new($root) do
		textvariable $options.blastlen
		borderwidth 5
		font TkFont.new('times 12')
		place('x' => 480, 'y' => 360)
		width 10
	end
	evalue = TkCheckButton.new ($root) do
		variable $options.evalue_filter
		text 'Max. E-value'
		font TkFont.new('times 20')
		place('x' => 50, 'y' => 400)
		command '$evalueentry.state == "disabled" ? $evalueentry.state = "normal" : $evalueentry.state = "disabled"'
	end
	$evalueentry = TkEntry.new($root) do
		textvariable $options.evalue
		borderwidth 5
		font TkFont.new('times 12')
		place('x' => 220, 'y' => 410)
		width 10
	end
	$widgets.push(percid, percidentry, blastlen, blastlenentry, $evalueentry)
	configure_buttons([evalue])
	evalue.width = 15
	$widgets.push(evalue)
	$evalueentry.state = "disabled" if $options.evalue_filter == 0
end
#-----------------------------------------------------------------------------------------------
def sort_windows
	sort = TkCheckButton.new($root) do
		variable $options.sort
		text "Sort by pops"
		place('x' => 50, 'y' => 450)
		command 'update_sort'
	end
	$hwe = TkCheckButton.new ($root) do
		variable $options.hwe
		text "Sort by HWE"
		place('x' => 240, 'y' => 450)
		command '$alpha.state == "disabled" ? $alpha.state = $alphaentry.state = "normal" : $alpha.state = $alphaentry.state = "disabled"'
	end
	$alpha = TkLabel.new($root) do
		text 'Alpha'
		font TkFont.new('times 20')
		place('x' => 440, 'y' => 450)
		pady 10
	end
	$alphaentry = TkEntry.new($root) do
		textvariable $options.alpha
		borderwidth 5
		font TkFont.new('times 12')
		place('x' => 500, 'y' => 460)
		width 10
	end
	$widgets.push(sort, $hwe, $alpha, $alphaentry)
	configure_buttons([sort, $hwe])
	$hwe.state = "disabled" if $options.sort == 0
	$alpha.state = $alphaentry.state = "disabled" if $options.hwe == 0
end
#-----------------------------------------------------------------------------------------------
def update_sort
	if $hwe.state == "disabled"
		$hwe.state = "normal"
	else 
		$hwe.state = $alpha.state = $alphaentry.state = "disabled"
		$options.hwe.value = 0
	end
end
#-----------------------------------------------------------------------------------------------
def pad_window(winy = 200)
	$pad = TkLabel.new($root) do
		text 'Pad length'
		font TkFont.new('times 20')
		place('x' => 50, 'y' => winy)
		pady 10
	end
	$padentry = TkEntry.new($root) do
		textvariable $options.pad
		borderwidth 5
		font TkFont.new('times 12')
		place('x' => 180, 'y' => winy + 10)
		width 10
	end
	$widgets.push($pad, $padentry)
end
#-----------------------------------------------------------------------------------------------
def length_window(winy = 150)
	$baitlength = TkLabel.new($root) do
		text 'Bait length'
		font TkFont.new('times 20')
		place('x' => 50, 'y' => winy)
		pady 10
	end
	$baitlengthentry = TkEntry.new($root) do
		textvariable $options.baitlength
		borderwidth 5
		font TkFont.new('times 12')
		place('x' => 180, 'y' => winy + 10)
		width 10
	end
	$widgets.push($baitlength, $baitlengthentry)
end
#-----------------------------------------------------------------------------------------------
def offset_window(winy = 200, winx = 50)
	$offset = TkLabel.new($root) do
		text 'Tiling offset'
		font TkFont.new('times 20')
		place('x' => winx, 'y' => winy)
		pady 10
	end
	$offsetentry = TkEntry.new($root) do
		textvariable $options.tileoffset
		borderwidth 5
		font TkFont.new('times 12')
		place('x' => winx + 130, 'y' => winy + 10)
		width 10
	end
	 $widgets.push($offset, $offsetentry)
end
#-----------------------------------------------------------------------------------------------
def varqual_window
	varqual = TkCheckButton.new($root) do
		variable $options.varqual_filter
		text "Min variant Q"
		place('x' => 300, 'y' => 200)
		command '$varqualentry.state == "disabled" ? $varqualentry.state = "normal" : $varqualentry.state = "disabled"'
	end
	$varqualentry = TkEntry.new($root) do
		textvariable $options.varqual
		borderwidth 5
		font TkFont.new('times 12')
		place('x' => 450, 'y' => 210)
		width 10
	end
	configure_buttons([varqual])
	$widgets.push(varqual, $varqualentry)
	$varqualentry.state = "disabled" if $options.varqual_filter == 0
	varqual.width = 12
end
#-----------------------------------------------------------------------------------------------
def general_window
	clear_widgets
	$labelVar.value = "General Options"
	prefix = TkLabel.new($root) do
		text 'Output prefix'
		font TkFont.new('times 20')
		place('x' => 50, 'y' => 100)
		pady 10
	end
	prefixentry = TkEntry.new($root) do
		textvariable $options.outprefix
		borderwidth 5
		font TkFont.new('times 12')
		place('x' => 200, 'y' => 110)
		width 30
	end
	outdir = TkButton.new($root) do
		text 'Output directory'
		command '$options.outdir.value = Tk.chooseDirectory'
		place('x' => 50, 'y' => 150)
	end
	outdirlabel = TkLabel.new($root) do
		textvariable $options.outdir
		font TkFont.new('times 12')
		place('x' => 300, 'y' => 150)
		pady 10
	end
	bed = TkCheckButton.new($root) do
		variable $options.coords
		text "Output BED"
		place('x' => 50, 'y' => 200)
	end
	rbed = TkCheckButton.new($root) do
		variable $options.rbed
		text "Output relative BED"
		place('x' => 300, 'y' => 200)
	end
	shuffle = TkCheckButton.new($root) do
		variable $options.shuffle
		text "Shuffle baits"
		place('x' => 550, 'y' => 200)
	end
	log = TkCheckButton.new($root) do
		variable $options.log
		text "Output log"
		place('x' => 50, 'y' => 250)
	end
	ncbi = TkCheckButton.new ($root) do
		variable $options.ncbi
		text "NCBI headers"
		place('x' => 300, 'y' => 250)
	end
	rna = TkCheckButton.new ($root) do
		variable $options.rna
		text "Output RNA"
		place('x' => 550, 'y' => 250)
	end
	rc = TkCheckButton.new($root) do
		variable $options.rc
		text "Reverse complement"
		place('x' => 50, 'y' => 300)
	end
	phred64 = TkCheckButton.new($root) do
		variable $options.phred64
		text "Phred64 Input"
		place('x' => 300, 'y' => 300)
	end
	gzip = TkCheckButton.new($root) do
		variable $options.gzip
		text "Gzip Output"
		place('x' => 550, 'y' => 300)
	end
	gaps = TkLabel.new($root) do
		text 'Gap strategy'
		font TkFont.new('times 20')
		place('x' => 50, 'y' => 350)
		pady 10
	end
	gapselect = Tk::Tile::Combobox.new($root) do
		textvariable $options.gaps
		values ["include", "exclude", "extend"]
		state "readonly"
		width 10
		height 2
		place('x' => 170, 'y' => 360)
	end
	threads = TkLabel.new($root) do
		text 'Threads'
		font TkFont.new('times 20')
		place('x' => 310, 'y' => 350)
		pady 10
	end
	threadentry = TkEntry.new($root) do
		textvariable $options.threads
		borderwidth 5
		font TkFont.new('times 12')
		place('x' => 390, 'y' => 360)
		width 10
	end
	fillin = TkLabel.new($root) do
		text 'Fill-in Motif'
		font TkFont.new('times 20')
		place('x' => 520, 'y' => 350)
		pady 10
	end
	fillinentry = TkEntry.new($root) do
		textvariable $options.fillin
		borderwidth 5
		font TkFont.new('times 12')
		place('x' => 640, 'y' => 360)
		width 10
	end
	threeprime = TkLabel.new($root) do
		text '3\' Sequence'
		font TkFont.new('times 20')
		place('x' => 270, 'y' => 400)
		pady 10
	end
	threeprimeentry = TkEntry.new($root) do
		textvariable $options.threeprime
		borderwidth 5
		font TkFont.new('times 12')
		place('x' => 390, 'y' => 410)
		width 10
	end
	fiveprime = TkLabel.new($root) do
		text '5\' Sequence'
		font TkFont.new('times 20')
		place('x' => 50, 'y' => 400)
		pady 10
	end
	fiveprimeentry = TkEntry.new($root) do
		textvariable $options.fiveprime
		borderwidth 5
		font TkFont.new('times 12')
		place('x' => 170, 'y' => 410)
		width 10
	end
	noaddenda = TkCheckButton.new($root) do
		variable $options.noaddenda
		text "Exclude addenda from QC"
		place('x' => 520, 'y' => 400)
	end
	rng = TkLabel.new($root) do
		text 'Random number seed'
		font TkFont.new('times 20')
		place('x' => 50, 'y' => 450)
		pady 10
	end
	rngentry = TkEntry.new($root) do
		textvariable $options.rng
		borderwidth 5
		font TkFont.new('times 12')
		place('x' => 240, 'y' => 460)
		width 80
	end
	configure_buttons([outdir, bed, rbed, shuffle, log, ncbi, rna, rc, phred64, gzip, noaddenda])
	bed.state = shuffle.state = "disabled" if $options.algorithm == "checkbaits"
	ncbi.state = "disabled" if $options.algorithm == "pyrad2baits"
	rbed.state = "disabled" unless $options.algorithm == "annot2baits" or $options.algorithm == "bed2baits" or $options.algorithm == "tilebaits" or $options.algorithm == "aln2baits"
	outdir.width = rbed.width = rc.width = noaddenda.width = 20
	ncbi.width = phred64.width = 15
	$widgets.push(prefix, prefixentry, outdir, outdirlabel, bed, rbed, shuffle, log, ncbi, rna, rc, phred64, gzip, gaps, gapselect, threads, threadentry, rng, rngentry, threeprime, threeprimeentry, fiveprime, fiveprimeentry, noaddenda, fillin, fillinentry)
end
#-----------------------------------------------------------------------------------------------
def subcommand_window(subcommand)
	$root.configure("title", subcommand)
	subcommand == "pyrad2baits" ? $options.maxsnps = TkVariable.new(1) : $options.maxsnps = TkVariable.new(2) # Maximum SNPs per contig
	subcommand == "pyrad2baits" ? $options.distance = TkVariable.new(100) : $options.distance = TkVariable.new(10000) # Minimum distance between SNPs in a contig
	clear_widgets
	$back_btn.state = "normal"
	$options.algorithm = subcommand
	$labelVar.value =  subcommand + " Options"
	case subcommand
	when "aln2baits"
		length_window
		offset_window
		haplodef_window
		inputlabel = "Input FASTA/FASTQ"
	when "annot2baits"
		reference_window
		pad_window
		length_window(250)
		offset_window(300)
		feature_window
		inputlabel = "Input GFF/GTF"
	when "bed2baits"
		list_format_window
		reference_window(200)
		pad_window(250)
		length_window(300)
		offset_window(350)
		inputlabel = "Input BED/Interval List"
	when "blast2baits"
		reference_window
		pad_window
		length_window(250)
		offset_window(300)
		inputlabel = "Input BLAST table"
		blast_window
	when "checkbaits"
		length_window(200)
		inbed_window
		inputlabel = "Input FASTA/FASTQ"
	when "pyrad2baits"
		length_window
		offset_window(150,290)
		haplodef_window
		pyrad_windows
		inputlabel = "Input LOCI"
	when "stacks2baits"
		reference_window
		snp_windows
		sort_windows
		inputlabel = "Input Stacks TSV"
	when "tilebaits"
		length_window
		offset_window
		inputlabel = "Input FASTA/FASTQ"
	when "vcf2baits"
		reference_window
		snp_windows
		varqual_window
		inputlabel = "Input VCF"
	end
	inputfile = TkButton.new($root) do
		command '$options.infile.value = Tk.getOpenFile'
		place('x' => 50, 'y' => 100)
	end
	inputfile.text = inputlabel
	configure_buttons([inputfile])
	inputfile.width = 20
	filelabel = TkLabel.new($root) do
		textvariable $options.infile
		font TkFont.new('times 12')
		place('x' => 300, 'y' => 100)
		pady 10
	end
	$widgets.push(inputfile, filelabel)
	$next_btn.state = "normal"
end
#-----------------------------------------------------------------------------------------------
def create_root_menu
	set_defaults # Create new options when returning to main menu
	$root.configure("title", "BaitsTools")
	aln2baits_btn = TkButton.new($root) do
		text "aln2baits"
		command 'subcommand_window("aln2baits")'
		place('x' => 100, 'y' => 100)
	end
	annot2baits_btn = TkButton.new($root) do
		text "annot2baits"
		command 'subcommand_window("annot2baits")'
		place('x' => 300, 'y' => 100)
	end
	bed2baits_btn = TkButton.new($root) do
		text "bed2baits"
		command 'subcommand_window("bed2baits")'
		place('x' => 500, 'y' => 100)
	end
	blast2baits_btn = TkButton.new($root) do
		text "blast2baits"
		command 'subcommand_window("blast2baits")'
		place('x' => 100, 'y' => 150)
	end
	checkbaits_btn = TkButton.new($root) do
		text "checkbaits"
		command 'subcommand_window("checkbaits")'
		place('x' => 300, 'y' => 150)
	end
	pyrad2baits_btn = TkButton.new($root) do
  		text "pyrad2baits"
  		command 'subcommand_window("pyrad2baits")'
		place('x' => 500, 'y' => 150)
	end
	stacks2baits_btn = TkButton.new($root) do
  		text "stacks2baits"
  		command 'subcommand_window("stacks2baits")'
		place('x' => 100, 'y' => 200)
	end
	tilebaits_btn = TkButton.new($root) do
		text "tilebaits"
		command 'subcommand_window("tilebaits")'
		place('x' => 300, 'y' => 200)
	end
	vcf2baits_btn = TkButton.new($root) do
		text "vcf2baits"
		command 'subcommand_window("vcf2baits")'
   		place('x' => 500, 'y' => 200)
	end
	$widgets = [aln2baits_btn, annot2baits_btn, bed2baits_btn, blast2baits_btn, checkbaits_btn, pyrad2baits_btn, stacks2baits_btn, tilebaits_btn, vcf2baits_btn]
	configure_buttons($widgets)
#	poonheli = TkLabel.new($root) do
#		image TkPhotoImage.new(:file => "~/baitstools/poonheli.gif")
#		place('height' => 118, 'width' => 118, 'x' => 341, 'y' => 400)
#	end
#	$widgets.push(poonheli)
	$labelVar.value = "Choose Subcommand"
end
#-----------------------------------------------------------------------------------------------
def taxa_window
	clear_widgets
	$labelVar.value = "vcf2baits Taxa Options and Previous Baits"
	taxafile = TkButton.new($root) do
		text 'Taxa file (optional)'
		command '$options.taxafile.value = Tk.getOpenFile'
		place('x' => 50, 'y' => 100)
	end
	taxalabel = TkLabel.new($root) do
		textvariable $options.taxafile
		font TkFont.new('times 12')
		place('x' => 300, 'y' => 100)
		pady 10
	end
	taxacount = TkLabel.new($root) do
		text 'Category variants (order: all,between,within) (comma-separated list)'
		font TkFont.new('times 20')
		place('x' => 50, 'y' => 150)
		pady 10
	end
	taxacountentry = TkEntry.new($root) do
		textvariable $options.taxacount
		borderwidth 5
		font TkFont.new('times 12')
		place('x' => 50, 'y' => 210)
		width 50
	end
	popcat = TkLabel.new($root) do
		text 'Population-specific variants (order of taxa file) (comma-separated list)'
		font TkFont.new('times 20')
		place('x' => 50, 'y' => 250)
		pady 10
	end
	popcatentry = TkEntry.new($root) do
		textvariable $options.popcategories
		borderwidth 5
		font TkFont.new('times 12')
		place('x' => 50, 'y' => 310)
		width 50
	end
	previousbaits = TkButton.new($root) do
		text 'Baits BED (optional)'
		command '$options.previousbaits.value = Tk.getOpenFile'
		place('x' => 50, 'y' => 360)
	end
	baitslabel = TkLabel.new($root) do
		textvariable $options.previousbaits
		font TkFont.new('times 12')
		place('x' => 300, 'y' => 360)
		pady 10
	end
	configure_buttons([taxafile, previousbaits])
	previousbaits.width = taxafile.width = 20
	$widgets.push(taxafile, taxalabel, taxacount, popcat, taxacountentry, popcatentry, previousbaits, baitslabel)
end
#-----------------------------------------------------------------------------------------------
def qc_window
	clear_widgets
	$labelVar.value = "Filtration Options"
	params = TkCheckButton.new($root) do
		variable $options.params
		text "Output parameters file"
		place('x' => 50, 'y' => 100)
		command 'update_params'
	end
	$no_lc = TkCheckButton.new($root) do
		variable $options.no_lc
		text "No linguistic complexity for parameters file"
		font TkFont.new('times 20')
		place('x' => 315, 'y' => 100)
		width 40
		pady 10
	end
	complete = TkCheckButton.new ($root) do
		variable $options.completebait
		text "Require complete length"
		place('x' => 50, 'y' => 150)
	end
	noNs = TkCheckButton.new($root) do
		variable $options.no_Ns
		text "No Ns"
		place('x' =>390, 'y' => 150)
	end
	collapse = TkCheckButton.new($root) do
		variable $options.collapse_ambiguities
		text "Collapse ambiguities"
		place('x' => 50, 'y' => 200)
	end
	maxmask = TkCheckButton.new($root) do
		variable $options.maxmask_filter
		text "Max mask%"
		place('x' => 390, 'y' => 200)
		command '$maxmaskentry.state == "disabled" ? $maxmaskentry.state = "normal" : $maxmaskentry.state = "disabled"'
	end
	$maxmaskentry = TkEntry.new($root) do
		textvariable $options.maxmask
		borderwidth 5
		font TkFont.new('times 12')
		place('x' => 580, 'y' => 210)
		width 10
	end
	mingc = TkCheckButton.new($root) do
		variable $options.mingc_filter
		text "Min GC%"
		place('x' => 50, 'y' => 250)
		command '$mingcentry.state == "disabled" ? $mingcentry.state = "normal" : $mingcentry.state = "disabled"'
	end
	$mingcentry = TkEntry.new($root) do
		textvariable $options.mingc
		borderwidth 5
		font TkFont.new('times 12')
		place('x' => 240, 'y' => 260)
		width 10
	end
	maxgc = TkCheckButton.new($root) do
		variable $options.maxgc_filter
		text "Max GC%"
		place('x' => 390, 'y' => 250)
		command '$maxgcentry.state == "disabled" ? $maxgcentry.state = "normal" : $maxgcentry.state = "disabled"'
	end
	$maxgcentry = TkEntry.new($root) do
		textvariable $options.maxgc
		borderwidth 5
		font TkFont.new('times 12')
		place('x' => 580, 'y' => 260)
		width 10
	end
	$mingcentry.state = "disabled" if $options.mingc_filter == 0
	$maxgcentry.state = "disabled" if $options.maxgc_filter == 0
	homopoly = TkCheckButton.new($root) do
		variable $options.maxhomopoly_filter
		text "Max homopolymer"
		place('x' => 50, 'y' => 300)
		command '$homopolyentry.state == "disabled" ? $homopolyentry.state = "normal" : $homopolyentry.state = "disabled"'
	end
	$homopolyentry = TkEntry.new($root) do
		textvariable $options.maxhomopoly
		borderwidth 5
		font TkFont.new('times 12')
		place('x' => 240, 'y' => 310)
		width 10
	end
	$lcentry = TkEntry.new($root) do
		textvariable $options.lc
		borderwidth 5
		font TkFont.new('times 12')
		place('x' => 580, 'y' => 310)
		width 10
	end
	minlc = TkCheckButton.new($root) do
		variable $options.lc_filter
		text "Min complexity"
		place('x' => 390, 'y' => 300)
		command 'update_minlc'
	end
	$homopolyentry.state = "disabled" if $options.maxhomopoly_filter == 0
	mint = TkCheckButton.new($root) do
		variable $options.mint_filter
		text "Min Tm (°C)"
		place('x' => 50, 'y' => 350)
		command 'update_mint'
	end
	$mintentry = TkEntry.new($root) do
		textvariable $options.mint
		borderwidth 5
		font TkFont.new('times 12')
		place('x' => 240, 'y' => 360)
		width 10
	end
	maxt = TkCheckButton.new($root) do
		variable $options.maxt_filter
		text "Max Tm (°C)"
		place('x' => 390, 'y' => 350)
		command 'update_maxt'
	end
	$maxtentry = TkEntry.new($root) do
		textvariable $options.maxt
		borderwidth 5
		font TkFont.new('times 12')
		place('x' => 580, 'y' => 360)
		width 10
	end
	$typelab = TkLabel.new($root) do
		text 'Hybridization'
		font TkFont.new('times 20')
		place('x' => 50, 'y' => 400)
		pady 10
	end
	$typeselect = Tk::Tile::Combobox.new($root) do
		textvariable $options.bait_type
		values ["RNA-DNA", "DNA-DNA", "RNA-RNA"]
		state "readonly"
		width 10
		height 3
		place('x' => 180, 'y' => 410)
	end
	$sodium = TkLabel.new($root) do
		text 'Sodium (M)'
		font TkFont.new('times 20')
		place('x' => 320, 'y' => 400)
		pady 10
	end
	$sodiumentry = TkEntry.new($root) do
		textvariable $options.na
		borderwidth 5
		font TkFont.new('times 12')
		place('x' => 450, 'y' => 410)
		width 10
	end
	$formamide = TkLabel.new($root) do
		text 'Formamide (M)'
		font TkFont.new('times 20')
		place('x' => 540, 'y' => 400)
		pady 10
	end
	$formamideentry = TkEntry.new($root) do
		textvariable $options.formamide
		borderwidth 5
		font TkFont.new('times 12')
		place('x' => 690, 'y' => 410)
		width 10
	end
	$maxmaskentry.state = "disabled" if $options.maxmask_filter == 0
	$mintentry.state = "disabled" if $options.mint_filter == 0
	$maxtentry.state = "disabled" if $options.maxt_filter == 0
	if $maxtentry.state == "disabled" && $mintentry.state == "disabled"
		$formamideentry.state = $formamide.state = "disabled"
		$sodiumentry.state = $sodium.state = "disabled"
		$typeselect.state = $typelab.state = "disabled"
	end	
	meanqual = TkCheckButton.new($root) do
		variable $options.meanqual_filter
		text "Min mean base Q"
		place('x' => 50, 'y' => 450)
		command '$meanqualentry.state == "disabled" ? $meanqualentry.state = "normal" : $meanqualentry.state = "disabled"'
	end
	$meanqualentry = TkEntry.new($root) do
		textvariable $options.meanqual
		borderwidth 5
		font TkFont.new('times 12')
		place('x' => 230, 'y' => 460)
		width 10
	end
	minqual = TkCheckButton.new($root) do
		variable $options.minqual_filter
		text "Min base Q"
		place('x' => 320, 'y' => 450)
		command '$minqualentry.state == "disabled" ? $minqualentry.state = "normal" : $minqualentry.state = "disabled"'
	end
	$minqualentry = TkEntry.new($root) do
		textvariable $options.minqual
		borderwidth 5
		font TkFont.new('times 12')
		place('x' => 450, 'y' => 460)
		width 10
	end
	fastascore = TkLabel.new($root) do
		text 'Assumed base Q'
		font TkFont.new('times 20')
		place('x' => 540, 'y' => 450)
		pady 10
	end
	$fastascoreentry = TkEntry.new($root) do
		textvariable $options.fasta_score
		borderwidth 5
		font TkFont.new('times 12')
		place('x' => 690, 'y' => 460)
		width 10
	end
	$meanqualentry.state = "disabled" if $options.meanqual_filter == 0
	$minqualentry.state = "disabled" if $options.minqual_filter == 0
	$lcentry.state = "disabled" if $options.lc_filter == 0
	$no_lc.state = "disabled" if $options.params == 0 or $options.lc_filter == 1
	$widgets.push(params, complete, noNs, collapse, maxmask, mingc, maxgc, homopoly, minlc, mint, maxt, meanqual, minqual)
	configure_buttons($widgets)
	$widgets.push($maxmaskentry, $mingcentry, $maxgcentry, $homopolyentry, $no_lc, $lcentry, $mintentry, $maxtentry, $typelab, $typeselect, $sodium, $sodiumentry, $formamide, $formamideentry, $meanqualentry, $minqualentry, $fastascoreentry, fastascore)
	params.width = complete.width = collapse.width = 20
	meanqual.width = homopoly.width = minlc.width = 15
end
#-----------------------------------------------------------------------------------------------
def update_params
	if $options.params == 1 && $options.lc_filter == 0
		$no_lc.state = "normal"
	elsif $options.params == 0
		$no_lc.state = "disabled"
		$options.no_lc.value = 0
	end
end
#-----------------------------------------------------------------------------------------------
def update_minlc
	$lcentry.state == "disabled" ? $lcentry.state = "normal" : $lcentry.state = "disabled"
	if $options.params == 1 && $options.lc_filter == 1
		$no_lc.state = "disabled"
		$options.no_lc.value = 0
	elsif $options.params == 1
		$no_lc.state = "normal"
	end	
end
#-----------------------------------------------------------------------------------------------
def update_maxt
	$maxtentry.state == "disabled" ? $maxtentry.state = "normal" : $maxtentry.state = "disabled"
	if $maxtentry.state == "disabled" && $mintentry.state == "disabled"
		$formamideentry.state = $formamide.state = "disabled"
		$sodiumentry.state = $sodium.state = "disabled"
		$typeselect.state = $typelab.state = "disabled"
	else
		$formamideentry.state = $formamide.state = "normal"
		$sodiumentry.state = $sodium.state = "normal"
		$typeselect.state = $typelab.state = "normal"
	end
end
#-----------------------------------------------------------------------------------------------
def update_mint
	$mintentry.state == "disabled" ? $mintentry.state = "normal" : $mintentry.state = "disabled"
	if $maxtentry.state == "disabled" && $mintentry.state == "disabled"
		$formamideentry.state = $formamide.state = "disabled"
		$sodiumentry.state = $sodium.state = "disabled"
		$typeselect.state = $typelab.state = "disabled"
	else
		$formamideentry.state = $formamide.state = "normal"
		$sodiumentry.state = $sodium.state = "normal"
		$typeselect.state = $typelab.state = "normal"
	end
end
#-----------------------------------------------------------------------------------------------
def go_back
	case $labelVar.value
	when "Filtration Options"
		if $options.algorithm == "vcf2baits" && $options.every != 1
			taxa_window
		else
			subcommand_window($options.algorithm)
		end
	when "vcf2baits Taxa Options and Previous Baits"
		subcommand_window($options.algorithm)
	when "General Options"
		qc_window
	else
		clear_widgets
		create_root_menu
		$back_btn.state = "disabled"
		$next_btn.state = "disabled"
	end
end
#-----------------------------------------------------------------------------------------------
def go_forward
	case $labelVar.value
	when "aln2baits Options", "annot2baits Options", "bed2baits Options", "blast2baits Options", "checkbaits Options", "tilebaits Options"
		if $options.infile == ""
			Tk::messageBox :message => 'Please specify an input file.'
		elsif ($options.algorithm == "annot2baits" or $options.algorithm == "bed2baits" or $options.algorithm == "blast2baits") && $options.refseq == ""
			Tk::messageBox :message => 'Please specify a reference sequence.'
		elsif ($options.algorithm == "annot2baits" or $options.algorithm == "bed2baits" or $options.algorithm == "blast2baits") && $options.pad < 0
			Tk::messageBox :message => 'Pad length cannot be less than 0.'
		elsif $options.baitlength < 1
			Tk::messageBox :message => 'Bait length must be greater than 0.'
		elsif $options.tileoffset < 1 && $options.algorithm != "checkbaits"
			Tk::messageBox :message => 'Tiling offset must be greater than 0.'
		elsif $options.percid < 0.0 or $options.percid > 100.0
			Tk::messageBox :message => 'Percent identity must be between 0.0 and 100.0.'
		elsif $options.blastlen < 1
			Tk::messageBox :message => 'Minimum BLAST hit length must be greater than 0.'
		elsif $options.evalue < 0 && $options.evalue_filter == 1
			Tk::messageBox :message => 'Maximum E-value must be at least 0.'
		else
			qc_window
		end
	when "pyrad2baits Options"
		if $options.infile == ""
			Tk::messageBox :message => 'Please specify an input file.'
		elsif $options.baitlength < 1
			Tk::messageBox :message => 'Bait length must be greater than 0.'
		elsif $options.tileoffset < 1
			Tk::messageBox :message => 'Tiling offset must be greater than 0.'
		elsif $options.minind < 1
			Tk::messageBox :message => 'Minimum individuals must be greater than 0.'
		elsif $options.strategy != "alignment"
			if $options.totalsnps < 1
				Tk::messageBox :message => 'The total number of variants must be greater than 0.'
			elsif $options.maxsnps < 1
				Tk::messageBox :message => 'The maximum number of SNPs per locus must be greater than 0.'
			elsif $options.distance < 1
				Tk::messageBox :message => 'The minimum distance between variants must be greater than 0.'
			elsif $options.lenbef < 0
				Tk::messageBox :message => 'The number of bait bases before the variant must be at least 0.'
			elsif $options.tileoffset < 1
				Tk::messageBox :message => 'Tiling offset cannot be less than 1.'
			elsif $options.tiledepth < 1
				Tk::messageBox :message => 'Tiles per SNP cannot be less than 1.'
			else
				qc_window
			end
		else
			qc_window
		end
	when "vcf2baits Options","stacks2baits Options"
		if $options.infile == ""
			Tk::messageBox :message => 'Please specify an input file.'
		elsif $options.refseq == "" && $options.no_baits == 0
			Tk::messageBox :message => 'Please specify a reference sequence.'
		elsif $options.totalsnps.to_i < 1 && $options.every == 0
			Tk::messageBox :message => 'The total number of variants must be greater than 0.'		
		elsif $options.scale == 0 && $options.maxsnps.to_i < 1 && $options.every == 0
			Tk::messageBox :message => 'The maximum number of variants per contig must be greater than 0.'
		elsif $options.distance.to_i < 1 && $options.every == 0
			Tk::messageBox :message => 'The minimum distance between variants must be greater than 0.'
		elsif $options.baitlength < 1
			Tk::messageBox :message => 'Bait length must be greater than 0.'
		elsif $options.lenbef.to_i < 0
			Tk::messageBox :message => 'The number of bait bases before the variant must be at least 0.'
		elsif $options.tileoffset.to_i < 1
			Tk::messageBox :message => 'Tiling offset cannot be less than 1.'
		elsif $options.tiledepth.to_i < 1
			Tk::messageBox :message => 'Tiles per SNP cannot be less than 1.'
		elsif ($options.alpha.to_f < 0.0 or $options.alpha.to_f > 1.0) && $options.hwe == 1
			Tk::messageBox :message => 'Alpha must be between 0.0 and 1.0.'
		elsif $labelVar.value == "vcf2baits Options" and $options.every != 1
			taxa_window
		else
			qc_window
		end
	when "vcf2baits Taxa Options and Previous Baits"
		data_fails = false
		if $options.taxafile != ""
			taxacnt = $options.taxacount.value.split(",").map { |x| x.to_i }
			if taxacnt.size != 3
				Tk::messageBox :message => 'Three values needed for all categories (AllPopulations,BetweenPopulations,WithinPopulations)'
				data_fails = true
			elsif taxacnt[0] < 0
				Tk::messageBox :message => 'The number of variants variable across populations must be at least 0.'
				data_fails = true
			elsif taxacnt[1] < 0
				Tk::messageBox :message => 'The number of variants variable between populations must be at least 0.'
				data_fails = true
			elsif taxacnt[2] < 0
				Tk::messageBox :message => 'The number of variants variable within populations must be at least 0.'
				data_fails = true
			elsif taxacnt[0] + taxacnt[1] + taxacnt[2] < 1  # Ruby 2.0 does not have Array#sum method
				Tk::messageBox :message => 'The total number of variants must be greater than 0.'
				data_fails = true
			elsif $options.popcategories != "" # Check that population values are ok
				popcnt = $options.popcategories.value.split(",").map { |x| x.to_i }
				for pop in popcnt
					if pop < 0 or pop > taxacnt[2]
						data_fails = true
						Tk::messageBox :message => 'Population-specific variant numbers must be between 0 and the total number of within-population variants.'
						break
					end
				end
			end
		end
		qc_window unless data_fails
	when "Filtration Options"
		data_fails = false
		if $options.maxmask_filter == 1
			if $options.maxmask < 0.0 or $options.maxmask > 100.0
				Tk::messageBox :message => 'Max mask% must be between 0 and 100.'
				data_fails = true
				return
			end
		elsif $options.mingc_filter == 1
			if $options.mingc < 0.0 or $options.mingc > 100.0
				Tk::messageBox :message => 'Min GC% must be between 0 and 100.'
				data_fails = true
				return
			elsif $options.maxgc_filter == 1
				if $options.maxgc.to_f < $options.mingc.to_f
					Tk::messageBox :message => 'Max GC% must be greater than min GC%.'
					data_fails = true
					return
				elsif $options.maxgc > 100.0
					Tk::messageBox :message => 'Max GC% must between min GC% and 100.'
					data_fails = true
					return
				end
			end
		elsif $options.maxgc_filter == 1
			if $options.maxgc < 0.0 or $options.maxgc > 100.0
				Tk::messageBox :message => 'Max GC% must be between 0 and 100.'
				data_fails = true
				return
			end
		elsif $options.maxhomopoly_filter == 1
			if $options.maxhomopoly < 1
				Tk::messageBox :message => 'Maximum homopolymer length must be greater than 0.'
				data_fails = true
				return
			end
		elsif $options.lc_filter == 1
			if $options.lc < 0.0 or $options.lc > 1.0
				Tk::messageBox :message => 'Complexity must be between 0.0 and 1.0'
				data_fails = true
				return
			end
		end
		if $options.mint_filter && $options.maxt_filter && $options.maxt.to_f < $options.mint.to_f
			Tk::messageBox :message => 'Max Tm cannot be less than min Tm.'
			data_fails = true
			return
		end
		if ($options.mint_filter or $options.maxt_filter)
			if $options.na.to_f < 0.0
				Tk::messageBox :message => 'Sodium concentration cannot be negative.'
				data_fails = true
				return
			elsif $options.formamide.to_f < 0.0
				Tk::messageBox :message => 'Formamide concentration cannot be negative.'
				data_fails = true
				return
			end
		end
		if $options.meanqual_filter == 1
			if $options.meanqual < 0.0 or $options.meanqual > 93.0
				Tk::messageBox :message => 'Mean quality must be between 0.0 and 93.0.'
				data_fails = true
				return
			end
		end
		if $options.minqual_filter == 1
			if $options.minqual < 0 or $options.minqual > 93
				Tk::messageBox :message => 'Minimum quality must be between 0 and 93.'
				data_fails = true
				return
			end
		end
		if $options.fasta_score < 0 or $options.fasta_score > 93
			Tk::messageBox :message => 'Assumed quality must be between 0 and 93.'
			data_fails = true
			return
		end
		general_window if !data_fails
	when "General Options"
		if $options.threads < 1
			Tk::messageBox :message => 'Threads must be greater than 0.'
		else
			start_baitstools
		end
	end
end
#-----------------------------------------------------------------------------------------------
def configure_buttons(buttons)
	for button in buttons
		button.borderwidth = 5
  		button.state = "normal"
  		button.cursor = "arrow"
  		button.font = TkFont.new('times 20')
  		button.relief = "groove"
  		button.padx = 15
  		button.pady = 10
  		button.width = 10
  	end
end
#-----------------------------------------------------------------------------------------------
def clear_widgets
	for widget in $widgets
		widget.destroy()
	end
	$widgets = []
end
#-----------------------------------------------------------------------------------------------
def set_defaults
	$options = OpenStruct.new
	$options.algorithm = "" # BaitsTools subcommand
	# Algorithm Parameters
	$options.infile = TkVariable.new("") # Primary input file
	$options.refseq = TkVariable.new("") # Reference sequence file
	$options.taxafile = TkVariable.new("") # Taxa assignment file
	$options.taxacount = TkVariable.new("") # Variant category counts
	$options.popcategories = TkVariable.new("") # Population-specific variant counts
	$options.previousbaits = TkVariable.new("") # Previous baits BED file
	$options.baitlength = TkVariable.new(120) # Bait length
	$options.tileoffset = TkVariable.new(60) # Offset between tiled baits
	$options.bait_type = TkVariable.new("RNA-DNA") # Hybridization type
	$options.haplodef = TkVariable.new("haplotype") # Haplotype definition for aln2baits
	$options.list_format = TkVariable.new("BED") # Interval list file format
	$options.features = TkVariable.new("") # Desired features in comma-separated list
	$options.pad = TkVariable.new(0) # BP to pad ends of extracted regions
	$options.totalsnps = TkVariable.new(30000) # Maximum requested SNPs
	$options.no_baits = TkVariable.new(0) # Flag to omit generating baits
	$options.every = TkVariable.new(0) # Flag that baits will be generated from every SNP
	$options.alt_alleles = TkVariable.new(0) # Flag to apply alternate alleles
	$options.scale = TkVariable.new(0) # Flag to scale SNPs per contig by contig length
	$options.varqual_filter = TkVariable.new(0) # Flag to determine whether to filter vcf variants by QUAL scores
	$options.varqual = TkVariable.new(30) # Minimum vcf variant QUAL score
	$options.minind = TkVariable.new(1) # Minimum individuals to include locus
	$options.strategy = TkVariable.new("alignment") # Strategy for LOCI files
	$options.uncollapsed_ref = TkVariable.new(0) # Flag to output uncollapsed reference sequence
	$options.lenbef = TkVariable.new(60) # Length before SNP in bait
	$options.tiledepth = TkVariable.new(1) # Tiling depth	
	$options.sort = TkVariable.new(0) # Flag to sort stack2baits SNPs by between/within population variation
	$options.hwe = TkVariable.new(0) # Flag to sort stacks2baits SNPs by Hardy-Weinberg Equilibrium
	$options.alpha = TkVariable.new(0.05) # HWE test alpha value
	$options.percid = TkVariable.new(0.0) # Percent identity for BLAST hit
	$options.blastlen = TkVariable.new(1) # Minimum length for BLAST hit
	$options.evalue = TkVariable.new(0.1) # Maximum E-value to include BLAST hit
	$options.evalue_filter = TkVariable.new(0) # Flag to filter BLAST hits by E-value
	$options.inbed = TkVariable.new("") # Optional corresponding BED file for checkbaits
	# Filtration Parameters
	$options.params = TkVariable.new(0) # Flag to output filtration parameters
	$options.no_lc = TkVariable.new(0) # Flag to disable linguistic complexity calculation for filtration
	$options.completebait = TkVariable.new(0) # Flag to filter by complete bait length	
	$options.no_Ns = TkVariable.new(0) # Flag to omit bait sequences with Ns
	$options.noaddenda = TkVariable.new(0) # Flag to omit sequence addenda from bait parameter calculations
	$options.collapse_ambiguities = TkVariable.new(0) # Flag to collapse ambiguities to a single nucleotide
	$options.maxmask_filter = TkVariable.new(0) # Flag to filter by % masked sequence
	$options.maxmask = TkVariable.new(25.0) # Max % masked sequence
	$options.maxhomopoly_filter = TkVariable.new(0) # Flag to filter by homopolymer length
	$options.maxhomopoly = TkVariable.new(4) # Max homopolymer length
	$options.lc_filter = TkVariable.new(0) # Flag to filter by sequence complexity
	$options.lc = TkVariable.new(0.9) # Minimum complexity
	$options.mingc_filter = TkVariable.new(0) # Flag to filter by minimum gc content
	$options.mingc = TkVariable.new(30.0) # Minimum GC content
	$options.maxgc_filter = TkVariable.new(0) # Flag to filter by maximum gc content
	$options.maxgc = TkVariable.new(50.0) # Maximum GC content
	$options.mint_filter = TkVariable.new(0) # Flag to filter by minimum melting temperature
	$options.mint = TkVariable.new(0.0) # Minimum melting temperature
	$options.maxt_filter = TkVariable.new(0) # Flag to filter by maximum melting temperature
	$options.maxt = TkVariable.new(120.0) # Maximum melting temperature
	$options.na = TkVariable.new(0.9) # Sodium concentration
	$options.formamide = TkVariable.new(0.0) # Formamide concentration
	$options.meanqual_filter = TkVariable.new(0) # Flag to filter by minimum mean base quality
	$options.meanqual = TkVariable.new(20.0) # Minimum mean base quality
	$options.minqual_filter = TkVariable.new(0) # Flag to filter by minimum base quality
	$options.minqual = TkVariable.new(10) # Minimum base quality
	$options.fasta_score = TkVariable.new(0) # Asssumed base quality score for FASTA sequences
	$options.phred64 = TkVariable.new(0) # Use phred64 quality encoding
	$options.gzip = TkVariable.new(0) # Gzip output
	# General Parameters
	$options.outprefix = TkVariable.new("out") # Output prefix
	$options.outdir = TkVariable.new(File.expand_path("./")) # Output directory
	$options.coords = TkVariable.new(0) # Flag to output BED table of baits
	$options.rbed = TkVariable.new(0) # Flag to output relative BED table of baits
	$options.shuffle = TkVariable.new(0) # Flag to shuffle baits if extend beyond contig
	$options.log = TkVariable.new(0) # Flag to output detailed log
	$options.ncbi = TkVariable.new(0) # Flag whether FASTA/FASTQ headers include NCBI-style descriptors
	$options.rna = TkVariable.new(0) # Flag whether baits are output as RNA
	$options.rc = TkVariable.new(0) # Flag to output reverse complement baits
	$options.gaps = TkVariable.new("include") # Flag to omit bait sequences with gaps
	$options.threeprime = TkVariable.new("") # 3' sequence addendum
	$options.fiveprime = TkVariable.new("") # 5' sequence addendum
	$options.fillin = TkVariable.new("") # Fill-in sequence motif
	$options.threads = TkVariable.new(1) # Number of threads
	$options.rng = TkVariable.new(srand) # RNG seed
end
#-----------------------------------------------------------------------------------------------
$root = TkRoot.new { title "BaitsTools" }
$root.geometry = '800x600'
set_defaults
$labelVar = TkVariable.new
$label = TkLabel.new($root) do
	textvariable $labelVar
	borderwidth 5
	font TkFont.new('times 30 bold')
	pack("side" => "top",  "padx"=> "50", "pady"=> "10")
end
sepline = Tk::Tile::Separator.new($root) do
	orient 'horizontal'
	place('width' => 760, 'x' => 20, 'y' => 55)
end
create_root_menu
exit_btn = TkButton.new($root) do
	text "Quit"
	command 'exit'
	place('x' => 10, 'y' => 560)
end
$back_btn = TkButton.new($root) do
	text "Back"
	command 'go_back'
	place('x' => 660, 'y' => 560)
end
$next_btn = TkButton.new($root) do 
	text 'Next'
	command 'go_forward'
	place('x' => 660, 'y' => 520)
end
credit = TkLabel.new($root) do
	text "Michael G. Campana, 2017-2021\nSmithsonian Conservation Biology Institute"
	borderwidth 5
	font TkFont.new('times 12')
	pack("side" => "bottom",  "padx"=> "50", "pady"=> "10")
	pady 10
	padx 10
end
configure_buttons([exit_btn, $back_btn, $next_btn])
$back_btn.state = "disabled"
$next_btn.state = "disabled"
Tk.mainloop
