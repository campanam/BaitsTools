Gem::Specification.new do |s|
	s.name = 'baitstools'
	s.version = '1.7.7'
	s.required_ruby_version = '>= 2.4.1'
	s.date = '2023-07-31'
	s.summary = 'BaitsTools: Software for hybridization capture bait design'
	s.description = 'Software for hybridization capture bait design'
	s.authors = ["Michael G. Campana"]
	s.email = 'campanam@si.edu'
	s.files = ["lib/aln2baits.rb","lib/annot2baits.rb","lib/baitslib.rb","lib/bed2baits.rb",
				"lib/blast2baits.rb","lib/checkbaits.rb","lib/pyrad2baits.rb","lib/stacks2baits.rb",
				"lib/tilebaits.rb","lib/vcf2baits.rb"]
	s.executables = ["baitstools","baitstoolsgui"]
	s.homepage = 'https://github.com/campanam/BaitsTools'
	s.license = 'Nonstandard'
	s.add_runtime_dependency 'tk','0.4.0'
	s.add_runtime_dependency 'shell', '0.8.1'
	s.requirements << 'macOS Mojave (or higher): ActiveTcl-8.6'
end
