Gem::Specification.new do |s|
	s.name = 'baitstools'
	s.version = '1.8.1'
	s.required_ruby_version = '>= 2.4.1'
	s.date = '2023-08-04'
	s.summary = 'BaitsTools: Software for hybridization capture bait design'
	s.description = 'Software for hybridization capture bait design'
	s.authors = ["Michael G. Campana"]
	s.email = 'campanam@si.edu'
	s.metadata["github_repo"] = "https://github.com/campanam/BaitsTools"
	s.files = ["lib/aln2baits.rb","lib/annot2baits.rb","lib/baitslib.rb","lib/bed2baits.rb",
				"lib/blast2baits.rb","lib/checkbaits.rb","lib/pyrad2baits.rb","lib/stacks2baits.rb",
				"lib/tilebaits.rb","lib/vcf2baits.rb"]
	if ARGV.include? 'gui'
		s.executables = ["baitstools","baitstoolsgui"]
		s.add_runtime_dependency 'tk','0.4.0'
		s.requirements << 'macOS Mojave (or higher): ActiveTcl-8.6'
	else
		s.executables = ["baitstools"]
	end
	s.homepage = 'https://github.com/campanam/BaitsTools'
	s.license = 'Nonstandard'
	s.add_runtime_dependency 'shell', '0.8.1'
end
