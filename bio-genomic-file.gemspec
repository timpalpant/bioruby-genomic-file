# Generated by jeweler
# DO NOT EDIT THIS FILE DIRECTLY
# Instead, edit Jeweler::Tasks in Rakefile, and run 'rake gemspec'
# -*- encoding: utf-8 -*-

Gem::Specification.new do |s|
  s.name = %q{bio-genomic-file}
  s.version = "0.3.0"

  s.required_rubygems_version = Gem::Requirement.new(">= 0") if s.respond_to? :required_rubygems_version=
  s.authors = ["timpalpant"]
  s.date = %q{2011-08-31}
  s.description = %q{Allows iteration and querying of common genomic file formats, including Bed, BedGraph, GFF, Wiggle, and SAM/BAM}
  s.email = %q{tim@palpant.us}
  s.extra_rdoc_files = [
    "LICENSE.txt",
    "README.rdoc"
  ]
  s.files = [
    ".document",
    ".rspec",
    "Gemfile",
    "Gemfile.lock",
    "LICENSE.txt",
    "README.rdoc",
    "Rakefile",
    "VERSION",
    "bio-genomic-file.gemspec",
    "lib/.DS_Store",
    "lib/bio-genomic-file.rb",
    "lib/bio/.DS_Store",
    "lib/bio/affy.rb",
    "lib/bio/bed.rb",
    "lib/bio/bedgraph.rb",
    "lib/bio/entry_file.rb",
    "lib/bio/entry_file_sniffer.rb",
    "lib/bio/genomics/assembly.rb",
    "lib/bio/genomics/contig.rb",
    "lib/bio/genomics/data_set.rb",
    "lib/bio/genomics/index_error.rb",
    "lib/bio/genomics/interval.rb",
    "lib/bio/genomics/nucleosome.rb",
    "lib/bio/genomics/read.rb",
    "lib/bio/genomics/spot.rb",
    "lib/bio/gff.rb",
    "lib/bio/macs.rb",
    "lib/bio/nuke_calls.rb",
    "lib/bio/read_file.rb",
    "lib/bio/sam.rb",
    "lib/bio/sequence_file.rb",
    "lib/bio/spot_file.rb",
    "lib/bio/twobit.rb",
    "lib/bio/utils/samtools.rb",
    "lib/bio/utils/tabix.rb",
    "lib/bio/utils/ucsc.rb",
    "lib/bio/wig.rb",
    "lib/sparse_array.rb",
    "lib/stats.rb",
    "lib/utils/fixed_precision.rb",
    "lib/utils/numeric.rb",
    "lib/utils/parallelizer.rb",
    "lib/utils/unix.rb",
    "spec/bio/affy_spec.rb",
    "spec/bio/assembly_spec.rb",
    "spec/bio/bed_spec.rb",
    "spec/bio/bedgraph_spec.rb",
    "spec/bio/entry_file_sniffer_spec.rb",
    "spec/bio/genomics/contig_spec.rb",
    "spec/bio/genomics/data_set_spec.rb",
    "spec/bio/genomics/interval_spec.rb",
    "spec/bio/genomics/nucleosome_spec.rb",
    "spec/bio/genomics/read_spec.rb",
    "spec/bio/genomics/spot_spec.rb",
    "spec/bio/gff_spec.rb",
    "spec/bio/macs_spec.rb",
    "spec/bio/nuke_calls_spec.rb",
    "spec/bio/sam_spec.rb",
    "spec/bio/twobit_spec.rb",
    "spec/bio/utils/samtools_spec.rb",
    "spec/bio/utils/tabix_spec.rb",
    "spec/bio/utils/ucsc_spec.rb",
    "spec/bio/wig_spec.rb",
    "spec/fixtures/test.2bit",
    "spec/fixtures/test.affy",
    "spec/fixtures/test.bam",
    "spec/fixtures/test.bed",
    "spec/fixtures/test.bedGraph",
    "spec/fixtures/test.bw",
    "spec/fixtures/test.bw.tmp",
    "spec/fixtures/test.fa",
    "spec/fixtures/test.file1",
    "spec/fixtures/test.file2",
    "spec/fixtures/test.file3",
    "spec/fixtures/test.gff",
    "spec/fixtures/test.len",
    "spec/fixtures/test.macs",
    "spec/fixtures/test.nukes",
    "spec/fixtures/test.sam",
    "spec/fixtures/test.wig",
    "spec/sparse_array_spec.rb",
    "spec/spec_helper.rb",
    "spec/stats_spec.rb",
    "spec/utils/fixed_precision_spec.rb",
    "spec/utils/numeric_spec.rb",
    "spec/utils/parallelizer_spec.rb",
    "spec/utils/unix_spec.rb"
  ]
  s.homepage = %q{http://github.com/timpalpant/bio-genomic-file}
  s.licenses = ["MIT"]
  s.require_paths = ["lib"]
  s.rubygems_version = %q{1.3.7}
  s.summary = %q{Work with common genomic file formats in a consistent, simple and performant fashion}

  if s.respond_to? :specification_version then
    current_version = Gem::Specification::CURRENT_SPECIFICATION_VERSION
    s.specification_version = 3

    if Gem::Version.new(Gem::VERSION) >= Gem::Version.new('1.2.0') then
      s.add_runtime_dependency(%q<bio>, [">= 1.4.0"])
      s.add_runtime_dependency(%q<parallel>, ["~> 0.5.5"])
      s.add_development_dependency(%q<rspec>, ["~> 2.6.0"])
      s.add_development_dependency(%q<bundler>, ["~> 1.0.0"])
      s.add_development_dependency(%q<jeweler>, ["~> 1.6.2"])
      s.add_development_dependency(%q<simplecov>, ["~> 0.4.2"])
    else
      s.add_dependency(%q<bio>, [">= 1.4.0"])
      s.add_dependency(%q<parallel>, ["~> 0.5.5"])
      s.add_dependency(%q<rspec>, ["~> 2.6.0"])
      s.add_dependency(%q<bundler>, ["~> 1.0.0"])
      s.add_dependency(%q<jeweler>, ["~> 1.6.2"])
      s.add_dependency(%q<simplecov>, ["~> 0.4.2"])
    end
  else
    s.add_dependency(%q<bio>, [">= 1.4.0"])
    s.add_dependency(%q<parallel>, ["~> 0.5.5"])
    s.add_dependency(%q<rspec>, ["~> 2.6.0"])
    s.add_dependency(%q<bundler>, ["~> 1.0.0"])
    s.add_dependency(%q<jeweler>, ["~> 1.6.2"])
    s.add_dependency(%q<simplecov>, ["~> 0.4.2"])
  end
end

