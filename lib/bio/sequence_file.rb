require 'bio/genomics/assembly'
require 'bio/genomics/data_set'
require 'stats'

module Bio
  # Base class for getting genomic sequences from 2bit/Fasta/etc files
  # Loads the sequences as Bio::Sequence::NA
  class SequenceFile < Genomics::DataSet
    attr_reader :assembly
    
    def initialize(filename, assembly = nil)
      @data_file = File.expand_path(filename)
      @assembly = assembly
    end
    
    # Open, optionally with a block
    def self.open(filename, assembly = nil)
      if block_given?
        yield self.new(filename, assembly)
      else
        return self.new(filename, assembly)
      end
    end
    
    # Check if this genome includes the specified interval
    alias :super_include? :include?
    def include?(chr, start = nil, stop = nil)
      if start.nil?
        super_include?(chr)
      elsif stop.nil?
        super_include?(chr) and start >= 1 and start <= self[chr]
      else
        super_include?(chr) and start >= 1 and stop <= self[chr]
      end
    end
    
    # Get a specific stretch of sequence
    def sequence(chr, start = nil, stop = nil)
      raise SequenceFileError, "Do not now how to query for sequences from #{File.basename(@data_file)}!"
    end
    
    # Alias for sequence
    def query(chr, start = nil, stop = nil)
      sequence(chr, start, stop)
    end
    
    # Return the total number of base pairs in the file
    # Genome#length will return the number of chromosomes
    def num_bases
      self.map { |chr,n| n }.sum
    end
    
    # Reduce this genome to an Assembly object (just chromosome id's and their lengths)
    def to_assembly(name = @assembly)
      a = Assembly.new(name, nil)
      self.each { |chr, n| a[chr] = n }
      return a
    end

    def to_s
      str = "Sequence File #{@assembly}: containing #{num_bases} base pairs\n"
      self.each do |chr, n|
        str += "\tChromosome #{chr} (length: #{n})\n"
      end
      
      return str
    end
  end

  class SequenceFileError < StandardError
  end
end
