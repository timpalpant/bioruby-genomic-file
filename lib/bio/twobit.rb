require 'bio/sequence_file'
require 'bio/utils/ucsc'

module Bio
  # Get genomic sequences from UCSC 2bit files
  class TwoBit < SequenceFile
    
    # Initialize a 2bit file with a genomic reference sequence
    def initialize(filename, assembly = File.basename(filename, '.2bit'))
      super(filename, assembly)
      Utils::UCSC.twobit_info(@data_file).each { |chr,n| self[chr] = n }
    end
    
    # Get a specific stretch of sequence
    def sequence(chr, start = nil, stop = nil)
      raise TwoBitError, "Genome does not include sequences for the range #{chr}:#{start}-#{stop}" if not include?(chr, start, stop)
      Utils::UCSC.twobit_to_fa(@data_file, chr, start, stop)
    end
  end

  class TwoBitError < SequenceFileError
  end
end
