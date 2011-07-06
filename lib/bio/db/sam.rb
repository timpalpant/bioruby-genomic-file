require 'bio/db/entry_file'
require 'bio/db/read_file'
require 'bio/genomics/read'
require 'bio/utils/samtools'
require 'stringio'

# An entry in a SAM file
# From the specification, see: http://samtools.sourceforge.net/SAM-1.4.pdf
# Also helpful: http://chagall.med.cornell.edu/NGScourse/SAM.pdf
# TODO: Better handling of empty (*) values
class SAMEntry < Read
  attr_accessor :qname, :flag, :rname, :mapq, :cigar, :rnext, :pnext, :tlen
  
  def self.parse(line)
    begin
      record = line.chomp.split("\t")
      raise SAMError, "Invalid SAM Entry: SAM/BAM alignments must have 11 columns" if record.length < 11
    
      entry = self.new
      entry.qname = record[0]
      entry.flag = record[1].to_i
      entry.chr = record[2]
      entry.mapq = record[4].to_i
      entry.cigar = record[5]
      entry.rnext = record[6]
      entry.pnext = record[7].to_i
      entry.tlen = record[8].to_i
      entry.seq = Bio::Sequence::NA.new(record[9])
      entry.qual = record[10]
    
      extend = (entry.tlen == 0) ? entry.seq.length : entry.tlen.abs
    
      # According to the SAM specification, POS is the leftmost base
      # so 5' end for forward-mapping reads and
      # 3' end for reverse-complement mapping reads
      # Adjust start/stop appropriately
      if entry.watson?
        entry.start = record[3].to_i
        entry.stop = entry.start + extend - 1
      else
        entry.start = record[3].to_i + entry.seq.length - 1
        entry.stop = entry.start - extend + 1
        entry.seq = entry.seq.reverse_complement
        entry.qual = entry.qual.reverse
      end

      return entry
    rescue
      raise SAMError, "Not a valid SAM Entry"
    end
  end
  
  def rname
    @chr
  end

  def pos
    low
  end

  # FLAGS: See SAM specification
  def paired?
    (@flag & 0x0001) != 0
  end

  def unpaired?
    not paired?
  end

  def single?
    not paired?
  end

  def proper_pair?
    (@flag & 0x0002) != 0
  end

  def mapped?
    (@flag & 0x0004) == 0
  end

  def unmapped?
    not mapped?
  end

  def mate_mapped?
    (@flag & 0x0008) == 0
  end

  def mate_unmapped?
    not mate_mapped?
  end

  def watson?
    (@flag & 0x0010) == 0
  end

  def crick?
    not watson?
  end

  def mate_watson?
    (@flag & 0x0020) == 0
  end

  def mate_crick?
    not mate_watson?
  end

  def first?
    (@flag & 0x0040) != 0
  end

  def second?
    (@flag & 0x0080) != 0
  end

  def primary?
    (@flag & 0x0100) == 0
  end

  def failed?
    (@flag & 0x0200) != 0
  end

  def duplicate?
    (@flag & 0x0400) != 0
  end
end


# Access SAM files through stream operations
# Best option for iterating over large SAM files
class SAMFile < TextEntryFile
  extend ReadFile
  include ReadFileMethods
  
  CHR_COL = 3
  START_COL = 4
  END_COL = 4
  
  def initialize(filename)
    super(filename, CHR_COL, START_COL, END_COL)
  end
  
  private
  
  def parse(line)
    SAMEntry.parse(line)
  end
end

# Access BAM files through stream operations using samtools
class BAMFile < BinaryEntryFile
  extend ReadFile
  include ReadFileMethods
  
  def initialize(filename)
    # BAM index is the filename plus the extension .bai
    super(filename, filename+'.bai')
  end
  
  # Count the number of alignments in a given lookup
  def count(chr = nil, start = nil, stop = nil)
    index() if not indexed?
    SAMTools.count(@data_file, chr, start, stop)
  end
  
  # Use samtools to get the chromosomes efficiently
  def chromosomes
    # Cache for performance
    if @chromosomes.nil?
      index() if not indexed?
      @chromosomes = Array.new
      SAMTools.idxstats(@data_file).each do |line| 
        entry = line.chomp.split("\t")
        chr = entry[0]
        bases = entry[1].to_i
        mapped = entry[2].to_i
        unmapped = entry[3].to_i
        @chromosomes << chr if (mapped+unmapped) > 0
      end
    end
    
    return @chromosomes
  end

  private
  
  # Index sorted alignment for fast random access. Index file <aln.bam>.bai will be created.
  def index
    # Block if this file is being indexed by another thread/process
    if indexing?
      while indexing?
        sleep 1
      end
    else
      begin
        self.indexing = true
        SAMTools.index(@data_file)
      rescue
        raise SAMError, "Error generating index for BAM file #{File.basename(@data_file)}!"
      ensure
        self.indexing = false
      end
    end
  end
  
  def query_lines(chr = nil, start = nil, stop = nil, &block)
    index() if not indexed?
    SAMTools.view(@data_file, chr, start, stop, &block)
  end
  
  def parse(line)
    SAMEntry.parse(line)
  end
end

class SAMError < EntryFileError
end
