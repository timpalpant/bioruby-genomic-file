require 'enumerator'
require 'set'
require 'fileutils'
require 'utils/unix'
require 'bio/utils/tabix'
require 'bio/entry_file_sniffer'

##
# A line-oriented data file
# base class for Bed, BedGraph, GFF, SAM, etc.
# Idea is that it should behave somewhat like a File object,
# but that relevant parsed Entry objects (BedEntry, SAMEntry, etc.) 
# are returned instead of lines
#
# In general, caching is kept to a minimum to avoid rampant memory usage
# and indexing is maximized. If you can benefit from caching, query for results
# and cache them where appropriate.
##

module Bio
  class EntryFile
    include Enumerable
    
    # Maximum number of entries to pre-cache
    CACHE_SIZE = 10_000

    def initialize(filename, index_file = nil)
      @data_file = File.expand_path(filename)
      raise EntryFileError, "Cannot find data file #{File.basename(@data_file)}. Does it exist?" if not File.exist?(@data_file)
      @index_file = File.expand_path(index_file) unless index_file.nil?
      @indexing = @data_file + '.indexing'
    end
    
    # Perform any additional cleanup operations (deleting indexes, etc.)
    def close
      File.delete(@index_file) if indexed?
      File.delete(@indexing) if indexing?
    end
    
    # Open the EntryFile (optionally with a block)
    def self.open(filename)
      entry_file = self.new(filename)
      
      if block_given?
        yield entry_file
        entry_file.close
      else
        return entry_file
      end
    end

    # Autodetect a file format
    def self.autodetect(filename, &block)
      clazz = EntryFileSniffer.sniff(filename)
      clazz.open(filename, &block)
    end
    
    # Iterate over each of the entries in an EntryFile
    def self.foreach(filename, chr = nil, start = nil, stop = nil)
      self.open(filename) do |f|
        f.each(chr, start, stop) { |entry| yield entry }
      end
    end
    
    # Make entry files enumerable
    def each(chr = nil, start = nil, stop = nil)
      skipped = 0
      
      query_lines(chr, start, stop) do |line|
        # Skip comment and track lines
        next if line.start_with?('#','@','track') or line.chomp.empty?

        begin
          yield parse(line)
        rescue EntryFileError
          skipped += 1
        end
      end
     
      skipped
    end
    
    # Return all of the chromosomes available in this EntryFile
    # This has fairly bad performance (requires parsing all entries in the file)
    # So it can be overridden if there is a better way for a specific filetype
    def chromosomes
      # Cache for performance
      if @chromosomes.nil?
        s = Set.new
        self.each { |entry| s << entry.chr }
        @chromosomes = s.to_a
      end
      
      return @chromosomes
    end
    
    # Allow EntryFiles to be indexed like GenomicDatas, optionally with a block
    # to avoid loading all elements of a chromosome into memory
    def [](chr_id, &block)
      chr(chr_id, &block)
    end
    
    # Query for a specific chromosome
    def chr(chr_id)
      if block_given?
        self.each(chr_id) { |entry| yield entry }
      else
        # TODO: Return a custom Enumerable object that doesn't require
        # preloading all of the entries
        entries = Array.new
        self.each(chr_id) { |entry| entries << entry }
        return entries
      end
    end
    
    # Count the number of entries that will be returned
    # Fairly bad performance, so can be overridden
    def count(chr = nil, start = nil, stop = nil)
      num = 0
      self.each(chr, start, stop) { |entry| num += 1 }
      return num
    end
    
    def to_bed(output)
      to_disk(output) { |entry| entry.to_bed }
    end
    
    def to_bedgraph(output)
      to_disk(output) { |entry| entry.to_bedgraph }
    end
    
    def to_gff(output)
      to_disk(output) { |entry| entry.to_gff }
    end
    
    private
    
    # Returns true if the index file exists (i.e. assume it has been indexed)
    def indexed?
      @index_file and File.exist?(@index_file)
    end

    # Returns true if this file is currently being indexed by any process
    def indexing?
      File.exist?(@indexing)
    end

    # Flag this file as currently being indexed
    def indexing=(b)
      if b
        FileUtils.touch(@indexing)
      else
        File.delete(@indexing) if File.exist?(@indexing)
      end
    end
    
    # Should be overridden in subclasses to parse an line into an object
    def parse(line)
      raise EntryFileError, "Do not know how to parse the entries in #{File.basename(@data_file)}!"
    end
    
    # Should be overridden in subclasses to query the entry file for lines
    def query_lines(chr = nil, start = nil, stop = nil)
      raise EntryFileError, "Do not know how to query for lines in #{File.basename(@data_file)}!"
    end
    
    def to_disk(output)
      File.open(File.expand_path(output), 'w') do |f|
        self.each do |entry|
          f.puts yield(entry)
        end
      end
    end
  end

  class TextEntryFile < EntryFile
    def initialize(filename, chr_col, start_col, end_col)
      super(filename)
      
      @chr_col = chr_col
      @start_col = start_col
      @end_col = end_col
      
      @filtered_file = @data_file + '.filtered'
      @sorted_file = @data_file + '.sorted'
      @bgzipped_file = @data_file + '.bgz'
      @index_file = @bgzipped_file + '.tbi'
    end
    
    # Delete the Tabix index file and the BGZipped version, if it exists
    def close
      File.delete(@index_file) if indexed?
      File.delete(@bgzipped_file) if bgzipped?
    end
    
    # Use wc to count the number of entries (assume one entry per line)
    def count(chr = nil, start = nil, stop = nil)  
      if chr and start.nil?
        %x[ grep -w #{chr} #{@data_file} | wc -l ].chomp.to_i
      else
        num = 0
        self.each(chr, start, stop) { |entry| num += 1 }
        num
      end
    end
    
    private
    
    def bgzipped?
      File.exist?(@bgzipped_file)
    end
    
    def min_num_columns
      # Cache for performance
      @min_num_columns = [@chr_col, @start_col, @end_col].max if @min_num_columns.nil?
      return @min_num_columns
    end
    
    # By default, just parse intervals
    def parse(line)
      entry = line.chomp.split("\t")
      raise EntryFileError, "Invalid interval entry" if entry.length < min_num_columns
      Genomics::Interval.new(entry[@chr_col-1], entry[@start_col-1].to_i, entry[@stop_col-1].to_i)
    end
    
    # Get all lines in the file matching chr:start-stop
    def query_lines(chr = nil, start = nil, stop = nil)
      raise EntryFileError, "Tabix only supports queries with start AND stop" if start and stop.nil?
      
      # If we're getting all entries, just use File#foreach
      if chr.nil?
        File.foreach(@data_file) { |line| yield line }
      # If we're getting a specific chromosome, use grep to filter the entries
      elsif start.nil? or stop.nil?
        File.grep(@data_file, chr) { |line| yield line }
      # If we're querying for a specific region, use Tabix to index the file
      else
        index() if not indexed?
        low = [start, stop].min
        high = [start, stop].max
        Tabix.query(@bgzipped_file, chr, low, high) { |line| yield line }
      end
    end
    
    # Index all TextEntryFiles with Tabix
    def index
      # Block if this file is currently being indexed by another thread/process
      if indexing?
        while indexing?
          sleep 1
        end
      else
        self.indexing = true
        
        begin
          # Filter unparseable entries
          # TODO: Find a more efficient way to filter unparseable entries without
          # having to copy the entire file line by line
          filtered = 0
          File.open(@filtered_file, 'w') do |f|
            File.foreach(@data_file) do |line|
              begin
                parse(line)
                f.write line
              rescue
                filtered += 0
              end
            end
          end
          puts "Filtered #{filtered} unparseable entries" if filtered > 0 and ENV['DEBUG']
          
          # File must be sorted
          File.sort(@filtered_file, @sorted_file, "-k#{@chr_col},#{@chr_col} -k#{@start_col},#{@start_col}n")
          
          # and BGZipped
          BGZip.compress(@sorted_file, @bgzipped_file)
                    
          # Now Tabix can index it
          Tabix.index(@bgzipped_file, @chr_col, @start_col, @end_col)
        rescue
          raise EntryFileError, "Error indexing file #{File.basename(@data_file)} for lookup!"
        ensure
          # Delete the temporary filtered and sorted files since they are unneeded
          File.delete(@filtered_file) if File.exist?(@filtered_file)
          File.delete(@sorted_file) if File.exist?(@sorted_file)
          self.indexing = false
        end
      end
    end
  end

  class BinaryEntryFile < EntryFile
  end

  class EntryFileError < StandardError
  end
end
