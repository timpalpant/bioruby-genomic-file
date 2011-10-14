require 'stringio'
require 'utils/unix'
require 'bio/genomics/contig'
require 'bio/genomics/assembly'
require 'bio/utils/ucsc'
require 'bio/wig_index'
include Bio::Utils

module Bio
  # Base class for TextWigFile and BigWigFile
  class WigFile
    KEY_GRANULARITY = 10_000
    CHUNK_SIZE = 200_000
    
    attr_reader :track_header, :data_file
    
    ##
    # INSTANCE METHODS
    ##
    
    # Open a Wig file and parse its track/chromosome information
    def initialize(filename)
      @data_file = File.expand_path(filename)
      @track_header = UCSC::TrackHeader.new(:type => 'wiggle_0')
      @index = WigIndex.new
    end
    
    # Put any cleanup operations here
    def close
    end
    
    # Open a Wig file with an optional block
    def self.open(filename)
      wig = self.new(filename)
    
      if block_given?
        yield wig
      else
        return wig
      end
      
      wig.close
    end
    
    # Autodetect whether this is a Wig or BigWig, and return the correct one
    def self.autodetect(filename, &block)
      if File.binary?(filename)
        BigWigFile.open(filename, &block)
      else
        TextWigFile.open(filename, &block)
      end
    end
    
    # Return an array of all chromosomes in this Wig file
    def chromosomes
      @index.collect { |contig_info| contig_info.chr }.uniq
    end
    
    # Does this Wig file include data for a given locus?
    def include?(chr, start = nil, stop = nil)
      if start.nil?
        chromosomes.include?(chr)
      elsif stop.nil?
        chromosomes.include?(chr) and start >= chr_start(chr) and start <= chr_stop(chr)
      else
        chromosomes.include?(chr) and start >= chr_start(chr) and stop <= chr_stop(chr)
      end
    end
    
    # Get the lowest base pair of data for a chromosome
    def chr_start(query_chr)
      chr_contigs(query_chr).collect { |contig_info| contig_info.start }.min
    end
    
    # Get the highest base pair of data for a chromosome
    def chr_stop(query_chr)
      chr_contigs(query_chr).collect { |contig_info| contig_info.stop }.max
    end
    
    def num_bases
      @index.num_bases
    end
    
    def total
      @index.total
    end
    
    def mean
      @index.mean
    end
    
    def stdev
      @index.stdev
    end
    
    def min
      @index.min
    end
    
    def max
      @index.max
    end
    
    def to_s
      str = StringIO.new
      str << "WigFile: #{File.basename(@data_file)}\n"
      @index.each { |contig_info| str << "#{contig_info}\n" }
      str << "Base pairs of data:\t#{num_bases}\n"
      str << "Mean:\t#{mean}\n"
      str << "Standard deviation:\t#{stdev}"
      
      return str.string
    end
    
    def inspect
      to_s
    end
    
    ##
    # ABSTRACT METHODS
    ##
    
    # Return a Contig of data from the specified region
    def query(chr, start, stop)
      raise WigError, "Should be overridden in a base class (BigWigFile/WigFile)!"
    end
    
    ##
    # HELPER METHODS
    ##
    private
    
    # Get the contigs for a chromosome
    def chr_contigs(query_chr)
      result = @index.select { |contig_info| contig_info.chr == query_chr }
      raise WigError, "Wig does not include data for chromosome #{query_chr}" unless result.length > 0
      return result
    end

    class WigError < StandardError
    end
  end


  ##
  # For documentation, see: http://genome.ucsc.edu/goldenPath/help/bigWig.html
  # Analogous to WigFile, but for compressed BigWigs
  ##
  class BigWigFile < WigFile  
    def initialize(filename)
      super(filename)
      
      info = UCSC.bigwig_info(@data_file)
      raise WigError, "You must first convert your Wig file to BigWig" if info.length < 8
      
      info[7..-6].each do |line|
        entry = line.chomp.split(' ')
        
        chr = entry.first
        chr_length = entry.last.to_i
        start = find_start_base(chr, chr_length)
        stop = find_stop_base(chr, chr_length)

        @index << BigWigContigInfo.new(chr, start, stop)
      end
      
      # TODO: bigWigInfo doesn't calculate mean/stdev accurately enough when values are small (10^-7)
      @index.num_bases = info[-5].chomp.split(':').last.delete(',').to_i
      @index.mean = info[-4].chomp.split(':').last.to_f
      @index.total = @index.mean * @index.num_bases
      @index.stdev = info[-1].chomp.split(':').last.to_f
      @index.min = info[-3].chomp.split(':').last.to_f
      @index.max = info[-2].chomp.split(':').last.to_f
    end
    
    ##
    # QUERY METHODS
    ##
    
    # Return a Contig of data from the specified region
    def query(chr, start = nil, stop = nil, type = 'mean')
      # Don't query off the ends of chromosomes
      start = chr_start(chr) if start.nil?
      stop = chr_stop(chr) if stop.nil?
      raise WigError, "BigWig does not contain data for the interval #{chr}:#{start}-#{stop}" if not include?(chr, start, stop)

      
      # bigWigSummary segfaults if query is too big, so use nice, bite-sized chunks
      query_start = start
      contig = Genomics::Contig.new(chr)
      while query_start <= stop
        query_stop = [query_start+CHUNK_SIZE-1, stop].min
        num_values = query_stop - query_start + 1
        begin
          chunk = UCSC.bigwig_summary(@data_file, chr, query_start, query_stop, num_values, type)
        rescue UCSC::ToolsError
          raise WigError, "BigWig does not contain data for the interval #{chr}:#{start}-#{stop}"
        end
        
        # Store the chunk of values in the Contig
        chunk.each_with_index do |value,i|
          next if value == 'n/a' or value == 'NaN'
          contig.set(query_start+i, value.to_f)
        end
        
        query_start = query_stop + 1
      end
      
      return contig
    end

    # Return the average value for the specified region
    def query_average(chr, start, stop)
      # Don't query off the ends of chromosomes
      raise WigError, "BigWig does not contain data for the interval #{chr}:#{start}-#{stop}" if not include?(chr, start, stop)
      UCSC.bigwig_summary(@data_file, chr, start, stop, 1).first.to_f
    end
    
    ##
    # OUTPUT METHODS
    ##
    
    def to_wig(output_file)
      BigWig.to_wig(@data_file, output_file)
    end
    
    def to_bedgraph(output_file)
      BigWig.to_bedgraph(@data_file, output_file)
    end
    
    # Convert a BigWigFile to a Wig file
    def self.to_wig(input_file, output_file)
      puts "Converting BigWig file (#{File.basename(input_file)}) to Wig (#{File.basename(output_file)})" if ENV['DEBUG']
      
      header_file = File.expand_path(output_file + '.header')
      data_file = File.expand_path(output_file + '.data')
      
      begin
        # Write a track header
        File.open(header_file, 'w') do |f|
          f.puts UCSC::TrackHeader.new(:name => File.basename(output_file))
        end
        
        # Extract the data with UCSC tools
        UCSC.bigwig_to_wig(input_file, data_file)
        
        # Cat the two parts together
        File.cat([header_file, data_file], output_file)
      rescue
        raise WigError, "Error converting BigWig file to Wig"
      ensure
        # Delete the two temp files
        File.delete(header_file) if File.exist?(header_file)
        File.delete(data_file) if File.exist?(data_file)
      end
    end

    # Convert a BigWig to a BedGraph
    def self.to_bedgraph(input_file, output_file)
      puts "Converting BigWig file (#{File.basename(input_file)}) to BedGraph (#{File.basename(output_file)})" if ENV['DEBUG']
      UCSC.bigwig_to_bedgraph(input_file, output_file)
    end

    ##
    # HELPER METHODS
    ##

    private

    # Find the first base pair with data
    def find_start_base(chr, chr_length)
      # Start at the first base pair and look forwards
      bp = 1
      while bp <= chr_length
        start = bp
        stop = [start+CHUNK_SIZE-1, chr_length].min
        num_values = stop - start + 1
        begin
          result = UCSC.bigwig_summary(@data_file, chr, start, stop, num_values, 'coverage')
          return start + result.find_index('1')
        rescue UCSC::ToolsError
        end
        
        bp += CHUNK_SIZE
      end
      
      raise WigError, "Could not find start base pair in BigWig file #{File.basename(@data_file)} for chromosome #{chr}"
    end

    # Find the last base pair with data
    def find_stop_base(chr, chr_length)
      # Start at the last base pair and look backwards
      bp = chr_length
      while bp >= 1
        start = [1, bp-CHUNK_SIZE+1].max
        stop = bp
        num_values = stop - start + 1
        begin
          result = UCSC.bigwig_summary(@data_file, chr, start, stop, num_values, 'coverage')
          return stop - result.reverse.find_index('1')
        rescue UCSC::ToolsError
        end
        
        bp -= CHUNK_SIZE
      end
      
      raise WigError, "Could not find stop base pair in BigWig file #{File.basename(@data_file)} for chromosome #{chr}"
    end
  end


  ##
  # An ASCII text Wiggle file
  ##
  class TextWigFile < WigFile
    # Open a Wig file and parse its track/contig information
    def initialize(filename)
      super(filename)
      @f = File.open(@data_file)
      
      # Load the track information from the first line
      begin
        @track_header = UCSC::TrackHeader.parse(@f.gets)
      rescue UCSC::TrackHeaderError
        puts "Error parsing track header" if ENV['DEBUG']
      end
      
      # Index the contigs in this ASCII Wig file
      # Attempt to load the index from disk if it has previously been indexed and saved
      @index_file = @data_file+WigIndex::INDEX_EXTENSION
      is_indexed = false
      if File.exist?(@index_file)
        puts "Attempting to load and match index from file" if ENV['DEBUG']
        begin
          @index = WigIndex.load(@index_file)
          is_indexed = @index.matches?(@data_file)
        rescue
          puts "Error loading/matching index from file!"
        end
      end
      
      if not is_indexed
        index_contigs()
        # Save the index to disk if the KEEP_INDEX environment variable is set
        save_index if ENV['KEEP_INDEXES']
      end
    
      # Raise an error if no chromosomes were found
      raise WigError, "No fixedStep/variableStep headers found in Wig file!" if @index.length == 0
    end
    
    # Cleanup operations
    def close
      @f.close
    end
    
    # Manually force a save of the index
    def save_index
      @index.compute_digest(@data_file)
      @index.to_disk(@index_file)
    end
    
    ##
    # QUERY METHODS
    ##
    
    # Return a Contig of single-bp data from the specified region
    def query(chr, start = nil, stop = nil)
      start = chr_start(chr) if start.nil?
      stop = chr_stop(chr) if stop.nil?
      raise WigError, "Interval not found in Wig file #{File.basename(@data_file)}!" unless include?(chr, start, stop)
      
      relev_contigs = @index.select { |info| info.chr == chr and (info.stop >= start and info.start <= stop) }
      output = Genomics::Contig.new(chr)
      
      # Find the relevant contigs for the requested interval
      relev_contigs.each do |info|
        puts "Loading data from contig: #{info}" if ENV['DEBUG']
        # Clamp to bases that are covered by this Contig
        low = [start, info.start].max
        high = [stop, info.stop].min
        
        # Find the closest known upstream base-pair position in the index
        closest_indexed_bp = info.upstream_indexed_bp(low)
        @f.seek(info.get_index(closest_indexed_bp))
        puts "Found closest indexed bp #{closest_indexed_bp}" if ENV['DEBUG']
        
        if info.fixed_step?
          puts "Loading fixedStep data" if ENV['DEBUG']
          # Figure out what lines in the file we need to get those bases
          start_line = info.line_for_bp(low)
          stop_line = info.line_for_bp(high)
          puts "Need lines #{start_line}-#{stop_line}" if ENV['DEBUG']
          
          # Query the file for the lines and store them in the Contig
          current_line = info.line_for_bp(closest_indexed_bp)
          puts "At line #{current_line}, moving #{start_line-current_line} lines forward" if ENV['DEBUG']
          (start_line - current_line).times { @f.gets }
          current_line = start_line
          bp = info.bp_for_line(current_line)
          puts "Shifted to base pair: #{bp} on line #{current_line}" if ENV['DEBUG']
          while current_line <= stop_line
            line = @f.gets
            begin
              value = Float(line)
            rescue
              puts "Error parsing Float value: #{line}"
              value = line.to_f
            end
            (bp...bp+info.span).each { |base| output.set(base, value) if base >= start and base <= stop }
            bp += info.step
            current_line += 1
          end
        else
          puts "Loading variableStep data" if ENV['DEBUG']
          while (line = @f.gets)
            # Trim trailing newline character
            line.chomp!
            # Break if at the end of a Contig
            break if line.empty? or line.start_with?('fixedStep', 'variableStep')
            
            entry = line.split("\t")
            if entry.length != 2
              puts "Invalid variableStep line: #{line}" if ENV['DEBUG']
              next
            end
            begin
              bp = Integer(entry[0])
            rescue
              puts "Error parsing Integer: #{entry[0]}"
              bp = entry[0].to_i
            end
            # Skip until we've found the base we're interested in
            next if bp+info.span-1 < start
            # Quit if we've gone past our last base
            break if bp > stop

            begin
              value = Float(entry[1])
            rescue
              puts "Error parsing Float: #{entry[1]}"
              value = entry[1].to_f
            end
            (bp...bp+info.span).each { |base| output.set(base, value) if base >= start and base <= stop }
          end
        end
      end
      
      return output
    end
    
    ##
    # OUTPUT METHODS
    ##
    
    # Convert this TextWigFile to a BigWigFile
    def to_bigwig(output_file, assembly)
      TextWigFile.to_bigwig(@datafile, output_file, assembly)
    end

    # Convert this WigFile to a BedGraph
    def to_bedgraph(output_file)
      TextWigFile.to_bedgraph(@datafile, output_file)
    end
    
    # For converting wigs to BigWigs without having to load (index them) first
    def self.to_bigwig(input_file, output_file, assembly)
      puts "Converting Wig file (#{File.basename(input_file)}) to BigWig (#{File.basename(output_file)})" if ENV['DEBUG']
      UCSC.wig_to_bigwig(input_file, assembly.len_file, output_file)
    end
    
    # For converting wigs to BedGraph without having to load (index them) first
    # Also creates the most compact BedGraph possible by joining equal neighbors
    def self.to_bedgraph(input_file, output_file)
      puts "Converting Wig file #{File.basename(input_file)} to BedGraph #{File.basename(output_file)}" if ENV['DEBUG']
      raise "Not yet implemented"
    end
    
    ##
    # HELPER METHODS
    ##
    
    private
    
    def index_contigs
      puts 'Indexing Contig header lines' if ENV['DEBUG']
      
      num_bases = 0
      total = 0.0
      sum_of_squares = 0.0
      min, max = nil, nil
      
      @index = WigIndex.new
      @f.rewind
      line_num = 0
      num_contigs = 0
      bp = 0
      info = TextWigContigInfo.new
      until @f.eof?
        if num_contigs > 0 and (line_num+1 - info.line_start) % KEY_GRANULARITY == 0
          cursor = @f.pos
        end
        line = @f.gets.chomp
        line_num += 1
        
        # Track lines / blank lines
        if line.empty? or line.start_with?('track')
          next
        # Contig header lines
        elsif line.start_with?('fixedStep', 'variableStep')
          # If this signals the end of the previous Contig, store it
          if num_contigs > 0
            puts "...storing contig in the global index" if ENV['DEBUG']
            info.line_stop = line_num - 1
            info.pos_stop = @f.pos
            info.stop = bp + info.span - 1
            @index << info
          end
          
          info = TextWigContigInfo.parse(line)
          num_contigs += 1
          puts info.to_s if ENV['DEBUG']
          
          # Read until we find the first data value (there could, but shouldn't, be intervening whitespace)
          puts "...finding first data value" if ENV['DEBUG']
          line.clear
          while line.empty?
            cursor = @f.pos
            line = @f.gets.chomp
            line_num += 1
          end
          info.line_start = line_num
          info.pos_start = cursor
          # Move back one line so the first data line will be processed
          # on the next iteration
          @f.seek(cursor)
          line_num -= 1
          if info.variable_step?
            info.start = line.split("\t").first.to_i
          else
            bp = info.start - info.step
          end
          
          puts "...indexing and computing descriptive statistics" if ENV['DEBUG']
        # Data lines
        else
          # Compute stats since we're going through all of the data anyway
          if info.fixed_step?
            bp += info.step
            begin
              value = Float(line)
            rescue
              puts "Error parsing Float: #{line}"
              value = line.to_f
            end
          else
            entry = line.split("\t")
            begin
              bp = Integer(entry[0])
            rescue
              puts "Error parsing Integer: #{entry[0]}"
              bp = entry[0].to_i
            end
            
            begin
              value = Float(entry[1])
            rescue
              puts "Error parsing Float: #{entry[1]}"
              value = entry[1].to_f
            end
          end
          
          min = value if min.nil? or value < min
          max = value if max.nil? or value > max
          num_bases += info.span
          total += info.span * value
          sum_of_squares += info.span * value**2
          
          # Store the position of this line
          if (line_num-info.line_start) % KEY_GRANULARITY == 0
            info.store_index(bp, cursor)
          end
        end
      end
      
      # Store the last contig
      puts "...storing contig in the global index" if ENV['DEBUG']
      info.line_stop = line_num
      info.pos_stop = @f.pos
      info.stop = bp + info.span - 1
      @index << info
      
      @index.min = min
      @index.max = max
      @index.num_bases = num_bases
      @index.total = total
      @index.mean = total / num_bases
      variance = (sum_of_squares - total*@index.mean) / num_bases
      @index.stdev = Math.sqrt(variance)
    end
  end
end
