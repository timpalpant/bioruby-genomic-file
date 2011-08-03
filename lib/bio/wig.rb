require 'enumerator'
require 'stringio'
require 'utils/unix'
require 'bio/genomics/contig'
require 'bio/genomics/assembly'
require 'bio/utils/ucsc'
include Bio::Utils

module Bio
  # Base class for TextWigFile and BigWigFile
  class WigFile
    include Enumerable
    
    CHUNK_SIZE = 200_000
    
    attr_reader :track_header, :data_file
    
    ##
    # INSTANCE METHODS
    ##
    
    # Open a Wig file and parse its track/chromosome information
    def initialize(filename)
      @data_file = File.expand_path(filename)
      @track_header = UCSC::TrackHeader.new(:type => 'wiggle_0')
      @contigs_index = Array.new
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
    
    # By default, iterate over lines in the Wig file
    def each(&block)
      each_line(&block)
    end
    
    # Iterate over the 'lines' in this Wig file
    # Note that for BigWigFiles these do not correspond to actual lines
    # (since the file is binary)
    def each_line
      raise "Should be overridden in subclasses!"
    end
    
    # Iterate over each value in this WigFile
    def each_value
      info = nil
      self.each_line do |line|
        if line.chomp.empty? or line.start_with?('track')
          next
        elsif line.start_with?('fixedStep', 'variableStep')
          info = ContigInfo.parse(line)
        elsif line == 'n/a' or line == 'NaN'
          yield nil
        else
          value = if info.fixed_step?
            Float(line)
          else
            Float(line.split("\t").last)
          end
          
          info.span.times { yield value }
        end
      end
    end
    
    # Return an array of all chromosomes in this Wig file
    def chromosomes
      @contigs_index.collect { |contig_info| contig_info.chr }.uniq
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
    
    # Output a summary about this BigWigFile
    def summary
      str = StringIO.new
      @contigs_index.each { |contig_info| str << "chrom=#{contig_info.chr}, start=#{contig_info.start}, stop=#{contig_info.stop}\n" }
      str << "\nMean:\t#{mean}\n"
      str << "Standard deviation:\t#{stdev}"
      
      return str.string
    end
    
    def to_s
      str = StringIO.new
      str << "WigFile: #{File.basename(@data_file)}\n"
      @contigs_index.each { |contig_info| str << "#{contig_info}\n" }
      str << "Mean:\t#{mean}\n"
      str << "Standard deviation:\t#{stdev}"
      
      return str.string
    end
    
    def inspect
      to_s
    end
    
    # Number of values in the Wig file
    def num_bases
      # Cache for performance
      compute_stats if @num_bases.nil?
      return @num_bases
    end
    
    # The sum of all values
    def total
      # Cache for performance
      compute_stats if @total.nil?
      return @total
    end
    
    # The mean of all values
    def mean  
      # Cache for performance
      compute_stats if @mean.nil?
      return @mean unless num_bases == 0
    end
    
    # The standard deviation of all values
    def stdev
      compute_stats if @stdev.nil?
      return @stdev unless num_bases == 0
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
      result = @contigs_index.select { |contig_info| contig_info.chr == query_chr }
      raise WigError, "Wig does not include data for chromosome #{query_chr}" unless result.length > 0
      return result
    end
    
    # Compute the coverage, total, and stdev in a single iteration
    # since it is no more costly
    def compute_stats
      @num_bases = 0
      @total = 0
      sum_of_squares = 0.0
      
      self.each_value do |value|
        next if value.nil?
        @num_bases += 1
        @total += value
        sum_of_squares += value**2
      end
      
      @mean = @total.to_f / @num_bases
      variance = (sum_of_squares - @total*@mean) / @num_bases
      @stdev = Math.sqrt(variance)
    end
    
    ##
    # Holds info about a Contig in a WigFile
    ##
    class ContigInfo
      attr_accessor :type, :chr, :start, :stop, :step, :span, :header_line, :line_start, :line_stop
      
      FIXED_STEP = 'fixedStep'
      VARIABLE_STEP = 'variableStep'
      
      def initialize(type = FIXED_STEP, chr = 'unknown', start = 1, stop = 1, step = 1, span = 1, header_line = nil, line_start = nil, line_stop = nil)
        @type = type
        @chr = chr
        @start = start
        @stop = stop
        @step = step
        @span = span
        @header_line = header_line
        @line_start= line_start
        @line_stop = line_stop
      end
      
      def fixed_step?
        @type == FIXED_STEP
      end
      
      def variable_step?
        @type == VARIABLE_STEP
      end
      
      # Return which line contains data for a given base pair
      def line_for_bp(bp)
        raise WigError, "Cannot compute the line for a base pair with variableStep contigs" if @type == VARIABLE_STEP
        raise WigError, "Contig does not contain data for base pair #{bp}" if bp < start or bp > stop
        @line_start + (bp-@start)/@step
      end
      
      # Return the start base pair of a given line number
      def bp_for_line(line_num)
        raise WigError, "Cannot compute the base pair for a line with variableStep contigs" if @type == VARIABLE_STEP
        raise WigError, "Contig does not include line number #{line_num}" if line_num < @line_start or line_num > @line_stop
        @start + @step * (line_num-@line_start)
      end
      
      # Parse a fixedStep/variableStep line
      def self.parse(line)
        # Remove any leading/trailing whitespace
        line.strip!
        
        # Store the type of Contig (fixedStep / variableStep)
        info = self.new
        if line.start_with?(FIXED_STEP)
          info.type = FIXED_STEP
        elsif line.start_with?(VARIABLE_STEP)
          info.type = VARIABLE_STEP
        else
          raise WigError, "Not a valid fixedStep/variableStep line!"
        end
      
        # Parse the other tokens
        line.split(' ').each do |opt|
          keypair = opt.split('=')
          key = keypair.first
          value = keypair.last
          
          case key
            when 'chrom'
              info.chr = value
            when 'start'
              info.start = value.to_i
            when 'step'
              info.step = value.to_i
            when 'span'
              info.span = value.to_i
          end
        end

        return info
      end

      def to_s
        s = StringIO.new
        s << "#{@type} chrom=#{@chr}"
        s << " start=#{@start}" if @start and fixed_step?
        s << " step=#{@step}" if @step and fixed_step?
        s << " span=#{@span}" if @span
        return s.string
      end
    end

    class WigError < StandardError
    end
  end


  ##
  # For documentation, see: http://genome.ucsc.edu/goldenPath/help/bigWig.html
  # Analogous to WigFile, but for compressed BigWigs
  ##
  class BigWigFile < WigFile  
    attr_reader :min, :max

    def initialize(filename)
      super(filename)
      
      info = UCSC.bigwig_info(@data_file)
      raise WigError, "You must first convert your Wig file to BigWig" if info.length < 8
      
      info[7..-6].each do |line|
        entry = line.chomp.split(' ')
        
        contig_info = ContigInfo.new
        contig_info.chr = entry.first
        chr_length = entry.last.to_i
        contig_info.start = find_start_base(contig_info.chr, chr_length)
        contig_info.stop = find_stop_base(contig_info.chr, chr_length)
        
        @contigs_index << contig_info
      end
      
      # bigWigInfo doesn't calculate mean/stdev accurately enough when values are small (10^-7)
      #@mean = info[-4].chomp.split(':').last.to_f
      #@stdev = info[-1].chomp.split(':').last.to_f
      @min = info[-3].chomp.split(':').last.to_f
      @max = info[-2].chomp.split(':').last.to_f
    end
    
    ##
    # QUERY METHODS
    ##
    
    # Iterate over the BigWigFile as if it were a TextWigFile (by line)
    def each_line
      # Grab chunks at a time since it is more efficient than calling bigWigSummary for every single value
      @contigs_index.each do |contig|
        # Yield a header line (fixedStep) for each chromosome
        yield contig.to_s
        
        # Yield all of the values for that chromosome
        query_start = contig.start
        while query_start < contig.stop
          query_stop = [query_start+CHUNK_SIZE-1, contig.stop].min
          num_values = query_stop - query_start + 1
          begin
            UCSC.bigwig_summary(@data_file, contig.chr, query_start, query_stop, num_values, 'mean').each { |value| yield value }
          rescue UCSC::ToolsError
            puts "Interval #{contig.chr}:#{query_start}-#{query_stop} has no data" if ENV['DEBUG']
            num_values.times { yield 'n/a' }
          ensure
            query_start = query_stop + 1
          end
        end
      end
    end
    
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
      
      # Load the track information from the first line
      File.open(@data_file) do |f|
        begin
          @track_header = UCSC::TrackHeader.parse(f.gets)
        rescue UCSC::TrackHeaderError
          puts "Error parsing track header" if ENV['DEBUG']
        end
      end
      
      # Call grep to load the contig information
      puts 'Indexing Contig header lines' if ENV['DEBUG']
      File.grep_with_linenum(@data_file, 'chrom') do |line_num,line|      
        begin
          info = ContigInfo.parse(line)
          info.header_line = line_num
          line_start = line_num
          first_line = String.new
          while first_line.nil? or first_line.chomp.empty?
            line_start += 1
            first_line = File.lines(@data_file, line_start, line_start)
          end
          info.line_start = line_start
          @contigs_index << info
        rescue
          puts "Not a valid fixedStep/variableStep header" if ENV['DEBUG']
        end
      end
      
      # Now find the start and stop of each Contig
      puts 'Cataloging Contig start/stops' if ENV['DEBUG']
      @contigs_index.each do |contig_info|
        puts "\t#{contig_info}" if ENV['DEBUG']
        # fixedStep lines give the start, so we just need to find the stop
        if contig_info.fixed_step?
          # Get the line number of the next contig in the file
          next_contig_line = @contigs_index.select { |info| info.header_line > contig_info.header_line }.collect { |info| info.header_line }.sort.first
          
          # Get the line immediately before it (with data)
          offset = 0
          
          # If there is no contig after this one, then read to EOF
          if next_contig_line.nil?
            next_contig_line = File.num_lines(@data_file)
            offset -= 1
          end
          
          last_line = String.new
          while last_line.nil? or last_line.strip.empty?
            offset += 1
            last_line = File.lines(@data_file, next_contig_line-offset, next_contig_line-offset)
          end
          contig_info.line_stop = next_contig_line - offset
          
          # Calculate the stop based on the number of values and the step/span sizes
          num_values = contig_info.line_stop - contig_info.line_start
          contig_info.stop = contig_info.start + contig_info.step*(num_values) + contig_info.span - 1
        # for variableStep lines, we need to find the start and stop
        else
          # Find the start, i.e. the first base pair with data
          contig_info.start = File.lines(@data_file, contig_info.line_start, contig_info.line_start).split("\t").first.to_i
          
          # Find the stop, i.e. the last base pair with data
          # Get the line number of the next contig in the file
          next_contig_line = @contigs_index.select { |info| info.header_line > contig_info.header_line }.collect { |info| info.header_line }.sort.first
          
          # Get the line immediately before it (with data)
          offset = 0
          
          # If there is no contig after this one, then read to EOF
          if next_contig_line.nil?
            next_contig_line = File.num_lines(@data_file)
            offset -= 1
          end
          
          last_line = String.new
          while last_line.nil? or last_line.chomp.empty?
            offset += 1
            last_line = File.lines(@data_file, next_contig_line-offset, next_contig_line-offset)
          end
          contig_info.line_stop = next_contig_line - offset
          contig_info.stop = last_line.split("\t").first.to_i + contig_info.span - 1
        end
      end
    
      # Raise an error if no chromosomes were found
      raise WigError, "No fixedStep/variableStep headers found in Wig file!" if @contigs_index.length == 0
    end
    
    ##
    # QUERY METHODS
    ##
    
    # Iterate over each line in the WigFile
    def each_line
      File.foreach(@data_file) { |line| yield line }
    end
    
    # Override each_chunk to move through the file in a single pass for efficiency
    # rather than querying randomly throughout the file (which greatly increases seek time)
    def each_chunk
      
    end
    
    # Return a Contig of single-bp data from the specified region
    def query(chr, start = nil, stop = nil)
      start = chr_start(chr) if start.nil?
      stop = chr_stop(chr) if stop.nil?
      raise WigError, "Chromosome #{chr} not found in Wig file #{@data_file}!" unless include?(chr, start, stop)
      
      relev_contigs = @contigs_index.select { |info| info.chr == chr and (info.stop >= start and info.start <= stop) }.sort_by { |info| info.start }
      output = Genomics::Contig.new(chr)
      
      # Find the relevant contigs for the requested interval
      relev_contigs.each do |info|        
        if info.fixed_step?
          # Clamp to bases that are covered by this Contig
          low = [start, info.start].max
          high = [stop, info.stop].min
          
          # Figure out what lines in the file we need to get those bases
          start_line = info.line_for_bp(low)
          stop_line = info.line_for_bp(high)
          
          # Query the file for the lines and store them in the Contig
          bp = info.bp_for_line(start_line)
          File.lines(@data_file, start_line, stop_line) do |line|
            value = line.to_f
            (bp...bp+info.span).each { |base| output.set(base, value) if base >= start and base <= stop }
            bp += info.step
          end
        else
          # With variableStep, there's no way to know exactly what lines are needed
          # so iterate through them all until we find what we want
          File.lines(@data_file, info.line_start, info.line_stop) do |line|
            entry = line.split("\t")
            bp = entry.first.to_i
            # Skip until we've found the base we're interested in
            next if bp+info.span-1 < start
            # Quit if we've reached our last base
            break if bp > stop
            
            value = entry.last.to_f
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
  end
end
