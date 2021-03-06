require 'utils/unix'
require 'bio/genomics/contig'
require 'bio/utils/ucsc'

module Bio
  # Additional methods mixed in for EntryFile types that also
  # have numeric data/values
  module SpotFile
    # Get a Spot with a specific id
    def id(query_id)
      entries = Array.new
      skipped = 0
      File.grep(@data_file, query_id) do |line|
        begin
          entry = parse(line)
          entries << entry if entry.id == query_id
        rescue
          skipped += 1
        end
      end
      
      raise EntryFileError, "No spot with id #{query_id} in file #{File.basename(@data_file)}" if entries.length == 0
      raise EntryFileError, "More than one spot with id #{query_id} in file #{File.basename(@data_file)}" if entries.length > 1
      puts "Skipped #{skipped} invalid spots with id #{query_id}" if ENV['DEBUG']
      
      return entries.first
    end
    
    ##
    # STATISTICAL METHODS
    ##
    
    # The number of spots with values
    def num_values
      # Cache for performance
      compute_stats if @num_values.nil?    
      return @num_values
    end
    
    # The sum of the values of all spots
    def total
      # Cache for performance
      compute_stats if @total.nil? 
      return @total
    end
    
    # The mean value of all spots
    def mean
      # Cache for performance
      compute_stats if @mean.nil?
      return @mean unless num_values == 0
    end
    
    # The standard deviation of all spots
    def stdev
      compute_stats if @stdev.nil?
      return @stdev unless num_values == 0
    end
    
    ##
    # QUERY METHODS
    ##
    
    # Return a Contig of values for the given window
    def query(chr, start, stop)
      raise EntryFileError, "Invalid query! (#{start} > #{stop})" if start > stop
      contig = Genomics::Contig.new(chr)
      
      self.each(chr, start, stop) do |spot|
        # Get the high and low spot coordinates, and clamp to the ends of the window
        low = [start, spot.low].max
        high = [spot.high, stop].min
      
        for bp in low..high
          contig.set(bp, spot.value) unless spot.value.nil?
        end
      end
      
      return contig
    end
    
    ##
    # OUTPUT METHODS
    ##
    
    # Write this array to Wig format
    # By first writing to bedGraph, then calling bedGraphToBigWig
    # then calling bigWigToWig
    def to_wig(filename, assembly)    
      tmp_bigwig = filename + '.bw'
      begin
        self.to_bigwig(tmp_bigwig, assembly)
        Utils::UCSC.bigwig_to_wig(tmp_bigwig, filename)
      rescue
        raise "Error converting to Wig!"
      ensure
        # Delete the temporary intermediate files
        File.delete(tmp_bigwig) if File.exist?(tmp_bigwig)
      end
    end
    
    # Write this array to BigWig format
    # By first writing to bedGraph, then calling bedGraphToBigWig
    def to_bigwig(output_file, assembly)
      begin
        tmp_bedgraph = File.expand_path(output_file + '.bedGraph')
        self.to_bedgraph(tmp_bedgraph)
      
        # bedGraph must be sorted to call bedGraphToBigWig
        tmp_sorted = tmp_bedgraph + '.sorted'
        File.sort(tmp_bedgraph, tmp_sorted, '-k1,1 -k2,2n')
        Utils::UCSC.bedgraph_to_bigwig(tmp_sorted, assembly.len_file, output_file)
      rescue
        raise "Error converting to BigWig!"
      ensure
        # Delete the temporary intermediate files
        File.delete(tmp_bedgraph) if File.exist?(tmp_bedgraph)
        File.delete(tmp_sorted) if File.exist?(tmp_sorted)
      end
    end
    
    ##
    # HELPER METHODS
    ##
    
    private
    
    # Cache the sum of all values and the number of values
    # in a single iteration
    # See: http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
    def compute_stats
      @num_values = 0
      @total = 0
      sum_of_squares = 0.0
      
      self.each do |entry|
        next if entry.value.nil?
        
        # Update the count and sum
        @num_values += 1
        @total += entry.value
        
        # Update the sum of squares total (for computing variance)
        sum_of_squares += entry.value**2
      end
      
      @mean = @total.to_f / @num_values
      variance = (sum_of_squares - @total*@mean) / @num_values
      @stdev = Math.sqrt(variance)
    end
  end
end
