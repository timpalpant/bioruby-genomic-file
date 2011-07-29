require 'utils/fixed_precision'
require 'sparse_array'

##
# Represents a contiguous block of genomic data values
# Returned when querying a WigFile / SpotArray
#
##
module Bio
  module Genomics
    class Contig < SparseArray
      attr_accessor :chr

      def initialize(chr = 'unknown', span = 1)
        super(span)
        @chr = chr
      end
            
      # Does this Contig have a fixed step size?
      # Return the value if it does, false otherwise
      def fixed_step?
        prev_step, prev = nil, nil
        self.indices.sort.each do |i|
          if prev
            curr_step = i - prev
            if prev_step
              return false if prev_step != curr_step
            end
            
            prev_step = curr_step
          end
          
          prev = i
        end
        
        return prev_step
      end

      ##
      # OUTPUT METHODS
      ##
      
      # Output this Contig as a variableStep Wiggle block
      def to_variable_step
        @data.map { |i,v| "#{i}\t#{v.to_s(5)}" }.join("\n")
      end
      
      # Output this Contig as a fixedStep Wiggle block
      # NOTE: This will be invalid if the data is not actually fixed-step (irregularly spaced)
      def to_fixed_step
        raise ContigError, "Cannot output variable-step Contig to fixed-step format!" if not fixed_step?
        values.map { |v| v.to_s(5) }.join("\n")
      end
    end
    
    # Raised if something goes wrong with a Contig
    class ContigError < StandardError
    end
  end
end


# For converting an Array to a Contig
class Array
  def to_contig(chr = 'unknown', start = 1, step = 1, span = 1)
    contig = Bio::Genomics::Contig.new(chr, span)
    
    bp = start
    self.each do |value|
      contig.set(bp, value)
      bp += step
    end
    
    return contig
  end
end
