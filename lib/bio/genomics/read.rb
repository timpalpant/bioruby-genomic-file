require 'bio'
require 'bio/genomics/interval'

##
# A genomic interval with an associated sequence (read)
# i.e. for high-throughput, short-read sequencing
##
module Bio
  module Genomics
    class Read < Interval
      attr_accessor :seq, :qual
      
      def to_s
        "Read: #{chr},#{@start},#{@stop},#{@seq}"
      end
    end
  end
end
