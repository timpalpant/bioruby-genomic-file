# Encapsulates information about an individual nucleosome
module Bio
  module Genomics
    class Nucleosome < Spot
      attr_accessor :conditional_position, :dyad, :dyad_stdev, :dyad_mean

      def occupancy=(value)
        @value = value
      end
      
      def occupancy
        @value
      end
      
      # Use the dyad as the nucleosome position
      def position
        dyad
      end
      
      def to_s
        "#{@chr}\t#{@start}\t#{@stop}\t#{@dyad}\t#{@dyad_stdev}\t#{@conditional_position}\t#{@dyad_mean}\t#{occupancy}"
      end
    end
  end
end
