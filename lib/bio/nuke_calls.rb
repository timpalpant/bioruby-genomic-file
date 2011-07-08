require 'bio/entry_file'
require 'bio/spot_file'
require 'bio/genomics/nucleosome'

module Bio
  class NukeCallsFile < TextEntryFile
    include SpotFile
    
    HEADER = "#Chromosome\tNuke Start\tNuke Stop\tDyad\tDyad StDev\tSmoothed Position\tDyad Mean\tDyad Count"
   
    CHR_COL = 1
    # Use the dyad as the position for lookups with Tabix
    START_COL = 4
    END_COL = 4
    
    def initialize(filename)
      super(filename, CHR_COL, START_COL, END_COL)
    end
    
    # Can't query nucleosomes by id because they don't have any
    def id(query_id)
      raise NukeCallsError, "Cannot query nucleosome calls file by id!"
    end
    
    # Don't allow querying for values since it's not exactly a genomic data set
    def query(chr = nil, start = nil, stop = nil)
      raise NukeCallsError, "Cannot query nucleosome calls file for genomic values!"
    end

    private
    
    def parse(line)
      NukeCallsEntry.parse(line)
    end
    
    class NukeCallsEntry < Genomics::Nucleosome
      # Parse a NukeCallsFile entry
      def self.parse(line)
        begin
          entry = line.chomp.split("\t")
          raise NukeCallsError, "Invalid Nucleosome call: Nucleosome calls file must have 8 columns" if entry.length != 8
      
          nuke = self.new
          nuke.chr = entry[0]
          nuke.start = entry[1].to_i
          nuke.stop = entry[2].to_i
          nuke.dyad = entry[3].to_i
          nuke.dyad_stdev = entry[4].to_f
          nuke.conditional_position = entry[5].to_f
          nuke.dyad_mean = entry[6].to_i
          nuke.value = entry[7].to_f
          
          return nuke
        rescue
          raise NukeCallsError, "Invalid nucleosome call entry!"
        end
      end
    end
  end

  class NukeCallsError < EntryFileError
  end
end
