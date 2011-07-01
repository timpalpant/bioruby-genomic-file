require 'entry_file'
require 'spot_file'
require 'spot'

##
# An entry in a GFF file
# For the format spec, see: http://genome.ucsc.edu/FAQ/FAQformat.html#format3
##
class GFFEntry < Spot
  attr_accessor :source, :feature, :frame
    
  def self.parse(line)
    begin
      record = line.chomp.split("\t")
      raise GFFError, "Invalid GFF Entry: GFF must have 9 columns" if record.length < 9
      
      spot = self.new
      spot.chr = record[0]
      spot.source = record[1]
      spot.feature = record[2]
      spot.start = record[3].to_i
      spot.stop = record[4].to_i
      spot.value = record[5].to_f unless record[5] == '.'
      strand = record[6]
      spot.frame = record[7]
      spot.id = record[8].split(';').first[9..-1]
      
      # Ensure that the coordinates (start/stop) match the strand, if specified
      if strand == '+'
        tmp_low = spot.low
        tmp_high = spot.high
        spot.start = tmp_low
        spot.stop = tmp_high
      elsif strand == '-'
        tmp_low = spot.low
        tmp_high = spot.high
        spot.start = tmp_high
        spot.stop = tmp_low
      end
        
      return spot
    rescue
      raise GFFError, "Not a valid GFF Entry" 
    end
  end
  
  def seqname
    @chr
  end
  
  def group
    @id
  end
  
  def score
    @value
  end
  
  # Override #to_gff to use all fields
  def to_gff
    source = @source ? @source : 'SpotArray'
    feature = @feature ? @feature : 'feature'
    value = @value ? @value : '.'
    frame = @frame ? @frame : '.'
    id = @id ? @id : 'no_id'
    
    "#{@chr}\t#{source}\t#{feature}\t#{low}\t#{high}\t#{value}\t#{strand}\t#{frame}\tprobe_id=#{id};count=1"
  end
end


##
# Stream GFF files
##
class GFFFile < TextEntryFile
  include SpotFile
  
  CHR_COL = 1
  START_COL = 4
  END_COL = 5
  
  def initialize(filename)
    super(filename, CHR_COL, START_COL, END_COL)
  end

  private
  
  # Define how to parse GFF entries
  def parse(line)
    GFFEntry.parse(line)
  end
end


class GFFError < EntryFileError
end
