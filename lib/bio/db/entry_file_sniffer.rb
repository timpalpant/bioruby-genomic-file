#
#  entry_file_sniffer.rb
#  ruby-genomics
#
#  Attempt to autodetect entry file types
#
#  Created by Timothy Palpant on 6/23/11.
#  Copyright 2011 UNC. All rights reserved.
#

require 'bio/utils/unix_file_utils'

class EntryFileSniffer
  def initialize(filename)
    @data_file = File.expand_path(filename)
  end
  
  # Sniff a filename and return its type (Class)
  def self.sniff(filename)
    self.new(filename).type
  end
  
  # Attempt to auto-sniff the type of file
  # and return the appropriate Class
  # NOTE: Order is important because some types are
  # not distinguishable (e.g. Bed vs. BedGraph)
  def type
    if bed?
      BedFile
    elsif bedgraph?
      BedGraphFile
    elsif sam?
      SAMFile
    elsif bam?
      BAMFile
    elsif binary?
      BinaryEntryFile
    else
      TextEntryFile
    end
  end
  
  def binary?
    @binary = File.binary?(@data_file) if not defined? @binary
    return @binary
  end
  
  def ascii?
    not binary?
  end
  
  def bed?
    return false if binary?
    
    begin
      BedEntry.parse(first_line)
    rescue
      return false
    end
    
    return true
  end
  
  def bedgraph?
    return false if binary?
    
    begin
      BedGraphEntry.parse(first_line)
    rescue
      return false
    end
    
    return true
  end
  
  def sam?
    return false if binary?
    
    begin
      SAMEntry.parse(first_line)
    rescue
      return true
    end
    
    return false
  end
  
  def bam?
    return false unless binary?
    
    # TODO: Better checking for BAM files
    return true
  end
  
  private
  
  # Get the first line of the file for sniffing
  def first_line
    if not defined? @first_line
      File.open(@data_file) do |f|
        @first_line = f.gets.chomp
      end
    end
    
    return @first_line
  end
end

class EntryFileSnifferError < StandardError
end
