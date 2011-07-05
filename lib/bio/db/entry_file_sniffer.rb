#
#  entry_file_sniffer.rb
#  ruby-genomics
#
#  Attempt to autodetect entry file types
#
#  Created by Timothy Palpant on 6/23/11.
#  Copyright 2011 UNC. All rights reserved.
#

require 'ptools'

class EntryFileSniffer
  def initialize(filename)
    @data_file = File.expand_path(filename)
  end
  
  def binary?
    File.binary?(@data_file)
  end
  
  def ascii?
    not binary?
  end
  
  def bed?
  end
  
  def bedgraph?
  end
  
  # Sniff a filename and return its type (class)
  # if possible
  def self.sniff(filename)
    sniffer = self.new(filename)
    raise EntryFileSnifferError, "Could not autodetect file type!"
  end
end

class EntryFileSnifferError < StandardError
end