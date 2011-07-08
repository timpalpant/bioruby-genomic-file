#
#  Assembly.rb
#  BioRuby
#
#  Created by Timothy Palpant on 5/25/11.
#  Copyright 2011 UNC. All rights reserved.
#

require 'bio/genomics/data_set'

module Bio
  module Genomics
    class Assembly < DataSet  
      attr_reader :name, :len_file
      
      # Initialize a new assembly
      def initialize(name, len_file = nil)
        @name = name
        @len_file = len_file
      end
      
      # Write this assembly file to disk in len format (chr_id \t num_bases)
      def to_len(filename)
        File.open(File.expand_path(filename), 'w') do |f|
          self.each do |chr, chr_len|
            f.puts "#{chr}\t#{chr_len}"
          end
        end
      end
      
      # Summary info about this Assembly
      def to_s
        str = "Assembly: #{@name}"
        self.each do |chr, chr_length|
          str += "\n\t#{chr}: #{chr_length} bases"
        end
        
        return str
      end

      # Load an assembly specified by name (builtin) or chrom lengths file
      def self.load(len_file)
        assembly = self.new(name, File.expand_path(len_file))
        File.foreach(File.expand_path(len_file)) do |line|
          entry = line.chomp.split("\t")
          raise "Invalid entry in len file #{name}" if entry.length != 2
          assembly[entry[0]] = entry[1].to_i
        end

        return assembly
      end
      
    end
  end
end
