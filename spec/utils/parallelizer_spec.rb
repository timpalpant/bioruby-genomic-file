#
#  parallelizer_spec.rb
#  bioruby-genomic-file
#
#  Created by Timothy Palpant on 6/28/11.
#  Copyright 2011 UNC. All rights reserved.
#

require 'spec_helper'
require 'utils/parallelizer'

describe Enumerable do
  before do
    @test = [1, 2, 3, 4, 5, 6, 7, 8]
  end
  
  context "#p_each" do
    it "should iterate over all entries once" do
      output_file = File.expand_path(File.dirname(__FILE__) + '/../fixtures/iteration-test.txt')
      begin
        @test.p_each do |entry| 
          File.open(output_file, 'a+') do |f| 
            f.puts entry
          end
        end
        
        result = File.readlines(output_file)
        result.length.should == 8
        result.each { |line| @test.delete(line.chomp.to_i) }
        @test.length.should == 0
      ensure
        File.delete(output_file) if File.exist?(output_file)
      end
    end
    
    it "should allow the number of threads to be set" do
      @test.p_each(:in_threads => 4) { |e| e }
    end
    
    it "should allow the number of processes to be set" do
      @test.p_each(:in_processes => 2) { |e| e }
    end
  end
  
  context "#p_map" do
    it "should return all entries" do
      result = @test.p_map { |e| e }
      result.length.should == @test.length
      result.each { |e| @test.delete(e) }
      @test.length.should == 0
    end
    
    it "should allow the number of threads to be set" do
      @test.p_map(:in_threads => 4) { |e| e }
    end
    
    it "should allow the number of processes to be set" do
      @test.p_map(:in_processes => 2) { |e| e }
    end
  end
end
