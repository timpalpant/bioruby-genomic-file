require 'spec_helper'
require 'bio/genomics/nucleosome'

describe Genomics::Nucleosome do
  before do
    @test = Genomics::Nucleosome.new('chrI', 100, 247)
    @test.dyad = 170
    @test.dyad_stdev = 23.4
    @test.conditional_position = 0.00012
    @test.dyad_mean = 175
    @test.occupancy = 520
  end
  
  it "should correctly output to NukeCalls format" do
    @test.to_s.should == "chrI\t100\t247\t170\t23.4\t0.00012\t175\t520"
  end
end
