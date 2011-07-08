require 'spec_helper'
require 'bio/genomics/data_set'

describe Genomics::DataSet do
  TEST_STR = 'my data'
  TEST_STR2 = 'my other data'
  
  before do
    @test = Genomics::DataSet.new
    @test['I'] = TEST_STR
    @test['III'] = TEST_STR2
  end
  
  it "should return keys for chromosomes" do
    @test.chromosomes.should == ['I', 'III']
  end

  it "should return data for the correct chromosome" do
    @test.chr('I').should == TEST_STR
    @test.chr('III').should == TEST_STR2
    @test.chr('V').should == nil
  end
end
