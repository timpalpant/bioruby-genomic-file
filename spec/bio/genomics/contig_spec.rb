require 'spec_helper'
require 'bio/genomics/contig'

describe Genomics::Contig do
  context "variable step" do
    context "with default span (= 1)" do
      before do
        @test = Genomics::Contig.new('chrI')
        @test.set(100, 1.5)
        @test.set(101, 1.5)
        (105..111).each { |bp| @test.set(bp, 2.1234321) }
      end
    
      it "should output to variableStep format" do
        @test.to_variable_step.should == "100\t1.5\n101\t1.5\n105\t2.1234\n106\t2.1234\n107\t2.1234\n108\t2.1234\n109\t2.1234\n110\t2.1234\n111\t2.1234"
      end
      
      it "should raise an error if attempting to output to fixedStep format" do
        lambda { @test.to_fixed_step }.should raise_error
      end
      
      it "should not be fixed step" do
        @test.should_not be_fixed_step
      end
    end
  
    context "with span = 5" do
      before do
        @test = Genomics::Contig.new('chrII', 5)
        @test.set(100, 1.5)
        @test.set(105, 2.0)
        @test.set(114, 3.5)
      end
  
      it "should output to variableStep format" do
        @test.to_variable_step.should == "100\t1.5\n105\t2\n114\t3.5"
      end
      
      it "should raise an error if attempting to output to fixedStep format" do
        lambda { @test.to_fixed_step }.should raise_error
      end
      
      it "should not be fixed step" do
        @test.should_not be_fixed_step
      end
    end
  end
  
  context "fixed step" do
    context "with default span (= 1)" do
      before do
        @test = Genomics::Contig.new('chrI')
        (105..111).each { |bp| @test.set(bp, 2.1234321) }
      end
    
      it "should output to variableStep format" do
        @test.to_variable_step.should == "105\t2.1234\n106\t2.1234\n107\t2.1234\n108\t2.1234\n109\t2.1234\n110\t2.1234\n111\t2.1234"
      end
      
      it "should output to fixedStep format" do
        @test.to_fixed_step.should == "2.1234\n2.1234\n2.1234\n2.1234\n2.1234\n2.1234\n2.1234"
      end
      
      it "should be fixed step" do
        @test.should be_fixed_step
      end
    end
  
    context "with span = 5" do
      before do
        @test = Genomics::Contig.new('chrII', 5)
        @test.set(100, 1.5)
        @test.set(110, 2.0)
        @test.set(120, 3.5)
      end
  
      it "should output to variableStep format" do
        @test.to_variable_step.should == "100\t1.5\n110\t2\n120\t3.5"
      end
      
      it "should output to fixedStep format" do
        @test.to_fixed_step.should == "1.5\n2\n3.5"
      end
      
      it "should be fixed step" do
        @test.should be_fixed_step
      end
    end
  end
end

describe Array do
  before do
    @test = [0,3,4,9,0,6,44,3,5,7,8,9,5,6,3,1,1,2,3,4,5,6,13,15,18,22]
  end
  
  it "should convert to Contig with default parameters: chr=unknown, start=1" do
    contig = @test.to_contig
    contig.chr.should == 'unknown'
    contig.start.should == 1
    contig.stop.should == @test.length
    contig.length.should == @test.length
    
    @test.each_with_index do |value,i|
      contig[i+1].should == @test[i]
    end
  end
  
  it "should allow parameters to be set" do
    contig = @test.to_contig('chrI', 25, 5, 5)
    contig.chr.should == 'chrI'
    contig.start.should == 25
    contig.stop.should == 25 + 5*(@test.length-1)+5-1
    contig.length.should == 5*(@test.length-1)+5
  end
end
