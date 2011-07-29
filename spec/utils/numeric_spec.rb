require 'spec_helper'
require 'utils/numeric'

describe String do
  before do
    @int = '123'
    @float = '-123.567'
    @alpha = 'hello'
    @nan = '4ha444.6'
  end
  
  context "when checking if it is numeric" do
    it "should detect integers" do
      @int.should be_integer
      @float.should_not be_integer
    end
  
    it "should detect floats" do
      @float.should be_numeric
      @int.should be_numeric
    end
  
    it "should detect words" do
      @alpha.should_not be_integer
      @alpha.should_not be_numeric
    end
  
    it "should detect non-numbers" do
      @nan.should_not be_integer
      @nan.should_not be_numeric
    end
  end
end