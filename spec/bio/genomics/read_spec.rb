#
#  read_spec.rb
#  bioruby-genomic-file
#
#  Created by Timothy Palpant on 6/25/11.
#  Copyright 2011 UNC. All rights reserved.
#

require 'spec_helper'
require 'bio/genomics/read'

describe Genomics::Read do
  before do
    @test = Genomics::Read.new('chrI', 1, 36)
  end
  
  it "should correctly output to SAM format"
  
  it "should correctly output to Bowtie hits format"
end
