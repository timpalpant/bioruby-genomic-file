#
#  entry_file_sniffer_spec.rb
#  BioRuby
#
#  Created by Timothy Palpant on 4/8/11.
#  Copyright 2011 UNC. All rights reserved.
#

require 'spec_helper'
require 'bio/entry_file_sniffer'

describe EntryFileSniffer do
  BED = File.expand_path(File.dirname(__FILE__) + '/../fixtures/test.bed')
  BEDGRAPH = File.expand_path(File.dirname(__FILE__) + '/../fixtures/test.bedGraph')
  SAM = File.expand_path(File.dirname(__FILE__) + '/../fixtures/test.sam')
  BAM = File.expand_path(File.dirname(__FILE__) + '/../fixtures/test.bam')
  UNKNOWN = File.expand_path(File.dirname(__FILE__) + '/../fixtures/test.fa')
  
  it "should autodetect Bed files" do
    EntryFileSniffer.sniff(BED).should == BedFile
  end
  
  it "should autodetect BedGraph files" do
    EntryFileSniffer.sniff(BEDGRAPH).should == BedGraphFile
  end
  
  it "should autodetect SAM files" do
    EntryFileSniffer.sniff(SAM).should == SAMFile
  end
  
  it "should autodetect BAM files" do
    EntryFileSniffer.sniff(BAM).should == BAMFile
  end
  
  it "should raise an error if the filetype is unknown" do
    lambda { EntryFileSniffer.sniff(UNKNOWN) }.should raise_error
  end
end
