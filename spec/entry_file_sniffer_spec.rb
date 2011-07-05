#
#  entry_file_sniffer_spec.rb
#  BioRuby
#
#  Created by Timothy Palpant on 4/8/11.
#  Copyright 2011 UNC. All rights reserved.
#

require 'spec_helper'
require 'bio-genomic-file'

describe EntryFileSniffer do
  BED_FILE = File.expand_path(File.dirname(__FILE__) + '/fixtures/test.bed')
  BEDGRAPH_FILE = File.expand_path(File.dirname(__FILE__) + '/fixtures/test.bedGraph')
  SAM_FILE = File.expand_path(File.dirname(__FILE__) + '/fixtures/test.sam')
  BAM_FILE = File.expand_path(File.dirname(__FILE__) + '/fixtures/test.bam')
  UNKNOWN_FILE = File.expand_path(File.dirname(__FILE__) + '/fixtures/test.fa')
  
  it "should autodetect Bed files" do
    EntryFileSniffer.sniff(BED_FILE).should == BedFile
  end
  
  it "should autodetect BedGraph files" do
    EntryFileSniffer.sniff(BEDGRAPH_FILE).should == BedGraphFile
  end
  
  it "should autodetect SAM files" do
    EntryFileSniffer.sniff(SAM_FILE).should == SAMFile
  end
  
  it "should autodetect BAM files" do
    EntryFileSniffer.sniff(BAM_FILE).should == BAMFile
  end
  
  it "should raise an error if the filetype is unknown" do
    lambda { EntryFileSniffer.sniff(UNKNOWN_FILE) }.should raise_error
  end
end
