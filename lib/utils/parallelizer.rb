#
#  parallelizer.rb
#  BioRuby
#  Parallelize enumerable computations with Parallel gem
#  See: https://github.com/grosser/parallel
#
#  Created by Timothy Palpant on 6/1/11.
#  Copyright 2011 UNC. All rights reserved.
#

require 'tmpdir'
require 'parallel'

module Enumerable  
  # Parallel each: iterate over each element using multiple parallel processes
  # NOTE: Each process has its own variable-space, so changes must be
  # persisted to disk
  def p_each(opts = {}, &block)  
    Parallel.each(self, opts, &block)
  end
  
  # Parallel each: iterate over each element using multiple parallel processes
  # NOTE: Each process has its own variable-space, so changes must be
  # persisted to disk
  def p_map(opts = {}, &block)  
    Parallel.map(self, opts, &block)
  end
end
