require 'rspec'
require 'bio-genomic-file'
require 'simplecov'
include Bio

SimpleCov.start do
  add_filter "spec"
  add_group "Common", "common"
end
