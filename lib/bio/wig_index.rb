module Bio
  class WigIndex < Array    
    INDEX_EXTENSION = '.widx'
    
    attr_accessor :num_bases, :total, :mean, :stdev, :min, :max
        
    # Load a Wig index from disk
    def self.load(filename)
      obj = nil
      File.open(File.expand_path(filename)) { |f| obj = Marshal.load(f) }
      return obj
    end
        
    # Save this Wig index to disk
    def to_disk(filename)
      File.open(File.expand_path(filename), 'w') do |f|
        Marshal.dump(self, f)
      end
    end
  end
  
  ##
  # Holds info about a Contig in a WigFile
  ##
  class ContigInfo
    attr_accessor :chr, :start, :stop, :span
    
    FIXED_STEP = 'fixedStep'
    VARIABLE_STEP = 'variableStep'
    BIGWIG = 'bigWig'
    
    def initialize(type, chr = 'unknown', start = 1, stop = 1, span = 1)
      @type = type
      @chr = chr
      @start = start
      @stop = stop
      @span = span
    end
    
    def fixed_step?
      @type == FIXED_STEP
    end
    
    def variable_step?
      @type == VARIABLE_STEP
    end
    
    def bigwig?
      @type == BIGWIG
    end
  end
  
  class BigWigContigInfo < ContigInfo
    def initialize(chr = 'unknown', start = 1, stop = 1)
      super(BIGWIG, chr, start, stop)
      @start = start
      @stop = stop
    end
    
    def to_s
      "#{@type} chrom=#{@chr} start=#{@start} span=#{@span}"
    end
  end
  
  class TextWigContigInfo < ContigInfo
    attr_accessor :line_start, :line_stop, :pos_start, :pos_stop, :index
        
    def initialize(type = FIXED_STEP, chr = 'unknown', start = 1, stop = 1, span = 1)
      super(type, chr, start, stop, span)
      @index = Hash.new
    end
    
    # Parse a fixedStep/variableStep line
    def self.parse(line)
      # Remove any leading/trailing whitespace
      line.strip!
      
      # Store the type of Contig (fixedStep / variableStep)
      info = if line.start_with?(FIXED_STEP)
        FixedStepContigInfo.new
      elsif line.start_with?(VARIABLE_STEP)
        VariableStepContigInfo.new
      else
        raise ContigError, "Not a valid fixedStep/variableStep line!"
      end
    
      # Parse the other tokens
      line.split(' ').each do |opt|
        keypair = opt.split('=')
        key = keypair.first
        value = keypair.last
        
        case key
          when 'chrom'
            info.chr = value
          when 'start'
            info.start = value.to_i
          when 'step'
            info.step = value.to_i
          when 'span'
            info.span = value.to_i
        end
      end

      return info
    end
    
    def store_index(bp, pos)
      @index[bp] = pos
    end
    
    def get_index(bp)
      @index[bp]
    end
    
    def upstream_indexed_bp(b)
      @index.select { |bp,pos| bp <= b }.keys.max
    end
  end
  
  class FixedStepContigInfo < TextWigContigInfo
    attr_accessor :step
    
    def initialize(chr = 'unknown', start = 1, stop = 1, span = 1, step = 1)
      super(FIXED_STEP, chr, start, stop, span)
      @step = step
    end
    
    # Return which line contains data for a given base pair
    def line_for_bp(bp)
      raise WigError, "Contig does not contain data for base pair #{bp}" if bp < start or bp > stop
      @line_start + (bp-@start)/@step
    end
    
    # Return the start base pair of a given line number
    def bp_for_line(line_num)
      raise WigError, "Contig does not include line number #{line_num}" if line_num < @line_start or line_num > @line_stop
      @start + @step * (line_num-@line_start)
    end
    
    def to_s
      "#{@type} chrom=#{@chr} start=#{@start} step=#{@step} span=#{@span}"
    end
  end
  
  class VariableStepContigInfo < TextWigContigInfo
    def initialize(chr = 'unknown', start = 1, stop = 1, span = 1)
      super(VARIABLE_STEP, chr, start, stop, span)
    end
    
    def to_s
      "#{@type} chrom=#{@chr} span=#{@span}"
    end
  end
  
  class ContigError < StandardError
  end
end