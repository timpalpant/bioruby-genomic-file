require 'fileutils'

# Encapsulate methods in the File class 
# that call native Unix command-line utilities
# such as wc, grep, head, tail...
class File
  # Get the number of lines in a file using wc -l
  def self.num_lines(filename)
    if File.directory?(filename)
      # Recurse directories
      total = 0
      Dir.glob(filename + '/*').each do |f2|
        total += File.num_lines(f2)
      end
      return total
    end
    
    # Add 1 because wc counts line breaks?
    %x[ wc -l #{filename} ].split(' ').first.to_i + 1
  end

  # Get the number of characters in a file using wc -c
  def self.num_chars(filename)
    %x[ wc -c #{filename} ].split(' ').first.to_i
  end

  # Get the number of words in a file using wc -c
  def self.num_words(filename)
    %x[ wc -w #{filename} ].split(' ').first.to_i
  end
  
  # Returns whether or not +file+ is a binary file.  Note that this is
  # not guaranteed to be 100% accurate.  It performs a "best guess" based
  # on a simple test of the first +File.blksize+ characters.
  #
  # Example:
  #
  #   File.binary?('somefile.exe') # => true
  #   File.binary?('somefile.txt') # => false
  #--
  # Based on code originally provided by Ryan Davis (which, in turn, is
  # based on Perl's -B switch).
  #
  # Adapted from ptools gem, and fixed to also count newline characters
  #
  def self.binary?(filename, threshold = 0.30)
    expanded = File.expand_path(filename)
    # Get the first block of data from the file and split to characters
    s = (File.read(expanded, File.stat(expanded).blksize) || "").split(//)
    # Is the proportion of non-ASCII characters greater than the threshold?
    ((s.size - s.grep(/[\s\w\b]/).size) / s.size.to_f) > threshold
  end
  
  # Return an array of strings resulting from the output of grep -v
  # Alternatively, get all lines and then use Array#select
  # Takes an optional block
  def self.inverse_grep(filename, search_str)
    if block_given?
      IO.popen("grep -v -w '#{search_str}' #{filename}") do |output|
        output.each { |line| yield line }
      end
    else
      return %x[ grep -v '#{search_str}' #{filename} ]
    end
  end
  
  # Return an array of strings resulting from the output of grep -n
  def self.grep_with_linenum(filename, search_str)
    if block_given?
      IO.popen("grep -n -w '#{search_str}' #{filename}") do |output|
        output.each do |line|
          entry = line.split(':')
          line_num = entry.first.to_i
          content = entry.last
          yield(line_num, content)
        end
      end
    else
      return %x[ grep -n -w '#{search_str}' #{filename} ].split("\n").map do |line|
        entry = line.split(':')
        line_num = entry.first.to_i
        content = entry.last
        [line_num, content]
      end
    end
  end
  
  # Return an array of strings resulting from the output of grep
  def self.grep(filename, search_str)
    if block_given?
      IO.popen("grep -w '#{search_str}' #{filename}") do |output|
        output.each { |line| yield line }
      end
    else
      return %x[ grep -w '#{search_str}' #{filename} ]
    end
  end
  
  # Get lines m..n of a file
  def self.lines(filename, start_line, end_line)
    # If head and tail are available, use them (better performance)
    if has_head? and has_tail?
      num_lines = end_line - start_line + 1
      if block_given?
        IO.popen("tail -n +#{start_line} #{filename} 2>&1 | head -n #{num_lines}") do |pipe|
          pipe.each { |line| yield line }
        end
      else
        return %x[ tail -n +#{start_line} #{filename} 2>&1 | head -n #{num_lines} ]
      end
    # Otherwise use native Ruby implementation
    else
      buffer = File.stat(filename).blksize
      count = 0
      lines = Array.new
      
      File.open(filename) do |f|
        while count < start_line and 
          chunk = f.read(buffer)
          raise "There are only #{count} lines in #{File.basename(filename)}! (cannot get #{start_line}..#{end_line})" if chunk.nil?
          lines_in_chunk = chunk.count("\n")
          count += lines_in_chunk
        end
        
        # Move back one chunk
        idx = [f.pos - buffer, 0].max
        f.seek(idx)
        # Move to the next line unless we're at the beginning of the file
        f.gets unless idx == 0
        line_num = count - lines_in_chunk
        
        while line_num < end_line
          line = f.gets
          line_num += 1
          raise "There are only #{line_num} lines lines in #{File.basename(filename)}! (cannot get #{start_line}..#{end_line})" if line.nil?
          next if line_num < start_line
          
          if block_given?
            yield line
          else
            lines << line.chomp
          end
        end
      end

      return lines unless lines.length == 0
    end
  end
  
  # Get the first n lines of a file
  def self.head(filename, num_lines)
    return if num_lines <= 0
    
    count = 0
    lines = Array.new
    File.foreach(filename) do |line|
      if block_given?
        yield line
      else
        lines << line
      end
      
      count += 1
      break if count == num_lines
    end
    
    return lines.join unless lines.length == 0
  end
    
  # Get the last n lines of a file
  def self.tail(filename, num_lines)
    return if num_lines <= 0
    
    # Use tail command if it's available (better performance)
    if has_tail?
      if block_given?
        IO.popen("tail -n #{num_lines} #{File.expand_path(filename)}") do |pipe|
          pipe.each { |line| yield line }
        end
      else
        return %x[ tail -n #{num_lines} #{File.expand_path(filename)} ]
      end
    # Otherwise use native code
    else
      chunks = Array.new
      
      File.open(filename) do |f|
        buffer = f.stat.blksize
        idx = [f.size - buffer, 0].max
        count = 0

        begin
          # Seek to the desired position and read a chunk
          f.seek(idx)
          chunk = f.read(buffer)
          
          # Count the number of lines in the chunk
          count += chunk.count("\n")
          
          # If we're returning an Array of lines,
          # store the chunks now so we don't have to re-read them
          chunks.unshift(chunk) unless block_given?
          
          # Move to the next chunk
          idx -= buffer
        end while count < num_lines and f.pos > 0
      
        if block_given?
          skip = [count - num_lines, 0].max
          i = 0
          f.seek(idx+buffer)
          f.each_line do |line|
            yield line unless i <= skip
            i += 1
          end
        else
          return chunks.join.lines.reverse_each.take(num_lines).reverse
        end
      end
    end
  end

  # Find the location of an executable in the $PATH
  def self.which(cmd)
    # Use which command if it's available
    if has_which?
      result = %x[ which #{cmd} ].chomp
      return nil if result.empty?
      return result
    # Otherwise, search through the $PATH manually
    else
      exts = ENV['PATHEXT'] ? ENV['PATHEXT'].split(';') : ['']
      ENV['PATH'].split(File::PATH_SEPARATOR).each do |path|
        exts.each do |ext|
          exe = "#{path}/#{cmd}#{ext}"
          return exe if File.executable?(exe)
        end
      end
      return nil
    end
  end

  # Add a newline to the end of a file only if there isn't one already
  def self.newlinify(input_file, output_file)
    FileUtils.copy(input_file, output_file)
    if not %x[ tail -n 1 #{output_file} ].end_with?("\n")
      File.open(output_file, 'a') { |f| f.write "\n" }
    end
  end

  # Add a newline to the file in place
  def self.newlinify!(file)
    if not %x[ tail -n 1 #{file} ].end_with?("\n")
      File.open(file, 'a') { |f| f.write "\n" }
    end
  end
  
  # Concatenate files
  def self.cat(input_files, output_file)
    raise "Less than 2 input files passed to cat!" if input_files.length < 2

    # Ensure that there are newline characters at the end of each file
    input_files.each { |f| File.newlinify!(f) }

    # Cat all of the files together
    %x[ cat #{input_files.join(' ')} > #{output_file} ]
  end
  
  # Sort files
  def self.sort(input_file, output_file, options)
    %x[ sort #{options} -o '#{output_file}' '#{input_file}' ]
  end
  
  # Diff two files
  def self.diff(file1, file2)
    %x[ diff #{File.expand_path(file1)} #{File.expand_path(file2)} ]
  end
  
  ##
  # HELPER METHODS
  ##
  
  private
  
  # Look for which executable in the $PATH
  def self.has_which?
    # Cache for performance
    if not defined? @@has_which
      @@has_which = false
      exts = ENV['PATHEXT'] ? ENV['PATHEXT'].split(';') : ['']
      ENV['PATH'].split(File::PATH_SEPARATOR).each do |path|
        exts.each do |ext|
          exe = "#{path}/which#{ext}"
          if File.executable?(exe)
            @@has_which = true
            return @@has_which
          end
        end
      end
    end
    
    return @@has_which
  end
  
  # Look for head executable in the $PATH
  def self.has_head?
    # Cache for performance
    @@has_head = !self.which('head').nil? if not defined? @@has_head
    return @@has_head
  end
  
  # Look for tail executable in the $PATH
  def self.has_tail?
    # Cache for performance
    @@has_tail = !self.which('tail').nil? if not defined? @@has_tail
    return @@has_tail
  end
end
