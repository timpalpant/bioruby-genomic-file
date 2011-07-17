#
#  samtools.rb
#  ruby-genomics
#  Wrapper for the samtools executable
#
#  Created by Timothy Palpant on 5/30/11.
#  Copyright 2011 UNC. All rights reserved.
#

# For documentation, see: http://samtools.sourceforge.net/samtools.shtml

module Bio
  module Utils
    module SAMTools
      
      # Extract/print all or sub alignments in SAM or BAM format. 
      # If no region is specified, all the alignments will be printed; 
      # otherwise only alignments overlapping the specified regions will be output. 
      # An alignment may be given multiple times if it is overlapping several regions. 
      # A region can be presented, for example, in the following format: �chr2� (the whole chr2), 
      # chr2:1000000� (region starting from 1,000,000bp) or �chr2:1,000,000-2,000,000 
      # (region between 1,000,000 and 2,000,000bp including the end points). The coordinate is 1-based.
      # For iterating over reads on a stream, use BAMFile#foreach
      def self.view(bam_file, chr = nil, start = nil, stop = nil)
        if block_given?
          IO.popen("samtools view #{File.expand_path(bam_file)} #{query_string(chr, start, stop)}") do |pipe|
            pipe.each { |line| yield line }
          end
        else
          %x[ samtools view #{bam_file} #{query_string(chr, start, stop)} ]
        end
      end

      # Count the number of alignments overlapping the region
      def self.count(bam_file, chr, start = nil, stop = nil)
        %x[ samtools view -c #{File.expand_path(bam_file)} #{query_string(chr, start, stop)} ].chomp.to_i
      end

      # Sort alignments by leftmost coordinates. File <out.prefix>.bam will be created. 
      # This command may also create temporary files <out.prefix>.%d.bam when the whole 
      # alignment cannot be fitted into memory (controlled by option -m).
      def self.sort(input, output, max_mem = 500_000_000)
        output_prefix = File.basename(output, '.bam')
        %x[ samtools sort -m #{max_mem} #{File.expand_path(input)} #{output_prefix} ]
      end

      # Index sorted alignment for fast random access. Index file <aln.bam>.bai will be created.
      def self.index(bam_file)
        puts "Generating index for BAM file #{bam_file}" if ENV['DEBUG']
        %x[ samtools index #{File.expand_path(bam_file)} ]
      end

      # Retrieve and print stats in the index file. The output is TAB delimited with each 
      # line consisting of reference sequence name, sequence length, # mapped reads and # unmapped reads.
      def self.idxstats(bam_file)
        %x[ samtools idxstats #{File.expand_path(bam_file)} ]
      end

      # Merge multiple sorted alignments. The header reference lists of all the input BAM files, 
      # and the @SQ headers of inh.sam, if any, must all refer to the same set of reference sequences. 
      # The header reference list and (unless overridden by -h) ‘@’ headers of in1.bam will be copied 
      # to out.bam, and the headers of other files will be ignored.
      def self.merge(input_bams, output_bam)
        %x[ samtools merge #{File.expand_path(output_bam)} #{input_bams.map { |f| File.expand_path(f) }.join(' ')} ]
      end

      # Index reference sequence in the FASTA format or extract subsequence from indexed reference 
      # sequence. If no region is specified, faidx will index the file and create <ref.fasta>.fai 
      # on the disk. If regions are speficified, the subsequences will be retrieved and printed to 
      # stdout in the FASTA format. The input file can be compressed in the RAZF format.
      def self.faidx(fasta_file, intervals=nil)
        if regions.nil?
          %x[ samtools faidx #{File.expand_path(fasta_file)} ]
        else
          %x[ samtools faidx #{File.expand_path(fasta_file)} #{intervals.join(' ')} ]
        end
      end

      # Generate BCF or pileup for one or multiple BAM files. 
      # Alignment records are grouped by sample identifiers in @RG header lines. 
      # If sample identifiers are absent, each input file is regarded as one sample.
      def self.mpileup
        raise "Not yet implemented"
      end

      # Replace the header in in.bam with the header in in.header.sam. 
      # This command is much faster than replacing the header with a BAM->SAM->BAM conversion.
      def self.reheader(bam_file, header_sam_file)
        %x[ samtools reheader #{File.expand_path(header_sam_file)} #{File.expand_path(bam_file)} ]
      end

      # Remove potential PCR duplicates: if multiple read pairs have identical external coordinates, 
      # only retain the pair with highest mapping quality. In the paired-end mode, this command ONLY 
      # works with FR orientation and requires ISIZE is correctly set. It does not work for unpaired 
      # reads (e.g. two ends mapped to different chromosomes or orphan reads).
      def self.rmdup(input_bam, output_bam, single_end = false)
        %x[ samtools rmdup #{single_end ? '-S' : ''} #{File.expand_path(input_bam)} #{File.expand_path(output_bam)} ]
      end

      def self.query_string(chr = nil, start = nil, stop = nil)
        query = StringIO.new
        query << chr if chr
        query << ':' << start if start
        query << '-' << stop if stop
        
        return query.string
      end
      
    end
  end
end
