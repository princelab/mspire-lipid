
module Mspire
  class SDF
    Atom = Struct.new

    # (4th line) an array of count data
    attr_accessor :count_data
    # an array: the first 3 lines
    attr_accessor :header

    # the string for one compound
    def initialize(string=nil)
      parse(string) if string
    end

    def parse(string)
      string.chomp!
      lines = string.split(/\r?\n/)
      @header = lines[0,3]
      count_line = lines[3]
      @count_data = count_line.split(/\s+/)
      @coordinates = lines[4,num_atoms]

    end

    def num_atoms
      @count_data[0].to_i
    end

    def num_bonds
      @count_data[1].to_i
    end

  end
end

=begin
 benzene
 ACD/Labs0812062058
 
  6  6  0  0  0  0  0  0  0  0  1 V2000
    1.9050   -0.7932    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.9050   -2.1232    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7531   -0.1282    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7531   -2.7882    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3987   -0.7932    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3987   -2.1232    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  0  0  0  0
  3  1  2  0  0  0  0
  4  2  2  0  0  0  0
  5  3  1  0  0  0  0
  6  4  1  0  0  0  0
  6  5  2  0  0  0  0
 M  END
=end

## note $$$$ is the split for multiple compounds
