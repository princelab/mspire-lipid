
require 'matrix'

module Mspire
  class SDF
    # nearly identical to Bio::PDB::Coordinate
    class Coordinate < Vector
      #http://bioruby.org/rdoc/Bio/PDB/Coordinate.html

      def self.[](_x,_y,_z)
        super
      end

      def distance(obj)
        self.-(obj).r
      end

      def x() self[0] end
      def x=(val) self[0]=val end

      def y() self[1] end
      def y=(val) self[1]=val end

      def z() self[2] end
      def z=(val) self[2]=val end

      def inspect
        "#{self.class}[#{x}, #{y}, #{z}]"
      end
    end

    COMMON_VALENCE = {
      :c => 4,
      :o => 2,
      :h => 1,
      :n => 3,
      :p => 5, # 3 ??
    }

    Atom = Struct.new(:element, :bonds, :coordinates, :valence) do
      # element should be all lowercase symbol
      def initialize(_element=:h, _bonds=[], _coordinates=Vector[0.0, 0.0, 0.0], _valence=nil)
        super( _element, _bonds, _coordinates, _valence || COMMON_VALENCE[_element])
      end

      # returns valence - bonds.size
      def empty_bonds
        valence - bonds.size
      end

      # if break_count is nil, breaks every bond between the atoms, otherwise
      # breaks break_count number of bonds between the atoms.  Returns the
      # number of bonds broken.
      def break(other_atom, break_count=nil)
        break_count ||= 1e10 
        new_bonds = []
        cnt = 0
        bonds.each do |bond|
          if bond.include?(other_atom)
            cnt += 1
          else
            new_bonds << bond
          end
        end
        bonds.replace(new_bonds)
        cnt
      end

      # returns self
      def add(other_atom, bond_cnt=1)
        bond = Bond.new([self, other_atom], bond_cnt)
        self.bonds << bond
        other_atom.bonds << bond
        self
      end

      def inspect
        "#<struct #{self.class} element=#{element.inspect}, bonds.size=#{bonds.size}, coordinates=#{coordinates.inspect}, valence=#{valence.inspect}>"
      end
    end

    class Bond < Array
      # count is the type of bond (1 single, 2 double, etc)
      def initialize(atoms=[], count=1)
        super(_atoms)
        @count = count
      end

      def smiles_type
        case count
        when 1
          '-'
        when 2
          '='
        when 3
          '#'
        else
          'U'
        end
      end

      def inspect
        "<#{self.map(&:object_id).join('--')} (#{self.map(&:element).join(smiles_type)})>"
      end
    end

    # (4th line) an array of count data
    attr_accessor :count_data
    # an array: the first 3 lines
    attr_accessor :header

    attr_accessor :atoms
    attr_accessor :bonds

    # the string for one compound
    def initialize(string=nil)
      parse!(string) if string
    end

    def parse!(string)
      lines = string.split(/\r?\n/)
      @header = lines[0,3]
      count_line = lines[3]
      @count_data = count_line.split(/\s+/)
      @atoms = lines[4,num_atoms].map do |line|
        data = line.split(/\s+/)
        Atom.new(data[3].downcase.to_sym, [], Coordinate[*data[0,3].map(&:to_f)])
      end
      atom_index_to_bonds = Hash.new {|h,k| h[k] = [] }

      @bonds = lines[num_atoms+4,num_bonds].map do |line|
        data = line.split(/\s+/).map(&:to_i)

        # using indices instead of the atom itself because the atom is a
        # struct and the hashing algorithm is quite extensive for structs
        atom_indices = data[0,2].map {|num| num-1 }
        bond = Bond.new(atom_indices.map {|i| @atoms[i] }, data[2] )

        atom_indices.each do |ai|
          atom_index_to_bonds[ai] << bond
        end
        bond
      end

      atom_index_to_bonds.each do |atom_i, bonds|
        atom = @atoms[atom_i]
        bonds.each do |bond|
          bond.atoms.delete(atom)
          @atoms[atom_i].atoms << bond.atoms.first
        end
      end
      self
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
