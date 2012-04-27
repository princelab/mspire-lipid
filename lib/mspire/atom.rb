
module Mspire
  class Atom

    def self.[](*args)
      new(*args)
    end

    COMMON_VALENCE = {
      :c => 4,
      :o => 2,
      :s => 2, # can do 4 or 6 also
      :h => 1,
      :n => 3,
      :p => 5, # can do 3 also
    }
    attr_accessor :element, :bonds, :valence

    def initialize(element=:h, bonds=[], valence=nil)
      @element = element
      @bonds = bonds
      @valence = valence || COMMON_VALENCE[@element] 
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

  end
end
