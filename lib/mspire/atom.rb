
module Mspire
  class Atom
    COMMON_VALENCE = {
      :c => 4,
      :o => 2,
      :h => 1,
      :n => 3,
      :p => 5, # 3 ??
    }

    attr_accessor :element
    attr_accessor :valence
    attr_accessor :bonds

    def initialize(element, bonds=[], valence=nil)
      @element = element.downcase.to_sym
      @bonds = bonds
      @valence = 
        if valence
          valence
        else
          COMMON_VALENCE[element]
        end
    end
  end
end
