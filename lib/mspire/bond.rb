
module Mspire
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
end
