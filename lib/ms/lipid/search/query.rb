
module MS
  module Lipid
    class Search
      # a Lipid Search
      class PossibleLipid
        # an MS::Lipid object
        attr_accessor :lipid
        # an MS::Lipid::Modifications object
        attr_accessor :modifications
        # the key attribute of a query

        def mz
          return @mz if @mz
          mass = @lipid.mass
          charge = 0
          @modifications.each do |mod|
            mass = mod.gain ? mass + mod.mass : mass - mod.mass
            charge += mod.charge
          end
          @mz = mass / charge
        end
      end
    end
  end
end

