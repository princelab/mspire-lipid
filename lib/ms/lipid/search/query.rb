
module MS
  class Lipid
    class Search
      # a Lipid Query (yields a specific m/z value)
      class Query
        # an MS::Lipid object
        attr_accessor :lipid
        # an MS::Lipid::Modifications object
        attr_accessor :modifications
        # the key attribute of a query

        def initialize(lipid, mods=[])
          @lipid = lipid
          @modifications = mods
        end

        def mz
          return @mz if @mz
          mass = @lipid.mass
          charge = 0
          @modifications.each do |mod|
            mass += mod.massdiff 
            charge += mod.charge
          end
          if charge == 0
            @mz = nil
          else
            @mz = mass / charge
          end
        end
      end
    end
  end
end

