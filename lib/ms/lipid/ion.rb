
module MS
  class Lipid
    # a lipid with modifications (typically the mods give it a charge so that
    # it can be seen in the mass spec)
    class Ion
      # an MS::Lipid object
      attr_accessor :lipid
      # an MS::Lipid::Modifications object
      attr_accessor :modifications
      # the key attribute of a query

      def initialize(lipid, mods=[])
        @lipid = lipid
        @modifications = mods
        @mz = nil
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

      def inspect
        "<|| Ion mz=#{mz} #{lipid.inspect} + #{modifications.map(&:inspect).join(', ')} ||>"
      end
    end
  end
end
