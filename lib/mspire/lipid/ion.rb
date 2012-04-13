require 'mspire/lipid/ion/fragment'

module Mspire
  class Lipid
    # a lipid with modifications (typically the mods give it a charge so that
    # it can be seen in the mass spec)
    class Ion
      # an Mspire::Lipid object
      attr_accessor :lipid
      # an Mspire::Lipid::Modifications object
      attr_accessor :modifications
      # the key attribute of a query

      def initialize(lipid, mods=[])
        @lipid = lipid
        @modifications = mods
        @mz = nil
      end

      def charge
        z = 0
        @modifications.each do |mod|
          z -= mod.charge
        end
        z
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
