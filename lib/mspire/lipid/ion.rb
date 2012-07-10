require 'mspire/lipid/ion/fragment'
require 'mspire/molecular_formula'

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
          z += mod.charge
        end
        z
      end

      # a MolecularFormula object
      def formula
        _formula = @lipid.formula
        _formula = Mspire::MolecularFormula.from_any(_formula) unless _formula.is_a?(Mspire::MolecularFormula)
        modifications.each do |mod|
          if mod.gain?
            _formula += mod.formula 
          else
            _formula -= mod.formula
          end
        end
        _formula
      end

      # value is cached
      def mz_signed
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

      # the unsigned m/z value
      def mz
        _mz_signed = mz_signed
        _mz_signed >= 0 ? _mz_signed : -_mz_signed
      end

      def inspect
        "<|| Ion mz=#{mz} #{lipid.inspect} + #{modifications.map(&:inspect).join(', ')} ||>"
      end

    end
  end
end
