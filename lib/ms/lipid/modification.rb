
module MS
  class Lipid
    class Modification
      FORMULAS = {
        :proton => 'H+',
        :ammonium => 'NH3H+',
        :lithium => 'Li',
        :water => 'H2O',
      }
      # determined by running formulas through MS::Mass.formula_to_mass
      MASSDIFFS = {}
      FORMULAS.each do |name, formula|
        MASSDIFFS[name] = MS::Mass.formula_to_mass(FORMULAS[name])
      end
      CHARGE = {
        :proton => 1,
        :ammonium => 1,
        :lithium => 1,
        :water => 0,
      }
      GAIN = {
        :proton => true, 
        :ammonium => true,
        :lithium => true,
        :water => false,
      }
      attr_accessor :name
      attr_accessor :formula
      attr_accessor :mass
      # the charge 
      attr_accessor :charge
      # if no mass or formula is given then it searches command mods for the name
      # @param [Symbol] name the name of the mod
      # A number of opts are allowed:
      #
      #     :formula = the chemical formula, lipidmaps style ("C2H4BrO")
      #     :mass = Float
      #     :charge = +/- Integer
      #     :gain = boolean # is this a mass lost or gain
      def initialize(name, opts={})
        @name = name
        @formula = opts[:formula] || FORMULAS[name]
        @massdiff = opts[:massdiff] || MASSDIFFS[name]
        @charge = opts[:charge] || CHARGE[name]
        @gain = opts[:gain] || GAIN[name]
      end
    end
  end
end


