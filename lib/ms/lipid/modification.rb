require 'ms/mass'

module MS
  class Lipid
    # the convention is all mods are gains unless the name ends in an
    # underscore
    class Modification
      FORMULAS = {
        :proton => 'H+',
        :ammonium => 'NH3H+',
        :lithium => 'Li+',
        :water => 'H2O',
      }
      # determined by running formulas through MS::Mass.formula_to_mass
      MASSDIFFS = {}
      FORMULAS.each do |name, formula|
        MASSDIFFS[name] = MS::Mass.formula_to_exact_mass(FORMULAS[name])
      end

      # as a symbol
      attr_accessor :name
      # as a molecular formula
      attr_accessor :formula
      # negative indicates a loss
      attr_accessor :massdiff
      # the charge 
      attr_accessor :charge

      # if no mass or formula is given then it searches command mods for the name
      # @param [Symbol] name the name of the mod
      # A number of opts are allowed:
      #
      #     :formula = the chemical formula, lipidmaps style ("C2H4BrO")
      #     :mass = Float
      #     :charge = +/- Integer
      #     :loss = true   will flip the sign of the mass difference
      def initialize(name, opts={})
        @name = name
        @formula = opts[:formula] || FORMULAS[name]
        @massdiff = opts[:massdiff] || MASSDIFFS[name]
        @massdiff = -@massdiff if opts[:loss]
      end

      def charge
        @formula.rindex .as.dfklas.dfasdf
      end

      def gain?
        massdiff > 0
      end

      def loss?
        !gain?
      end

    end
  end
end


