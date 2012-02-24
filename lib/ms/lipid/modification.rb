require 'ms/mass'

module MS
  class Lipid


    # the convention is all mods are gains unless the name ends in an
    # underscore
    class Modification

      # given a string with a formula and charge, returns the formula portion
      # and the charges (as a signed integer)
      def self.formula_and_charge(string)
        md = string.match(/([^+]*)(\+*)$/)
        charges_string = md[2]
        if charges_string.nil?
          0
        else
          charges_string.count(charges_string[0])
          int = -int if charges_string[0] == '-'
        end
        [md[1], int]
      end

      # calculates the mass diff.  For every positive charge the mass of an
      # electron is subtracted; for every negative charge the mass of an
      # electron is added.  If gain is false, then the mass diff will be
      # negative.
      def self.massdiff(formula, charge, gain=true)
        MS::Mass.formula_to_exact_mass(formula)
        massdiff = MS::Mass.formula_to_exact_mass(formula)
        massdiff -= (charge * MS::Mass::ELECTRON) # + charge subtracts, - charge adds
        massdiff = -massdiff unless gain
        massdiff
      end

      # the charge on the mod should be represented by the number of plusses
      # or minuses after the formula (Li+ for a +1 charge Lithium or H2++, 2
      # protons with a total of 2 charges)
      FORMULAS = {
        :proton => 'H',
        :ammonium => 'NH3H',
        :lithium => 'Li',
        :water => 'H2O',
      }
      CHARGE = {
        :proton => 1,
        :ammonium => 1,
        :lithium => 1,
        :water => 0,
      }

      # determined by running formulas through MS::Mass.massdiff
      MASSDIFFS = {}
      FORMULAS.each do |name, formula|
         MASSDIFFS[name] = self.massdiff(formula, CHARGE[name])
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
      # A number of opts are expected if they are not found in the FORMULAS,
      # CHARGE, or MASSDIFFS hashes:
      #
      #     attributes:
      #     :formula = the chemical formula, lipidmaps style ("C2H4BrO")
      #     :massdiff = +/-Float
      #     :charge = +/- Integer
      #
      #     instruction:
      #     :loss = true   flips the mass diff sign during initialization
      #                    necessary to get negative massdiff on named molecule
      #                    (unnecessary if you input massdiff manually)
      def initialize(name, opts={})
        @name = name
        @formula = opts[:formula] || FORMULAS[name]
        @massdiff = opts[:massdiff] || MASSDIFFS[name]
        @charge = opts[:charge] || CHARGE[name]
        # necessary if you are using a named molecule and you want its loss
        # rather than gain (i.e., you want a negative massdiff)
        @massdiff = -@massdiff if opts[:loss]
      end

      def charged_formula
        @formula + @charge.abs.times.map { (@charge > 0) ? '+' : '-' }.join
      end

      def gain?
        massdiff > 0
      end

      def loss?
        !gain?
      end

      def inspect
        "<Mod: #{charged_formula}>"
      end

    end
  end
end


