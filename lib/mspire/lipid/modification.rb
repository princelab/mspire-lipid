require 'mspire/mass'
require 'mspire/molecular_formula'

module Mspire
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
        Mspire::Mass.formula_to_exact_mass(formula)
        massdiff = Mspire::Mass.formula_to_exact_mass(formula)
        massdiff -= (charge * Mspire::Mass::ELECTRON) # + charge subtracts, - charge adds
        massdiff = -massdiff unless gain
        massdiff
      end

      # the charge on the mod should be represented by the number of plusses
      # or minuses after the formula (Li+ for a +1 charge Lithium or H2++, 2
      # protons with a total of 2 charges)
      FORMULAS = {
        :proton => 'H',
        :ammonium => 'NH4',
        :lithium => 'Li',
        :sodium => 'Na',
        :water => 'H2O',
        :ammonia => 'NH3',
        :carbon_dioxide => 'CO2',
      }
      CHARGE = {
        :proton => 1,
        :ammonium => 1,
        :lithium => 1,
        :sodium=> 1,
        :water => 0,
        :ammonia => 0,
        :carbon_dioxide => 0,
      }

      # determined by running formulas through Mspire::Mass.massdiff
      MASSDIFFS = {}
      FORMULAS.each do |name, formula|
         MASSDIFFS[name] = self.massdiff(formula, CHARGE[name])
      end

      # as a symbol
      attr_accessor :name
      # a MolecularFormula object
      attr_accessor :formula
      # negative indicates a loss
      attr_accessor :massdiff
      # the charge 
      attr_accessor :charge

      # if no mass or formula is given then it searches command mods for the name
      # @param [Symbol] name the name of the mod
      # A number of opts are expected if they are not found in the FORMULAS,
      # CHARGE, or MASSDIFFS hashes.  However, the massdiff will be inferred
      # from the formula if it is not given:
      #
      #     attributes:
      #     :formula = the chemical formula, lipidmaps style ("C2H4BrO") or
      #                any valid argument to MolecularFormula.from_any
      #     :massdiff = +/-Float
      #     :charge = +/- Integer
      #
      #     instruction:
      #     :loss = true   negates the mass diff sign and charge during initialization
      #                    this option is typically only done for molecules
      #                    already present in the FORMULA hash (e.g.)
      #
      #     proton_loss = Mspire::Lipid::Modification.new(:proton, :loss => true)
      #     water_loss = Mspire::Lipid::Modification.new(:water, :loss => true)
      #
      def initialize(name, opts={})
        @name = name
        @formula = 
          if ( form_string = (opts[:formula] || FORMULAS[name]) )
            Mspire::MolecularFormula.from_any( form_string )
          end
        @massdiff = opts[:massdiff] || MASSDIFFS[name]
        @charge = opts[:charge] || CHARGE[name]

        if opts[:loss]
          @charge = -@charge
          # necessary if you are using a named molecule and you want its loss
          # rather than gain (i.e., you want a negative massdiff)
          @massdiff = -@massdiff 
        end
      end

      def charged_formula_string
        @formula.to_s + @charge.abs.times.map { (@charge > 0) ? '+' : '-' }.join
      end

      alias_method :to_s, :charged_formula_string

      def gain?
        massdiff > 0
      end

      def loss?
        !gain?
      end

      def inspect
        "<Mod: #{to_s}>"
      end

    end
  end
end


