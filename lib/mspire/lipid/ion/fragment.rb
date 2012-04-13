
module Mspire
  class Lipid

    # goes from 1 to 99
    CHAIN_PREFIXES = {
      'meth' => 1,
      'eth' => 2,
      'prop' => 3,
      'but' => 4,
      'pent' => 5,
      'hex' => 6,
      'hept' => 7,
      'oct' => 8,
      'non' => 9,
      'dec' => 10,
      'undec' => 11,
      'dodec' => 12,
      'tridec' => 13,
      'tetradec' => 14,
      'pentadec' => 15,
      'hexadec' => 16,
      'heptadec' => 17,
      'octadec' => 18,
      'nonadec' => 19,
      'eicos' => 20,
      'heneicos' => 21,
      'docos' => 22,
      'tricos' => 23,
      'tetracos' => 24,
      'pentacos' => 25,
      'hexacos' => 26,
      'heptacos' => 27,
      'octacos' => 28,
      'nonacos' => 29
    }

    consistent = { 
      0 => '',
      1 => 'hen',
      2 => 'do',
      3 => 'tri',
      4 => 'tetra',
      5 => 'penta',
      6 => 'hexa',
      7 => 'hepta',
      8 => 'octa',
      9 => 'nona',
    }

    (3..9).each do |tens_place|
      (0..9).each do |ones_place| 
        key = consistent[ones_place] + consistent[tens_place] + "cont"
        CHAIN_PREFIXES[key] = 10*tens_place + ones_place
      end
    end

    class Ion
      module Fragment
        # predicts the MS/MS fragments for this ion
        def predict_fragment_mzs
          lipid.smiles
        end

      end

      include Fragment
    end
  end
end
