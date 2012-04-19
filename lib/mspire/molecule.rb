require 'mspire/sdf'
require 'mspire/molecular_formula'

module Mspire
  class Molecule
    # an array of atoms that know how they are bound together
    attr_accessor :atoms
    attr_accessor :charge

    def initialize(atoms=[], charge=0)
      @atoms = atoms
      @charge = charge
    end

    def self.from_lipidmaps_sdf_string(sdf_string)
      sdf = Mspire::SDF.new(sdf_string.gsub(/\s*\|\s*/, "\n"))
      self.new( sdf.atoms )
    end

    def molecular_formula
      mf = Hash.new {|h,k| h[k] = 0 }
      atoms.each do |atom|
        mf[atom.element] += 1
      end
      Mspire::MolecularFormula.new( mf, charge )
    end

    def mass
      molecular_formula.mass
    end

  end
end

=begin
"LMFA01010012|  LIPDMAPS01251212252D|| 14 13  0  0  0  0  0  0  0  0999 V2000|    6.4289    5.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|    7.8578    5.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|    7.1433    5.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|    8.5723    5.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|    9.2868    5.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|   10.0013    5.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|   10.7158    5.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|   11.4302    5.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|   12.1448    5.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|   12.8592    5.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|   13.5737    5.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0|   12.8592    6.2375    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0|    5.7144    5.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|    5.0000    5.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|  2  3  1  0  0  0  0|  3  1  1  0  0  0  0|  4  2  1  0  0  0  0|  5  4  1  0  0  0  0|  6  5  1  0  0  0  0|  7  6  1  0  0  0  0|  8  7  1  0  0  0  0|  9  8  1  0  0  0  0| 10  9  1  0  0  0  0| 11 10  1  0  0  0  0| 10 12  2  0  0  0  0| 13  1  1  0  0  0  0| 14 13  1  0  0  0  0|M  END"
=end

