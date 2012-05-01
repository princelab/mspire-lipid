#require 'mspire/sdf'
require 'mspire/molecular_formula'
require 'openbabel'

module Mspire
  class Molecule
    include Enumerable

    #http://openbabel.org/docs/2.3.1/FileFormats/Overview.html
    FORMAT_DESCRIPTIONS = {
      CONTCAR:     'VASP format',
      POSCAR:      'VASP format',
      acr:         'ACR format',
      adfout:      'ADF output format',
      alc:         'Alchemy format',
      arc:         'Accelrys/MSI Biosym/Insight II CAR format',
      bgf:         'MSI BGF format',
      box:         'Dock 3.5 Box format',
      bs:          'Ball and Stick format',
      c3d1:        'Chem3D Cartesian 1 format',
      c3d2:        'Chem3D Cartesian 2 format',
      caccrt:      'Cacao Cartesian format',
      can:         'Canonical SMILES format.',
      car:         'Accelrys/MSI Biosym/Insight II CAR format',
      ccc:         'CCC format',
      cdx:         'ChemDraw binary format',
      cdxml:       'ChemDraw CDXML format',
      cif:         'Crystallographic Information File',
      ck:          'ChemKin format',
      cml:         'Chemical Markup Language',
      cmlr:        'CML Reaction format',
      crk2d:       'Chemical Resource Kit diagram(2D)',
      crk3d:       'Chemical Resource Kit 3D format',
      ct:          'ChemDraw Connection Table format',
      cub:         'Gaussian cube format',
      cube:        'Gaussian cube format',
      dat:         'Generic Output file format',
      dmol:        'DMol3 coordinates format',
      dx:          'OpenDX cube format for APBS',
      ent:         'Protein Data Bank format',
      fch:         'Gaussian formatted checkpoint file format',
      fchk:        'Gaussian formatted checkpoint file format',
      fck:         'Gaussian formatted checkpoint file format',
      feat:        'Feature format',
      fract:       'Free Form Fractional format',
      fs:          'FastSearching',
      g03:         'Gaussian Output',
      g09:         'Gaussian Output',
      g92:         'Gaussian Output',
      g94:         'Gaussian Output',
      g98:         'Gaussian Output',
      gal:         'Gaussian Output',
      gam:         'GAMESS Output',
      gamess:      'GAMESS Output',
      gamin:       'GAMESS Input',
      gamout:      'GAMESS Output',
      gpr:         'Ghemical format',
      gukin:       'GAMESS-UK Input',
      gukout:      'GAMESS-UK Output',
      gzmat:       'Gaussian Z-Matrix Input',
      hin:         'HyperChem HIN format',
      inchi:       'InChI format',
      inp:         'GAMESS Input',
      ins:         'ShelX format',
      jout:        'Jaguar output format',
      log:         'Generic Output file format',
      mcdl:        'MCDL format',
      mcif:        'Macromolecular Crystallographic Information',
      mdl:         'MDL MOL/SDF format',
      ml2:         'Sybyl Mol2 format',
      mmcif:       'Macromolecular Crystallographic Information',
      mmd:         'MacroModel format',
      mmod:        'MacroModel format',
      mol:         'MDL MOL/SDF format',
      mol2:        'Sybyl Mol2 format',
      mold:        'Molden format',
      molden:      'Molden format',
      moo:         'MOPAC Output format',
      mop:         'MOPAC Cartesian format',
      mopcrt:      'MOPAC Cartesian format',
      mopin:       'MOPAC Internal',
      mopout:      'MOPAC Output format',
      mpc:         'MOPAC Cartesian format',
      mpo:         'Molpro output format',
      mpqc:        'MPQC output format',
      msi:         'Accelrys/MSI Cerius II MSI format',
      nwo:         'NWChem output format',
      out:         'Generic Output file format',
      outmol:      'DMol3 coordinates format',
      output:      'Generic Output file format',
      pc:          'PubChem format',
      pcm:         'PCModel Format',
      pdb:         'Protein Data Bank format',
      png:         'PNG files with embedded data',
      png2:         'PNG files with embedded data (VERSION 2??)',
      pqr:         'PQR format',
      pqs:         'Parallel Quantum Solutions format',
      prep:        'Amber Prep format',
      qcout:       'Q-Chem output format',
      res:         'ShelX format',
      rsmi:        'Reaction SMILES format',
      rxn:         'MDL RXN format',
      sd:          'MDL MOL/SDF format',
      sdf:         'MDL MOL/SDF format',
      smi:         'SMILES format',
      smiles:      'SMILES format',
      sy2:         'Sybyl Mol2 format',
      t41:         'ADF TAPE41 format',
      tdd:         'Thermo format',
      therm:       'Thermo format',
      tmol:        'TurboMole Coordinate format',
      txt:         'Title format',
      unixyz:      'UniChem XYZ format',
      vmol:        'ViewMol format',
      xml:         'General XML format',
      xtc:         'XTC format',
      xyz:         'XYZ cartesian coordinates format',
      yob:         'YASARA.org YOB format',
      svg:         'scalable vector graphics', 
      pov:         'pov-ray format', 
      fix:         'SMILES FIX format',
      fasta:         'FASTA format',
      fa:         'FASTA format',
      fsa:         'FASTA format',
    }

    def self.from_lipidmaps_sdf_string(string)
      mol = self.new
      mol.read_string(string.split('|').join("\n"), :sdf)
      mol
    end

    def self.type_from_file(file)
      type = File.extname(file)[1..-1].to_sym
      unless FORMAT_DESCRIPTIONS.keys.include?(type)
        raise(ArgumentError, "unsupported format: #{type}") 
      end
      type
    end

    # the low-level molecule object
    attr_accessor :obmol
    # the object controlling input and output
    attr_accessor :obconv

    def initialize
      @obconv = OpenBabel::OBConversion.new
      @obmol = OpenBabel::OBMol.new
    end

    def each_atom(&block)
      block or return enum_for(__method__)
      (1..@obmol.num_atoms).each do |n|
        block.call( @obmol.get_atom(n) )
      end
    end
    alias_method :each, :each_atom

    def each_bond(&block)
    end

    def mass
      @obmol.get_exact_mass
    end

    def charge
      @obmol.get_total_charge
    end

    def molecular_formula
      MolecularFormula.new(@obmol.get_formula, charge)
    end

    def fragment_on(type)
      send("fragment_#{type}")
    end

    # returns an array OpenBabel::OBAtom objects that match the string
    def smart(smart_string)
    end

    def fragment_ether
    end

    # level 1 enumerates every possible fragment with only one break point
    def fragment(level=1)

    end

    # sets the molecule from file and returns self
    def read_file(filename)
      # consider using ext_from_filename in future
      type = self.class.type_from_file(filename)
      @obconv.set_in_format(type.to_s)
      @obconv.read_file(@obmol, filename)
      self
    end

    # sets the molecule from string and returns self
    def read_string(string, type=:smi)
      @obconv.set_in_format(type.to_s)
      @obconv.read_string(@obmol, string)
      self
    end

    # outputs canonical SMILES format by default (with no name!)
    def to_s(type=:can)
      string = write_string(type)
      case type
      when :smi, :smiles, :can
        # remove name with options in the future
        string.split(/\s+/).first
      else
        string
      end
    end

    def write_string(type=:can)
      @obconv.set_out_format(type.to_s)
      @obconv.write_string(@obmol)
    end

    # uses the extension to determine the type, returns num bytes written
    def write_file(filename)
      type = self.class.type_from_file(filename)
      # later, could do this using the write_file method
      IO.write(filename, to_s(type))
    end

    def method_missing(*args, &block)
      if @obmol.respond_to?(args.first)
        @obmol.send(*args, &block)
      else
        super(*args, &block)
      end
    end
  end
end

=begin
    # an array of atoms that know how they are bound together
    attr_accessor :atoms

    # charge that is not localized to a particular atom
    attr_accessor :delocalized_charge

    # the element used to fill any remaining valence locations (typically
    # hydrogen (:h))
    attr_accessor :fill_valence

    def initialize(atoms=[], delocalized_charge=0, fill_valence=:h)
      @atoms = atoms
      @delocalized_charge = delocalized_charge
      @fill_valence = fill_valence
    end

    # sum of the charge on individual atoms + any delocalized charge
    def charge
      atoms.reduce(0) {|sum,atom| sum + atom.charge } + delocalized_charge
    end

    def self.from_lipidmaps_sdf_string(sdf_string)
      sdf = Mspire::SDF.new(sdf_string.gsub(/\s*\|\s*/, "\n"))
      self.new( sdf.atoms )
    end

    # if fill_valence is nil, then no extra atoms are intuited.
    def molecular_formula
      mf = Hash.new {|h,k| h[k] = 0 }
      _charge = 0
      atoms.each do |atom|
        # calculate the charge during the determination of elements 
        _charge += atom.charge
        mf[atom.element] += 1
        if fill_valence
          mf[fill_valence] += (1 * atom.empty_bonds)
        end
      end
      Mspire::MolecularFormula.new( mf, _charge + @delocalized_charge )
    end

    def mass
      molecular_formula.mass
    end
=end


=begin
"LMFA01010012|  LIPDMAPS01251212252D|| 14 13  0  0  0  0  0  0  0  0999 V2000|    6.4289    5.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|    7.8578    5.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|    7.1433    5.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|    8.5723    5.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|    9.2868    5.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|   10.0013    5.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|   10.7158    5.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|   11.4302    5.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|   12.1448    5.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|   12.8592    5.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|   13.5737    5.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0|   12.8592    6.2375    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0|    5.7144    5.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|    5.0000    5.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|  2  3  1  0  0  0  0|  3  1  1  0  0  0  0|  4  2  1  0  0  0  0|  5  4  1  0  0  0  0|  6  5  1  0  0  0  0|  7  6  1  0  0  0  0|  8  7  1  0  0  0  0|  9  8  1  0  0  0  0| 10  9  1  0  0  0  0| 11 10  1  0  0  0  0| 10 12  2  0  0  0  0| 13  1  1  0  0  0  0| 14 13  1  0  0  0  0|M  END"
=end

