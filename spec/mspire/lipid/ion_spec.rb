require 'spec_helper'


require 'mspire/lipid'
require 'mspire/lipid/modification'
require 'mspire/lipid/ion'

module MSS
  CO2 = ( %w(c o o).map {|e| Mspire::Mass::MONO[e] }.reduce(:+) )
  H2O = Mspire::Mass::MONO['h2o']
  PROTON_LOSS = Mspire::Lipid::Modification.new(:proton, :loss => true)
end


describe Mspire::Lipid::Ion do

  before do
    lipid = Mspire::Lipid.new
    lipid.mass = 300.2
    proton = Mspire::Lipid::Modification.new(:proton)
    h2o_loss = Mspire::Lipid::Modification.new(:water, :loss => true)
    @plus1 = Mspire::Lipid::Ion.new(lipid, [proton, h2o_loss])
    @plus2 = Mspire::Lipid::Ion.new(lipid, [proton, proton, h2o_loss])
  end

  it 'calculates the correct m/z' do
    @plus1.mz.should be_within(1e5).of(283.196711735)
    @plus2.mz.should be_within(1e5).of(142.101994085)
  end

  describe 'predicting ms/ms fragments' do
    describe 'predicting simple fatty acyls' do
      xit 'hexadecanoic acid' do
        lipid = Mspire::Lipid.new('LMFA01010001', 'Palmitic acid', 'hexadecanoic acid',	'C16H32O2', 256.24, 'Fatty Acyls [FA]', 'Fatty Acids and Conjugates [FA01]', 'Straight chain fatty acids [FA0101]')

        ion = Mspire::Lipid::Ion.new(lipid, [MSS::PROTON_LOSS])
        frags = [ion.mz - MSS::H2O]
        frags << ion.mz - MSS::CO2
        frags << ion.mz - (MSS::H2O + MSS::CO2)

        mzs = ion.predict_fragment_mzs
        mzs.should be_an(Array)
        mzs.sort.should == frags.sort
      end

      it '11-methyl-9S-hydroxy-13E-hexadecenoic acid' do
        lipid = Mspire::Lipid.new(nil, 'Made Up', '11-methyl-9S-hydroxy-13E-hexadecenoic acid',	'C17H32O3', 284.23514488492003, 'Fatty Acyls [FA]', 'Fatty Acids and Conjugates [FA01]', 'Branched fatty acids [FA0102]')
        ion = Mspire::Lipid::Ion.new(lipid, [MSS::PROTON_LOSS])

        # one water loss
        frags = [ion.mz - MSS::H2O]

        # two water loss
        frags << (ion.mz - (2*MSS::H2O))

        # CO2 loss
        frags << (ion.mz - MSS::CO2)

        # CO2 + H20 loss
        frags << ion.mz - (MSS::H2O + MSS::CO2)

        oh_frag_mass = 140.12011513268
        frag1 = (ion.mz - oh_frag_mass)

        frags << frag1 - MSS::H2O

        frags << frag1 - MSS::CO2

        frags << frag1 - (MSS::H2O + MSS::CO2)

        mzs = ion.predict_fragment_mzs
        p mzs
        #mzs.sort.should == frags.sort
      end

    end

    it 'predicts simple glycerophosopholipids' do
    end

  end

end
