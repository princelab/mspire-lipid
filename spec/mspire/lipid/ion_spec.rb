require 'spec_helper'

require 'mspire/lipid'
require 'mspire/lipid/modification'
require 'mspire/lipid/ion'
require 'mspire/mass/element'
require 'mspire/mass/common'

module MSS
  CO2 = ( %w(C O O).map {|e| Mspire::Mass::Element::MONO_STRING[e] }.reduce(:+) )
  H2O = Mspire::Mass::Common::MONO_STRING['H2O']
  PROTON_LOSS = Mspire::Lipid::Modification.new(:proton, :loss => true)
end


describe Mspire::Lipid::Ion do

  before do
    lipid = Mspire::Lipid.new("LMGP02010009", "PE(16:0/18:1(9Z))", "1-hexadecanoyl-2-(9Z-octadecenoyl)-sn-glycero-3-phosphoethanolamine", 'C39H76NO8P', 717.5308, 'Glycerophospholipids [GP]', 'Glycerophosphoethanolamines [GP02]', 'Diacylglycerophosphoethanolamines [GP0201]', '7850611', 'FHQVHHIBKUMWTI-OTMQOFQLSA-N')
    proton = Mspire::Lipid::Modification.new(:proton)
    proton_loss = Mspire::Lipid::Modification.new(:proton, :loss => true)
    h2o_loss = Mspire::Lipid::Modification.new(:water, :loss => true)
    @plus1_less_h20 = Mspire::Lipid::Ion.new(lipid, [proton, h2o_loss])
    @plus2_less_h20 = Mspire::Lipid::Ion.new(lipid, [proton, proton, h2o_loss])
    @minus1_less_h20 = Mspire::Lipid::Ion.new(lipid, [proton_loss, h2o_loss])
  end

  it 'calculates the correct m/z' do
    @plus1_less_h20.mz.should be_within(1e-5).of(700.52751178307)
    @plus2_less_h20.mz.should be_within(1e-5).of(350.76739412492003)
    @minus1_less_h20.mz.should be_within(1e-5).of(698.5129588842301)
  end

  it 'calculates the correct formula' do
    @plus1_less_h20.formula.to_s.should == "C39H75NO7P"
    @plus1_less_h20.charge.should == 1

    @plus2_less_h20.formula.to_s.should == "C39H76NO7P"
    @plus2_less_h20.charge.should == 2
    @minus1_less_h20.formula.to_s.should == "C39H73NO7P"
    @minus1_less_h20.charge.should == -1
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

      # WORKING THIS GUY OUT!!!!!!
      xit '11-methyl-9S-hydroxy-13E-hexadecenoic acid' do
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
        #mzs.sort.should == frags.sort
      end

    end

    it 'predicts simple glycerophosopholipids' do
    end

  end

end
