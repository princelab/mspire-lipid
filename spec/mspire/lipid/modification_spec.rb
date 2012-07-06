require 'spec_helper'

require 'mspire/lipid/modification'
require 'mspire/lipid'
require 'mspire/lipid/ion'

describe Mspire::Lipid::Modification do
  Mod = Mspire::Lipid::Modification

  describe 'creates common mods easily' do

    it 'water loss' do
      water_loss = Mod.new(:water, :loss => true)
      water_loss.loss?.should be_true
      water_loss.massdiff.<(0).should be_true
      water_loss.charge.should == 0
      water_loss.charged_formula.should == 'H2O'
    end

    it 'proton gain' do
      prot = Mod.new(:proton)
      prot.gain?.should be_true
      prot.massdiff.>(0).should be_true
      prot.charge.should == 1
      prot.charged_formula.should == 'H+'
      ion = Mspire::Lipid::Ion.new(@lipid, [prot])
      p ion.formula
      p ion.mass
    end

    it 'proton loss' do
      prot_loss = Mod.new(:proton, :loss => true)
      prot_loss.gain?.should be_false
      prot_loss.loss?.should be_true
      prot_loss.massdiff.<(0).should be_true
      prot_loss.charge.should == -1
      prot_loss.charged_formula.should == 'H-'
      ion = Mspire::Lipid::Ion.new(@lipid, [prot_loss])
      p ion.formula
      p ion.mass
    end

    it 'ammonium gain' do
      ammon = Mod.new(:ammonium)
      ammon.gain?.should be_true
      ammon.massdiff.>(0).should be_true
      ammon.charge.should == 1
      ammon.charged_formula.should == 'H4N+'
    end
  end

  it 'can create custom mods' do
    mymod = Mod.new(:super_snazzy)
    mymod.formula.should be_nil
    mymod.massdiff.should be_nil
    mymod.charge.should be_nil

    mymod.formula = 'CH4'
    mymod.charge = 2
    mymod.massdiff = Mspire::Lipid::Modification.massdiff(mymod.formula, mymod.charge)
    mymod.massdiff.should be_within(1e4).of(16.030202)
  end
end
