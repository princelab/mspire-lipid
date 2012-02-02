require 'spec_helper'

require 'ms/lipid/modification'

describe MS::Lipid::Modification do
  Mod = MS::Lipid::Modification

  it 'can create common mods easily' do
    # water loss
    water_loss = Mod.new(:water, :loss => true)
    water_loss.loss?.should be_true
    water_loss.massdiff.<(0).should be_true
    water_loss.charge.should == 0
    water_loss.charged_formula.should == 'H2O'

    # proton gain
    prot = Mod.new(:proton)
    prot.gain?.should be_true
    prot.massdiff.>(0).should be_true
    prot.charge.should == 1
    prot.charged_formula.should == 'H+'

    ammon = Mod.new(:ammonium)
    ammon.gain?.should be_true
    ammon.massdiff.>(0).should be_true
    ammon.charge.should == 1
    ammon.charged_formula.should == 'NH3H+'
  end

  it 'can create custom mods' do
    mymod = Mod.new(:super_snazzy)
    mymod.formula.should be_nil
    mymod.massdiff.should be_nil
    mymod.charge.should be_nil

    mymod.formula = 'CH4'
    mymod.charge = 2
    mymod.massdiff = MS::Lipid::Modification.massdiff(mymod.formula, mymod.charge)
    mymod.massdiff.should be_within(1e4).of(16.030202)
  end
end
