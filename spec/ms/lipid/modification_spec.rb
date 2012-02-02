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

    # proton gain
    prot = Mod.new(:proton)
    p prot
    prot.gain?.should be_true
    prot.massdiff.>(0).should be_true
    prot.charge.should == 1

    ammon = Mod.new(:ammonium)
    p ammon
    prot.gain?.should be_true
    prot.massdiff.>(0).should be_true
    prot.charge.should == 1

  end
  it 'behaves' do
  end
end
