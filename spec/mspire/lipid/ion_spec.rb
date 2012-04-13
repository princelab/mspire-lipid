require 'spec_helper'


require 'mspire/lipid'
require 'mspire/lipid/modification'
require 'mspire/lipid/ion'

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

end
