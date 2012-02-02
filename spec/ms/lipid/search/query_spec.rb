require 'spec_helper'

require 'ms/lipid/search/query'
require 'ms/lipid/modification'
require 'ms/lipid'

describe MS::Lipid::Search::Query do
  before do
    lipid = MS::Lipid.new
    lipid.mass = 300.2
    proton = MS::Lipid::Modification.new(:proton)
    h2o_loss = MS::Lipid::Modification.new(:water, :loss => true)
    @plus1 = MS::Lipid::Search::Query.new(lipid, [proton, h2o_loss])
    @plus2 = MS::Lipid::Search::Query.new(lipid, [proton, proton, h2o_loss])
  end

  it 'calculates the correct m/z' do
    @plus1.mz.should be_within(1e5).of(283.196711735)
    @plus2.mz.should be_within(1e5).of(142.101994085)
  end

end
