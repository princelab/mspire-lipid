require 'spec_helper'

require 'ms/lipid_maps'
require 'ms/lipid/search'
require 'ms/lipid/search/query'
require 'ms/lipid/modification'

describe MS::Lipid::Search do
  before do
    @lipids = MS::LipidMaps.parse_file(TESTFILES + '/lipidmaps_short.tsv')
    @queries = @lipids.map do |lipid| 
      proton = MS::Lipid::Modification.new(:proton)
      h2o_loss = MS::Lipid::Modification.new(:water, :loss => true)
      [[proton], [proton, h2o_loss]].map do |mods|
        MS::Lipid::Search::Query.new(lipid, mods)
      end
    end.flatten(1)
  end

  it 'creates a search spectrum' do
    spec = subject.create_search_spectrum(@queries)
    spec.mzs.any? {|mz| mz.nil? }.should be_false
    spec.mzs.size.should == 56
    spec.intensities.map(&:size).count(2).should == 4
    spec.intensities.map(&:size).count(1).should == 52
  end

  it 'creates a probability distribution' do
    subject.create_probability_function(@queries, :prob_bincnt => 20)
  end

end
