require 'spec_helper'

require 'ms/lipid_maps'
require 'ms/lipid/search'
require 'ms/lipid/search/query'
require 'ms/lipid/modification'

describe MS::Lipid::Search do
  before do
    @proton = MS::Lipid::Modification.new(:proton)
    @h2o_loss = MS::Lipid::Modification.new(:water, :loss => true)
  end
  describe 'searching a section of lipid maps' do
    before do
      @lipids = MS::LipidMaps.parse_file(TESTFILES + '/lipidmaps_short.tsv')
      @queries = @lipids.map do |lipid| 
        [[@proton], [@proton, @h2o_loss]].map do |mods|
          MS::Lipid::Search::Query.new(lipid, mods)
        end
      end.flatten(1)
    end

    xit 'creates a search spectrum' do
      spec = subject.create_search_spectrum(@queries)
      spec.mzs.any? {|mz| mz.nil? }.should be_false
      spec.mzs.size.should == 56
      spec.intensities.map(&:size).count(2).should == 4
      spec.intensities.map(&:size).count(1).should == 52
    end

    xit 'creates a probability distribution' do
      subject.create_probability_function(@queries, :prob_bincnt => 20)
    end
  end

  describe 'searching a full lipid maps' do

    before do
      # this will be specific to your install since it's not part of install
      path_to_lipidmaps_db = "#{ENV['HOME']}/tmp/tamil/lipidmaps_20120103_classes_1_2_3_4_5_6_7_8.exact_mass.tsv"
      @lipids = MS::LipidMaps.parse_file(path_to_lipidmaps_db)
      @queries = @lipids.map do |lipid| 
        [[@proton], [@proton, @proton], [@proton, @h2o_loss]].map do |mods|
          MS::Lipid::Search::Query.new(lipid, mods)
        end
      end.flatten(1)
    end

    it 'creates a probability distribution' do
      @queries.size
      subject.create_probability_function(@queries, :prob_bincnt => 1000)
    end

  end

end
