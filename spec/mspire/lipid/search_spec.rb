require 'spec_helper'

require 'mspire/lipid_maps'
require 'mspire/lipid/search'
require 'mspire/lipid/search/query'
require 'mspire/lipid/modification'

describe Mspire::Lipid::Search do
  before do
    @proton = Mspire::Lipid::Modification.new(:proton)
    @h2o_loss = Mspire::Lipid::Modification.new(:water, :loss => true)
  end
  describe 'searching a section of lipid maps' do
    before do
      @lipids = Mspire::LipidMaps.parse_file(TESTFILES + '/lipidmaps_programmatic_short.tsv')
      @ions = @lipids.map do |lipid| 
        [[@proton], [@proton, @h2o_loss]].map do |mods|
          Mspire::Lipid::Ion.new(lipid, mods)
        end
      end.flatten(1)
      @samples = Hash[ {
        :sample1 => [[187.1633, 244.22, 616.51, 717.50], 
          [100, 200, 100, 200]],
        :sample2 => [[187.164, 396.15, 244.24, 347.28, 618.502],
          [110, 210, 110, 210, 110]],
        :sample3 => [[187.160, 396.28, 244.24, 347.263, 618.511],
          [120, 220, 120, 220, 120]],
        :sample4 => [[187.157, 396.20, 244.30, 618.22, 933.01],
          [30, 33, 38, 99, 22]],
      }.map {|key,data| [key, Mspire::Spectrum.new(data)] } ]
      @pretend_search_mzs = [187.157, 396.20, 244.30, 618.22, 933.01]
    end

    xit 'creates a query search spectrum' do
      #spec = .create_query_search_spectrum(@ions)
      #spec.mzs.any? {|mz| mz.nil? }.should be_false
      #spec.mzs.size.should == 56
      #spec.intensities.map(&:size).count(2).should == 4
      #spec.intensities.map(&:size).count(1).should == 52
    end

    xit 'creates a probability function' do
      #subject.create_search_function(@ions, :prob_min_bincnt => 20)
    end

    xit 'searches mz values' do
      searcher = Mspire::Lipid::Search.new(@ions, :query_min_count_per_bin => 8, :num_rand_samples_per_bin => 1000, :ppm => false)
      num_nearest_hits = 3
      (hit_groups, qvals) = searcher.search(@pretend_search_mzs, 3)
      p hit_groups.map(&:first).map(&:pvalue)
    end
  end

  describe 'searching a full lipid maps' do

    before do
      # this will be specific to your install since it's not part of install
      path_to_lipidmaps_db = "#{ENV['HOME']}/tmp/tamil/lipidmaps_20120103_classes_1_2_3_4_5_6_7_8.exact_mass.tsv"
      @lipids = Mspire::LipidMaps.parse_file(path_to_lipidmaps_db)
      @ions = @lipids.map do |lipid| 
        [[@proton], [@proton, @proton], [@proton, @h2o_loss]].map do |mods|
          Mspire::Lipid::Search::Query.new(lipid, mods)
        end
      end.flatten(1)
      @pretend_search_mzs = [187.157, 396.20, 244.30, 618.22, 933.01]
    end

    xit 'returns hit groups parallel with input m/zs' do
      searcher = Mspire::Lipid::Search.new(@ions, :query_min_count_per_bin => 1000, :ppm => false)
      hit_groups = searcher.search(@pretend_search_mzs, 3)
      best_hits = hit_groups.map(&:best_hit)
      best_hits.map {|hit| hit.observed_mz }.should == @pretend_search_mzs
    end

    it 'works with :ppm => true'

  end

end
