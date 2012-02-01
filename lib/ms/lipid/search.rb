
module MS
  module Lipid
    class Search
      STANDARD_MODIFICATIONS = {
        :proton => [1,2],
        :ammonium => [1],
        :lithium => [1],
        :water => [1,2],
      }
      STANDARD_SEARCH = {
        :tolerance => 5,
        :units => :ppm,
        :start_mz => 300,
        :end_mz => 2000,
        :prob_binsize => 100,
        :prob_bincnt => 100,
      }

      attr_accessor :options

      # will generate PossibleLipid objects and return a new search object
      # uses only one kind of loss at a time and one type of gain at a time
      # will also do the combination of a gain and a loss if gain_and_loss is
      # true
      def self.simple_search(lipids, mods=STANDARD_MODIFICATIONS, gain_and_loss=false)
        possible_lipids = []
        real_mods_and_cnts = mods.map {|name, cnts| [MS::Lipid::Modification.new(name), cnts] }
        # one of each
        real_mods_and_cnts.each do |mod, counts|
          counts.each do |cnt|
            possible_lipids << MS::Lipid::Search::PossibleLipid.new(lipid, Array.new(cnt, mod))
          end
        end
        if gain_and_loss
          # one of each gain + one of each loss
          (gain_mod_cnt_pairs, loss_mod_cnt_pairs) = real_mods_and_cnts.partition {|mod, count| mod.gain }
          gain_mod_cnt_pairs.each do |mod, cnt|
            lipids.each do |lipid|
              #### need to implement still (use combinations or something...)
              get_this_working!
            end
          end
        end
        self.new(possible_lipids)
      end

      # lipids is an array of Lipid objects.
      # block yields(lipids, modifications
      def initialize(possible_lipids)
        create_probability_function(possible_lipids)
      end

      def search(mzs, opts={})
      end

      def create_probability_function(possible_lipids, opts)
        rng = Random.new
        mzs = possible_lipids.map(&:mz)
        # create a mock spectrum of intensity 1
        spec = Spectrum.new(mzs, Array.new(mzs.size,1))
        mzs.each_slice(opts[:prob_bincnt]) do |mzs|
          random_mzs = opts[:prob_bincnt].times.map { rng.rand(mzs.first mzs.last) }
          deltas = random_mzs.map {|random_mz| (random_mz - spec.find_all_nearest(random_mz).first).abs }
          File.write("data.txt", deltas.join("\n"))
        end
      end
    end
  end
end




