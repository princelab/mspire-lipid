require 'ms/spectrum'

module MS
  class Lipid
    class Search
      STANDARD_MODIFICATIONS = {
        :proton => [1,2],
        :ammonium => [1],
        :lithium => [1],
        :water => [1,2],
      }
      STANDARD_SEARCH = {
        :units => :ppm,
        :start_mz => 300,
        :end_mz => 2000,
        :prob_bincnt => 100,
        :num_rand_samples_per_bin => 100,
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

      # queries are MS::Lipid::Search::Query objects
      # each one should give a non-nil m/z value
      def initialize(queries=[])
        create_probability_function(queries) if queries.size > 0
      end

      # returns a funny kind of search spectrum where the m/z values are all
      # the m/z values to search for and the intensities each an array
      # corresponding to all the lipid queries matching that m/z value
      def create_search_spectrum(queries)
        mzs = [] ; query_groups = []
        pairs = queries.group_by(&:mz).sort_by(&:first)
        pairs.each {|mz, ar| mzs << mz ; query_groups << ar }
        MS::Spectrum.new([mzs, query_groups])
      end

      def create_probability_function(queries, opts={})
        opts = STANDARD_SEARCH.merge( opts )
        rng = Random.new
        spec = create_search_spectrum(queries)
        # create a mock spectrum of intensity 1
        spec.peaks.each_slice(opts[:prob_bincnt]) do |peaks|
          mz_range = Range.new( peaks.first.first, peaks.last.first )
          random_mzs = opts[:num_rand_samples_per_bin].times.map { rng.rand(mz_range) }
          deltas = random_mzs.map {|random_mz| (random_mz - spec.find_all_nearest(random_mz).first).abs }
        end
      end
    end
  end
end




