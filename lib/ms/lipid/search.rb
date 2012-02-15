require 'ms/spectrum'
require 'rserve/simpler'  # TODO: move to integrated interface with rserve when available
require 'core_ext/array/in_groups'
require 'ms/lipid/search/hit'
require 'ms/lipid/search/bin'
require 'ms/lipid/modification'
require 'ms/lipid/search/probability_distribution'

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
        :query_min_count_per_bin => 500,  # min number of peaks per bin
        :num_rand_samples_per_bin => 1000,
      }

      attr_accessor :options
      attr_accessor :search_function

      # will generate PossibleLipid objects and return a new search object
      # uses only one kind of loss at a time and one type of gain at a time
      # will also do the combination of a gain and a loss if gain_and_loss is
      # true
      def self.generate_simple_queries(lipids, mods=STANDARD_MODIFICATIONS, gain_and_loss=false)
        possible_lipids = []
        real_mods_and_cnts = mods.map {|name, cnts| [MS::Lipid::Modification.new(name), cnts] }
        # one of each
        real_mods_and_cnts.each do |mod, counts|
          counts.each do |cnt|
            possible_lipids << MS::Lipid::Search::Query.new(lipid, Array.new(cnt, mod))
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

      # ions are MS::Lipid::Ion objects
      # each one should give a non-nil m/z value
      def initialize(ions=[], opts={})
        @options = opts
        @db_isobar_spectrum = create_db_isobar_spectrum(ions)
        @search_function = create_search_function(ions, opts)
      end

      # returns an array of HitGroup and a parallel array of BH derived
      # q-values (will switch to Storey soon enough).  The HitGroups are
      # returned in the order in which the mz_values are given.
      def search(mz_values, num_nearest=3)

        # associate each mz_value with a bin


        # takes an array rather than single m/z values because we can
        # eventually do a ratchet search and it allows us to calculate FDR for
        # the whole group.
        hit_groups = mz_values.map do |mz|
          @search_function.call(mz, num_nearest)
        end

        # from http://stats.stackexchange.com/questions/870/multiple-hypothesis-testing-correction-with-benjamini-hochberg-p-values-or-q-va
        # but I've already coded this up before, too, in multiple ways...
        prev_bh_value = 0
        num_total_tests = hit_groups.size

        #hit_groups.each {|hg| p [hg.first.pvalue, hg] }

        # calculate Q-values BH style for now:
        # first hit is the best hit in the group
        pval_hg_index_triplets = hit_groups.each_with_index.map {|hg,i| [hg.first.pvalue, hg, i] }

        if pval_hg_index_triplets.any? {|pair| pair.first.nan? }
          $stderr.puts "pvalue of NaN!"
          $stderr.puts ">>> Consider increasing query_min_count_per_bin or setting ppm to false <<<"
          raise
        end

        pval_hg_index_triplets.sort.each_with_index do |triplet,i|
          pval = triplet.first
          bh_value = pval * num_total_tests / (i + 1)
          # Sometimes this correction can give values greater than 1,
          # so we set those values at 1
          bh_value = [bh_value, 1].min

          # To preserve monotonicity in the values, we take the
          # maximum of the previous value or this one, so that we
          # don't yield a value less than the previous.
          bh_value = [bh_value, prev_bh_value].max
          prev_bh_value = bh_value
          triplet[1].first.qvalue = bh_value # give the top hit the q-value
        end

        # return hit groups
        pval_hg_index_triplets.sort_by(&:last).map {|trip| trip[1] }
      end


      def create_search_function(ions, opts={})
        opt = STANDARD_SEARCH.merge( opts )

        db_isobar_spectrum = create_db_isobar_spectrum(ions)

        search_bins = create_search_bins(db_isobar_spectrum, opt[:query_min_count_per_bin])

        create_probability_distribution_for_search_bins!(search_bins, db_isobar_spectrum, opt[:num_rand_samples_per_bin], opt[:ppm])
        search_bins
        # create the actual search function
        # always returns a hit group
        lambda do |mz, num_nearest_hits|
          bin = search_bins.find {|bin| bin === mz } 
          # returns a HitGroup
          bin.nearest_hits(mz, num_nearest_hits)
        end
      end

      #####################################################
      # Ancillary to create_search_function:
      #####################################################

      # returns a DB isobar spectrum where the m/z values are all the m/z
      # values to search for and the intensities each an array corresponding
      # to all the lipid ions matching that m/z value
      def create_db_isobar_spectrum(ions)
        mzs = [] ; query_groups = []
        pairs = ions.group_by(&:mz).sort_by(&:first)
        pairs.each {|mz, ar| mzs << mz ; query_groups << ar }
        MS::Spectrum.new([mzs, query_groups])
      end

      # use_ppm uses ppm or amu if false
      # returns the search_bins
      def create_probability_distribution_for_search_bins!(search_bins, db_isobar_spectrum, num_rand_samples_per_bin, use_ppm=true)
        search_bins.each do |search_bin| 
          rng = Random.new
          random_mzs = num_rand_samples_per_bin.times.map { rng.rand(search_bin.to_range)  }
          # find the deltas
          diffs = random_mzs.map do |random_mz| 
            nearest_random_mz = db_isobar_spectrum.find_nearest(random_mz)
            delta = (random_mz - nearest_random_mz).abs
            use_ppm ? delta./(nearest_random_mz).*(1e6) : delta
          end
          search_bin.probability_distribution = ProbabilityDistribution.deviations_to_probability_distribution((use_ppm ? :ppm : :amu), diffs)
        end
        search_bins
      end

      def create_search_bins(db_isobar_spectrum, min_n_per_bin)

        # make sure we get the right bin size based on the input
        ss = db_isobar_spectrum.mzs.size ; optimal_num_groups = 1
        (1..ss).each do |divisions|
          if  (ss.to_f / divisions) >= min_n_per_bin
            optimal_num_groups = divisions
          else ; break
          end
        end

        mz_ranges = []
        prev = nil

        groups = db_isobar_spectrum.points.in_groups(optimal_num_groups,false).to_a

        case groups.size
        when 0
          raise 'I think you need some data in your query spectrum!'
        when 1
          group = groups.first
          [ MS::Lipid::Search::Bin.new( Range.new(group.first.first, group.last.first), db_isobar_spectrum ) ]
        else
          search_bins = groups.each_cons(2).map do |points1, points2|
            bin = MS::Lipid::Search::Bin.new( Range.new(points1.first.first, points2.first.first, true), db_isobar_spectrum )
            prev = points2
            bin
          end
          _range = Range.new(prev.first.first, prev.last.first)
          search_bins << MS::Lipid::Search::Bin.new(_range, db_isobar_spectrum) # inclusive
        end
      end
    end
  end
end




