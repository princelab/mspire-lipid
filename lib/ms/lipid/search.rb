require 'ms/spectrum'
require 'rserve/simpler'  # TODO: move to integrated interface with rserve when available
require 'core_ext/array/in_groups'
require 'ms/lipid/search/hit'

module MS
  class Lipid
    class Search

      # A Search::Bin is a range that contains the *entire* query spectrum
      # (not just the portion covered by the range).  the query spectrum, and
      # an EVD that describes the probability that a peak's delta to nearest
      # peak is that small by chance.
      class Bin < Range
        # the intensity value of the query spectrum should be a query
        attr_accessor :query_spectrum
        attr_accessor :evd

        def initialize(range_obj, query_spectrum)
          super(range_obj.begin, range_obj.end, range_obj.exclude_end?)
          @query_spectrum = query_spectrum
        end

        # returns the nearest num_hits MS::Lipid::Search::Hits sorted by delta
        # [with tie going to the lower m/z]
        def nearest_hits(mz, num_hits=1)
          mzs = @query_spectrum.mzs
          queries = @query_spectrum.intensities
          index = @query_spectrum.find_nearest_index(mz)
          _min = index - (num_nearest-1)
          (_min >= 0) || (_min = 0)
          _max = index + (num_nearest-1)
          (_max < @query_spectrum.size) || (_max = @query_spectrum - 1)
          delta_index_pairs = (_min.._max).map {|i| [mz.-(mzs[i]).abs, i] }.sort[0, num_to_return]
          # this could be improved by updating the evd if it happens to jump
          # to the next major Search::Bin
          # of course, this has the advantage that the next nearest peaks will
          # have the same probability distribution as the nearest peak.

          delta_index_pairs.map do |delta, index|
            hit = Hit.new( :query => queries[index], :observed_mz => mz)
            dev = (evd.type == :ppm) ? hit.ppm : hit.delta.abs
            hit.pvalue = evd.pvalue(dev)
            hit
          end
        end

        def to_range
          Range.new( self.begin, self.end, self.exclude_end? )
        end
      end

      class EVD
        DEFAULT_TYPE = :ppm
        EVD_R = Rserve::Simpler.new
        # takes location, scale and shape parameters
        attr_accessor :location, :scale, :shape
        # type is :ppm or :amu
        attr_accessor :type
        def initialize(location, scale, shape, type=DEFAULT_TYPE)
          @location, @scale, @shape = location, scale, shape
          @type = type
        end

        # takes a deviation and returns the pvalue
        def pvalue(dev)
          EVD_R.converse "pgev(log(#{dev}), #{@location}, #{@scale}, #{@shape})"
        end

        def self.require_r_library(lib)
          reply = EVD_R.converse "library(#{lib})"
          puts "REAPLY: #{reply}"
          unless reply.size > 4  # ~roughly
            $stderr.puts "The libraries ismev and evd must be installed in your R env!"
            $stderr.puts "From within R (works best if R is started with sudo or root for installing):"
            $stderr.puts %Q{install.packages("ismev") ; install.packages("evd")}
            raise "must have R (rserve) and ismev and evd installed!"
          end
        end

        # returns an EVD object
        def self.deviations_to_evd(type, devs)
          %w(ismev evd).each {|lib| require_r_library(lib) }
          params = EVD_R.converse("m <- gev.fit(log(devs_r))\n c(m$mle[1], m$mle[2], m$mle[3])", :devs_r => devs )
          EVD.new(*params, type)
        end
      end


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
        :prob_min_bincnt => 500,  # min number of peaks per bin (spread out over all)
        :num_rand_samples_per_bin => 1000,
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
        @prob_function = create_probability_function(queries) if queries.size > 0
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

        query_spectrum = create_search_spectrum(queries)

        # make sure we get the right bin size based on the input
        ss = query_spectrum.mzs.size ; optimal_num_groups = 1
        (1..ss).each do |divisions|
          if  (ss.to_f / divisions) >= opts[:prob_min_bincnt]
            optimal_num_groups = divisions
          else ; break
          end
        end

        mz_ranges = []
        last_peaks_group = nil
        search_bins = query_spectrum.peaks.in_groups(optimal_num_groups,false).each_cons(2).map do |peaks1, peaks2|
          bin = MS::Lipid::Search::Bin.new( Range.new(peaks1.first.first, peaks2.first.first, true), query_spectrum )
          last_peaks_group = peaks2
          bin
        end
        _range = Range.new(last_peaks_group.first.first, last_peaks_group.last.first)
        search_bins << MS::Lipid::Search::Bin.new(_range, query_spectrum) # inclusive

        search_bins.map do |search_bin| 
          rng = Random.new
          random_mzs = opts[:num_rand_samples_per_bin].times.map { rng.rand(search_bin.to_range)  }
          # find the deltas
          diffs = random_mzs.map do |random_mz| 
            nearest_random_mz = query_spectrum.find_nearest(random_mz)
            delta = (random_mz - nearest_random_mz).abs
            opts[:ppm] ? delta./(nearest_random_mz).*(1e6) : delta
          end
          EVD.deviations_to_evd(diffs, (opts[:ppm] ? :ppm : :amu))
        end
      end
    end
  end
end




