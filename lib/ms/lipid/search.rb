require 'ms/spectrum'
require 'rserve/simpler'  # TODO: move to integrated interface with rserve when available
require 'core_ext/array/in_groups'

module MS
  class Lipid
    class Search
      class EVD
        EVD_R = Rserve::Simpler.new
        attr_accessor :location, :scale, :shape
        # takes location, scale and shape parameters
        def initialize(location, scale, shape)
          @location, @scale, @shape = location, scale, shape
        end
        # takes a ppm and returns the pvalue
        # note that it will take the 
        def pvalue(ppm)
          EVD_R.converse "pgev(log(#{ppm}), #{@location}, #{@scale}, #{@shape})"
        end

        # returns an EVD object
        def self.ppms_to_evd(ppms)
          %w(ismev evd).each do |lib|
            reply = EVD_R.converse "library(#{lib})"
            puts "REAPLY: #{reply}"
            unless reply.size > 4  # ~roughly
              $stderr.puts "The libraries ismev and evd must be installed in your R env!"
              $stderr.puts "From within R:"
              $stderr.puts %Q{install.packages("ismev") ; install.packages("evd")}
              raise "must have R (rserve) and ismev and evd installed!"
            end
          end
          #ppmsl = ppms.map {|v| Math::log(v) }
          #p EVD_R.converse("cor(c(1,2,3), c(4,4,6))", :ppmsl_r => ppmsl )
          #params = EVD_R.converse("hist(ppms_r, 50)", :ppms_r => ppms )
          params = EVD_R.converse("m <- gev.fit(log(ppms_r)); c(m$mle[1], m$mle[2], m$mle[3])", :ppms_r => ppms )
          #EVD_R.command("svggev.diag(m)")
          #EVD_R.pause
          EVD.new(*params)
        end
            #ppmsl = log(ppms)
            #gev.fit(ppmsl)
        # pgev(log10(#{TEST_PVAL}), model$mle[1], model$mle[2], model$mle[3])

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

        spec = create_search_spectrum(queries)
        # make sure we get the right bin size based on the input
        ss = spec.mzs.size ; optimal_num_groups = 1
        (1..ss).each do |divisions|
          if  (ss.to_f / divisions) >= opts[:prob_min_bincnt]
            optimal_num_groups = divisions
          else ; break
          end
        end

        rng = Random.new
        # create a mock spectrum of intensity 1
        mz_ranges = []
        evds = spec.peaks.in_groups(optimal_num_groups,false).map do |peaks|
          mz_range = Range.new( peaks.first.first, peaks.last.first )
          mz_ranges << mz_range
          random_mzs = opts[:num_rand_samples_per_bin].times.map { rng.rand(mz_range) }
          puts "*****************************"
          puts "MZRANGE: "
          p mz_range
          # find the deltas
          diffs = random_mzs.map do |random_mz| 
            nearest_random_mz = spec.find_all_nearest(random_mz).first
            delta = (random_mz - nearest_random_mz).abs
            opts[:ppm] ? delta./(nearest_random_mz).*(1e6) : delta
          end
          #File.write("data_#{mz_range}.dataset", "deltas\n" + diffs.join("\n"))
          EVD.ppms_to_evd(diffs)
        end
        File.open("plot_evd_params_ppm_#{!!opt[:ppm]}.data",'w') do |out|
          out.puts %w(mzrange loc scale shape).join("\t")
          evds.zip(mz_ranges) do |evd, range|
            out.puts [range, evd.location, evd.scale, evd.shape].join("\t")
          end
        end
      end

    end
  end
end




