
require 'ms/mass'

module MS
  # the mass is the formula + the mass diffs of the modifications
  Lipid = Struct.new(:lm_id,:common_name,:systematic_name,:formula,:mass,:category,:main_class,:sub_class)

  class Lipid
    class Modification
      # mass of electron: 5.486 x 10-4 amu
      FORMULAS = {
        :proton => 'H+',
        :ammonium => 'NH3H+',
        :lithium => 'Li',
        :water => 'H2O',
      }
      # determined by running formulas through MS::Mass.formula_to_mass
      MASSDIFFS = {}
      FORMULAS.each do |name, formula|
        MASSDIFFS[name] = MS::Mass.formula_to_mass(FORMULAS[name])
      end
      CHARGE = {
        :proton => 1,
        :ammonium => 1,
        :lithium => 1,
        :water => 0,
      }
      GAIN = {
        :proton => true, 
        :ammonium => true,
        :lithium => true,
        :water => false,
      }
      attr_accessor :name
      attr_accessor :formula
      attr_accessor :mass
      # the charge 
      attr_accessor :charge
      # if no mass or formula is given then it searches command mods for the name
      # @param [Symbol] name the name of the mod
      # A number of opts are allowed:
      #
      #     :formula = the chemical formula, lipidmaps style ("C2H4BrO")
      #     :mass = Float
      #     :charge = +/- Integer
      #     :gain = boolean # is this a mass lost or gain
      def initialize(name, opts={})
        @name = name
        @formula = opts[:formula] || FORMULAS[name]
        @massdiff = opts[:massdiff] || MASSDIFFS[name]
        @charge = opts[:charge] || CHARGE[name]
        @gain = opts[:gain] || GAIN[name]
      end
    end

    class Search
      # a Lipid Search
      class PossibleLipid
        # an MS::Lipid object
        attr_accessor :lipid
        # an MS::Lipid::Modifications object
        attr_accessor :modifications
        # the key attribute of a query

        def mz
          return @mz if @mz
          mass = @lipid.mass
          charge = 0
          @modifications.each do |mod|
            mass = mod.gain ? mass + mod.mass : mass - mod.mass
            charge += mod.charge
          end
          @mz = mass / charge
        end
      end

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
        #:prob_binsize => 100,
        :prob_bincnt => 100,
        # if flexible_bin_size then overides prob_binsize
        :flexible_bin_size => true, 
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


  module LipidMaps
    # returns an array of Lipids
    # if high_res_mass is true (default), then the formula is used to calculate a higher
    # resolution mass than what is in lipidmaps
    def self.parse_file(lipidmaps_tsv, high_res_mass=true, skip_clas_defs=true)
      first_line = nil
      IO.foreach(lipidmaps_tsv).map do |line|
        line.chomp!
        pieces = line.split("\t")
        if pieces[3] !~ /[A-Z]/  # <- there is no formula!
          nil
        else
          pieces[4] = MS::Mass.formula_to_exact_mass(pieces[3]) if high_res_mass
          if first_line
            Lipid.new *pieces
          else
            first_line = pieces
            warn "lipidmaps column headers are not right!" unless first_line.map(&:downcase) == Lipid.members.map(&:to_s)
            nil
          end
        end
      end.compact
    end
  end

end


=begin
require 'bsearch'
ar = [1.0, 2.0, 3.0, 4.0, 5.0]

delta = 0.6
val = 2.5
reply = ar.bsearch_range do |x|
  if (val - x).abs <= delta
    0
  else
    x <=> val
  end
end
=end
