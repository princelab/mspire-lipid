
module MS
  class Lipid
    class Search
      class ProbabilityDistribution
        DEFAULT_TYPE = :ppm
        R = Rserve::Simpler.new
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
          R.converse "pgev(log(#{dev}), #{@location}, #{@scale}, #{@shape})"
        end

        def self.require_r_library(lib)
          reply = R.converse "library(#{lib})"
          unless reply.size > 4  # ~roughly
            $stderr.puts "The libraries ismev and evd must be installed in your R env!"
            $stderr.puts "From within R (works best if R is started with sudo or root for installing):"
            $stderr.puts %Q{install.packages("ismev") ; install.packages("evd")}
            raise "must have R (rserve) and ismev and evd installed!"
          end
        end

        # returns an EVD object
        def self.deviations_to_probability_distribution(type, devs)
          %w(ismev evd).each {|lib| require_r_library(lib) }
          params = R.converse("m <- gev.fit(log(devs_r))\n c(m$mle[1], m$mle[2], m$mle[3])", :devs_r => devs )
          self.new(*params, type)
        end
      end
    end
  end
end
