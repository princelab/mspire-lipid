
module MS
  class Lipid
    class Search
      class Hit
        # an MS::Lipid::Search::Query object
        attr_accessor :query
        # the experimental m/z value
        attr_accessor :observed_mz
        # the probability the hit is due to random chance
        attr_accessor :pvalue

        def initialize(hash={})
          hash.each {|k,v| instance_variable_set("@#{k}", v) }
        end

        # observed_mz - query m/z
        def delta
          qmz = @query.mz.to_f
          (@observed_mz - qmz) / qmz
        end

        # parts per million (divided by theoretical m/z)
        def ppm
          (delta / @query.mz) * 1e6
        end

        def theoretical_mz
          @query.mz
        end
      end
    end
  end
end


