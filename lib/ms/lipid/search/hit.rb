
module MS
  class Lipid
    class Search
      class Hit
        # an array of MS::Lipid::Search::Query objects
        attr_accessor :query_group
        # the experimental m/z value
        attr_accessor :observed_mz
        # the probability the hit is due to random chance
        attr_accessor :pvalue
        # the FDR if the threshold accepts this pvalue.  Note that this value
        # is relative to the number of tests performed and not completely 
        # intrinsic to the hit itself.
        attr_accessor :qvalue

        def initialize(hash={})
          hash.each {|k,v| instance_variable_set("@#{k}", v) }
        end

        # observed_mz - query m/z
        def delta
          puts @observed_mz
          @observed_mz - @query_group.first.mz.to_f
        end

        # parts per million (divided by theoretical m/z)
        def ppm
          (delta / @query_group.first.mz) * 1e6
        end

        def theoretical_mz
          @query_group.first.mz
        end

        def inspect
          "<<#{super} -- <ppm=#{ppm} delta=#{delta} theoretical_mz=#{theoretical_mz}>>"
        end
      end

      # A query that matched multiple items.  Each search returns a hit group
      # which consists of the best hits for that experimental m/z.  When
      # queried for values like delta or ppm, it will delegate to the first hit.
      # So, in many ways it can be used as a container for hits, but it puts
      # its best face forward.
      class HitGroup < Array

        # should implement with delegator obviously...
        # should allow setting ???

        def delta() first.delta end
        def ppm() first.ppm end
        def theoretical_mz() first.theoretical_mz end
        def query_group() first.query_group end
        def observed_mz() first.observed_mz end
        def pvalue() first.pvalue end
        def qvalue() first.pvalue end

        def best_hit() first end

      end
    end
  end
end
