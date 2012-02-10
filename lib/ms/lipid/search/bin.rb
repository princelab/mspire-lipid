
module MS
  class Lipid
    class Search

      # A Search::Bin is a range that contains the *entire* query spectrum
      # (not just the portion covered by the range).  the query spectrum, and
      # a ProbabilityDistribution -- the probability that a peak's delta to
      # nearest peak is that small by chance.
      class Bin < Range
        # the intensity value of the query spectrum should be a query
        attr_accessor :query_spectrum
        attr_accessor :probability_distribution

        def initialize(range_obj, query_spectrum)
          super(range_obj.begin, range_obj.end, range_obj.exclude_end?)
          @query_spectrum = query_spectrum
        end

        # searches many mzs and returns an array of HitGroup objects parallel
        # to mzs
        def hit_groups(mzs, hit_group_size=1)
              
        end

        # returns an array of MS::Lipid::Search::HitGroups sorted by delta
        # [with tie going to the lower m/z]
        def hit_group(mz, hit_group_size=1)
          # this could be improved by updating the probability_distribution if
          # it happens to jump to the next major Search::Bin of course, this
          # has the advantage that the next nearest peaks will have the same
          # probability distribution as the nearest peak.

          hits = delta_index_pairs.map do |delta, index|
            hit = Hit.new( :query_group => query_groups[index], :observed_mz => mz)
            dev = (probability_distribution.type == :ppm) ? hit.ppm : hit.delta.abs
            hit.pvalue = probability_distribution.pvalue(dev)
            hit
          end
          MS::Lipid::Search::HitGroup.new(hits)
        end

        # returns an array of doublets [delta, query_group] the nearest doublets to that m/z
        # and the associated QueryGroup object
        def nearest_points(mz, num_query_groups)
          mzs = @query_spectrum.mzs
          query_groups = @query_spectrum.intensities
          index = @query_spectrum.find_nearest_index(mz)
          _min = index - (hit_group_size-1)
          (_min >= 0) || (_min = 0)
          _max = index + (hit_group_size-1)
          (_max < mzs.size) || (_max = @query_spectrum - 1)
          delta_index_pairs = (_min.._max).map {|i| [mz.-(mzs[i]).abs, i] }.sort[0, hit_group_size]
        end

        def inspect
          "<(#{super}) @query_spectrum(points size)=#{query_spectrum.mzs.size} @probability_distribution=#{probability_distribution}>"
        end

        def to_range
          Range.new( self.begin, self.end, self.exclude_end? )
        end
      end
    end
  end
end


