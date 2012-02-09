
module MS
  class Lipid
    class Search

      # A Search::Bin is a range that contains the *entire* query spectrum
      # (not just the portion covered by the range).  the query spectrum, and
      # an EVD that describes the probability that a peak's delta to nearest
      # peak is that small by chance.
      class Bin < Range
        # the intensity value of the query spectrum should be a query
        attr_accessor :db_spectrum
        attr_accessor :probability_distribution
        attr_accessor :queries

        def initialize(range_obj, db_spectrum)
          super(range_obj.begin, range_obj.end, range_obj.exclude_end?)
          @db_spectrum = db_spectrum
          @queries = []
        end

        def <<(query)

        # returns the nearest num_hits MS::Lipid::Search::Hits sorted by delta
        # [with tie going to the lower m/z]
        # searches all queries and removes them from the queries queue
        def queries_to_nearest_hits!(num_hits=1)
          query_mzs = queries.map {|query| query.mz }

          query_mzs.each do |query_mz|

            db_mzs = @db_spectrum.mzs
            query_groups = @db_spectrum.intensities
            index = @db_spectrum.find_nearest_index(query_mz)
            _min = index - (num_hits-1)
            (_min >= 0) || (_min = 0)
            _max = index + (num_hits-1)
            (_max < db_mzs.size) || (_max = @db_spectrum - 1)
            delta_index_pairs = (_min.._max).map {|i| [query_mz.-(db_mzs[i]).abs, i] }.sort[0, num_hits]
            # this could be improved by updating the probability_distribution if
            # it happens to jump to the next major Search::Bin of course, this
            # has the advantage that the next nearest peaks will have the same
            # probability distribution as the nearest peak.

            hits = delta_index_pairs.map do |delta, index|
              Hit.new( :db_isobar_groups => db_isobar_groups[index], :observed_mz => mz)
            end
            #dev = (probability_distribution.type == :ppm) ? hit.ppm : hit.delta.abs
            #hit.pvalue = probability_distribution.pvalue(dev)

          # updates the pvalues for all the hits
          probability_distribution.pvalues!( hits )

          MS::Lipid::Search::HitGroup.new(hits)
        end

        def inspect
          "<(#{super}) @db_spectrum(points size)=#{db_spectrum.mzs.size} @probability_distribution=#{probability_distribution}>"
        end

        def to_range
          Range.new( self.begin, self.end, self.exclude_end? )
        end
      end
    end
  end
end


