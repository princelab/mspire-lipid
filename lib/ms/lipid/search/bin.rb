require 'bin'

module MS
  class Lipid
    class Search

      # A Search::Bin is a range that contains the *entire* query spectrum
      # (not just the portion covered by the range).  the query spectrum, and
      # a ProbabilityDistribution -- the probability that a peak's delta to
      # nearest peak is that small by chance.
      class Bin < ::Bin
        # the intensity value of the query spectrum should be a query
        attr_accessor :db_spectrum
        attr_accessor :probability_distribution

        def initialize(range_obj, db_spectrum)
          super(range_obj.begin, range_obj.end, range_obj.exclude_end?)
          @db_spectrum = db_spectrum
        end

        def <<(query)
          @data << query
        end

        # returns the nearest num_hits MS::Lipid::Search::Hits sorted by delta
        # [with tie going to the lower m/z]
        # searches all queries and removes them from the data queue
        def queries_to_hit_groups!(num_hits=1)
          queries = @data.dup
          @data.clear

          @db_isobar_groups_by_index = @db_spectrum.intensities

          hit_groups = queries.map do |query|
            best_hits(query, num_hits)
          end

          all_top_hits = hit_groups.map(&:first)

          # updates the pvalues for all the hits
          pvalues = probability_distribution.pvalues( all_top_hits )
          all_top_hits.zip(pvalues) {|hit, pvalue| hit.pvalue = pvalue }

          hit_groups
        end

        # returns a HitGroup object
        def best_hits(query, num_hits)
          query_mz = query.mz
          #puts "MZ: #{query_mz}"
          db_mzs = @db_spectrum.mzs
          index = @db_spectrum.find_nearest_index(query_mz)
          _min = index - (num_hits-1)
          (_min >= 0) || (_min = 0)
          _max = index + (num_hits-1)
          (_max < db_mzs.size) || (_max = @db_spectrum - 1)
          delta_index_pairs = (_min.._max).map {|i| [query_mz.-(db_mzs[i]).abs, i] }
          closest_delta_index_pairs = delta_index_pairs.sort
          top_num_hits_delta_index_pairs = closest_delta_index_pairs[0, num_hits]
          top_num_hit_indices = top_num_hits_delta_index_pairs.map(&:last)
          hit_group = top_num_hit_indices.map do |index|
            Hit.new( :db_isobar_group => @db_isobar_groups_by_index[index], :observed_mz => query_mz)
          end
          HitGroup.new(hit_group)
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


