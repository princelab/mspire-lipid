


module Mspire
  class Lipid
    class Search
      class Query

        # the experimentally observed lowest mz
        attr_accessor :mz

        # the index of search spectrum that the m/z was derived from
        # this allows for the creation of an isotope envelope starting from a
        # particular m/z value.
        attr_accessor :index

        def initialize(mz, index)
          @mz, @index = mz, index
        end
      end
    end
  end
end
