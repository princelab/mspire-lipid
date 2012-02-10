
module MS
  class Lipid
    class Search
      # this is a group of lipids and modifications that all have the same m/z
      class DBIsobarGroup < Array
        attr_accessor :mz
        def initialize( ar=[], mz=nil)
          @mz = mz if mz
          self.replace(ar)
        end
      end
    end
  end
end
