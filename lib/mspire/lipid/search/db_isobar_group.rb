
module Mspire
  class Lipid
    class Search
      # this is a group of Lipid::Ion objects that all have the same (or
      # possibly similar) m/z
      class DBIsobarGroup < Array
        # it is implemented like this so that the isobar group *could* have
        # individuals in it with slightly different m/z values and this coudl
        # still be used as a container.  In my current implementation they
        # have exactly the same m/z
        attr_accessor :mz
        def initialize( ar=[], mz=nil)
          @mz = mz if mz
          self.replace(ar)
        end
      end
    end
  end
end
