
module MS
  class Lipid
    def self.members
      [:lm_id,:common_name,:systematic_name,:formula,:mass,:category,:main_class,:sub_class]
    end

    members.each {|mem| attr_accessor mem }

    def initialize(*args)
      (@lm_id,@common_name,@systematic_name,@formula,@mass,@category,@main_class,@sub_class) = args
    end
  end
end
