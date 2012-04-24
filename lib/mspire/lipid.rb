
module Mspire
  class Lipid
    def self.members
      [:lm_id,:common_name,:systematic_name,:formula,:mass,:category,:main_class,:sub_class,:pubchem_id,:inchi_key,:kegg_id,:chebi_id,:structure]
    end

    members.each {|mem| attr_accessor mem }

    def initialize(*args)
      (@lm_id,@common_name,@systematic_name,@formula,@mass,@category,@main_class,@sub_class,@pubchem_sid, @inchi_key, @kegg_id, @chebi_id, @structure) = args
    end

    def inspect
      cut_common_name = (common_name.size <= 20) ? common_name : (common_name[0,20]+"...")
      "<#{lm_id}: #{formula}: #{mass} #{cut_common_name}>"
    end
  end
end
