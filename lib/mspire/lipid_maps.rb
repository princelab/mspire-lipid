require 'mspire/lipid'
require 'mspire/mass'

module Mspire
  module LipidMaps
    # returns an array of Lipids
    # if high_res_mass is true (default), then the formula is used to calculate a higher
    # resolution mass than what is in lipidmaps
    def self.parse_file(lipidmaps_tsv, high_res_mass=true, skip_clas_defs=true)
      seen_first_line = false
      IO.foreach(lipidmaps_tsv).map do |line|
        line.chomp!
        pieces = line.split("\t")
        if pieces[3] !~ /[A-Z]/  # <- there is no formula!
          nil
        else
          if seen_first_line
            pieces[4] = Mspire::Mass.formula_to_exact_mass(pieces[3]) if high_res_mass
            l = Mspire::Lipid.new *pieces
          else
            seen_first_line = true
            warn "lipidmaps column headers are not right!" unless pieces.map(&:downcase) == Mspire::Lipid.members.map(&:to_s)
            nil
          end
        end
      end.compact
    end
  end
end


