require 'mspire/lipid'
require 'mspire/mass'
require 'mspire/molecule'

module Mspire
  module LipidMaps

    DEFAULTS = {
      :high_res_mass => true,
      :sdf_objects => false,
      :molecular_formula_objects => true,
    }

    # returns an array of Lipids
    # if high_res_mass is true (default), then the formula is used to calculate a higher
    # resolution mass than what is in lipidmaps
    #
    #     :high_res_mass => true (ensures that a high res mass is present or calculated)
    def self.parse_file(lipidmaps_tsv, opts={})

      opts = DEFAULTS.merge(opts)

      io = File.open(lipidmaps_tsv)
      header = io.readline.split("\t")
      # the lipidmaps_filetype
      lm_ft = case header.size
              when 8
                :programmatic
              when 20
                :download
              when 21
                :download_sd
              end
      index_mapping = 
        case lm_ft
        when :programmatic
          (0...(Mspire::Lipid.members.size)).to_a
        when :download, :download_sd
          indices = {
            :lm_id => 0,
            :systematic_name => 1,
            :category => 3,
            :main_class => 4,
            :mass => 5,
            :formula => 6,
            :pubchem_id => 7,
            :inchi_key => 8,
            :common_name => 11,
            :kegg_id => 12,
            :chebi_id => 13,
            :sub_class => 14,
            :structure => 20,
          }
          Mspire::Lipid.members.map {|key| indices[key] }
        end
      
      formula_i = index_mapping[Mspire::Lipid.members.index(:formula)]

      lipids = io.each_line.map do |line|
        line.chomp!
        data = line.split("\t")
        if data[formula_i] =~ /[A-Z]/  # <- there is a formula!
          lipid = Mspire::Lipid.new( *index_mapping.map {|i| data[i] } )
          lipid.mass = lipid.mass.to_f
          lipid
        end
      end.compact

      if opts.values_at(:molecular_formula_objects, :sdf_objects).any? || (opts[:high_res_mass] && lm_ft == :programmatic)
        lipids.each do |lipid|
          if opts[:molecular_formula_objects]
            lipid.formula = Mspire::MolecularFormula.new(lipid.formula)
          end
          if lm_ft == :programmatic && opts[:high_res_mass]
            lipid.mass = Mspire::Mass.formula_to_exact_mass(lipid.formula)
          end
          if opts[:sdf_objects]
            lipid.structure = Mspire::Molecule.from_lipidmaps_sdf_string(lipid.structure)
          end
        end
      end
      lipids
    end
  end
end


