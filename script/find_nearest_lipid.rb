#!/usr/bin/env ruby

require 'trollop'

require 'bsearch'

require 'mspire/lipid/ion'
require 'mspire/lipid/modification'
require 'mspire/lipid_maps'

parser = Trollop::Parser.new do
  banner "usage: #{File.basename(__FILE__)} <lipidmaps.tsv> <m/z> ..."
  banner ""
  opt :top_n, "how many closest ions to print", :default => 3
  banner ""
  text "modifications: (at least 1 charged mod is required)"
  opt :lithium, "search for i down to 1 lithium adducts", :default => 0
  opt :sodium, "search for i down to 1 sodium adducts", :default => 0
  opt :ammonium, "search for i down to 1 ammonium adducts", :default => 0
  opt :proton_gain, "search for i down to 1 proton additions", :default => 0
  opt :proton_loss, "search for i down to 1 proton losses", :default => 0
  opt :water_loss, "if used, *all* mods are also considered with i down to 0 water losses", :default => 0
  opt :textfile, "a text file with m/z values, one per line", :type => String
  opt :lower_bound, "use lower bound searching (requires m/z's to be in sorted order)"
  opt :sort, "sorts the m/z's"
end

opts = parser.parse(ARGV)

if ARGV.size == 0
  parser.educate
  exit
end

CHARGED_MODS = [:lithium, :sodium, :ammonium, :proton_gain, :proton_loss]

unless CHARGED_MODS.any? {|key| opts[key] > 0 }
  puts "*" * 78
  puts "ArgumentError: need at least one charged mod!"
  puts "*" * 78
  parser.educate
  exit
end

(lipidmaps, *actual_mzs) = ARGV

if f=opts[:textfile]
  actual_mzs = IO.readlines(f).map(&:chomp)
end

actual_mzs.map!(&:to_f)

unless actual_mzs.size > 0
  STDERR.puts "NO m/z values given!!!"
  parser.educate
  exit
end

$VERBOSE = opts[:verbose]

LipidMod = Mspire::Lipid::Modification

mods = {
  proton_gain: LipidMod.new(:proton),
  water_loss: LipidMod.new(:water, :loss => true),
  lithium: LipidMod.new(:lithium),
  sodium: LipidMod.new(:sodium),
  ammonium: LipidMod.new(:ammonium),
  proton_loss: LipidMod.new(:proton, :loss => true)
}

lipids = Mspire::LipidMaps.parse_file(lipidmaps)

ions = []
lipids.each do |lipid| 
  CHARGED_MODS.each do |key|
    if opts[key] > 0
      opts[key].downto(1) do |num_charge_mod|
        mods_to_use = [mods[key]] * num_charge_mod
        opts[:water_loss].downto(0) do |i|
          ions << Mspire::Lipid::Ion.new(lipid, mods_to_use + ([mods[:water_loss]]*i)) 
        end
      end
    end
  end
end

ions.sort_by!(&:mz)

PPM = 1500

cats = %w(exp_mz rank ppm ppm_abs category lm_id common_name modifications)
puts cats.join("\t")

actual_mzs.sort! if opts[:sort]

starting_i = 0
len = ions.size

actual_mzs.each do |exp_mz|

  range = ions.bsearch_range(starting_i...len) do |ion| 
    part_diff = ((exp_mz - ion.mz)/ion.mz)*1e6
    if part_diff.abs <= PPM then 0
    elsif part_diff < PPM then 1
    else 
      -1
    end
  end

  starting_i = range.begin if opts[:lower_bound]
  
  closest = ions[range].sort_by {|ion| [(ion.mz - exp_mz).abs, ion.mz] }
  row = [exp_mz]
  closest[0,opts[:top_n]].each_with_index do |ion,i|
    rank = i + 1
    ppm = ((exp_mz - ion.mz) / ion.mz) * 1e6
    lipid = ion.lipid
    row.push( rank, ppm, ppm.abs, lipid.category, lipid.lm_id, lipid.common_name, ion.modifications.map(&:charged_formula_string).join(", ") )
  end
  puts row.join("\t")
end
