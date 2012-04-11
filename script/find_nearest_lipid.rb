#!/usr/bin/env ruby

require 'trollop'

require 'ms/lipid/search'
require 'ms/lipid/ion'
require 'ms/lipid_maps'

parser = Trollop::Parser.new do
  banner "usage: #{File.basename(__FILE__)} <lipidmaps.tsv> <m/z> ..."
  opt :top_n, "how many closest ions to print", :default => 3
  text "modifications (at least 1 charged mod is required):"
  opt :lithium, "search for lithium adducts"
  opt :ammonium, "search for ammonium adducts"
  opt :proton_gain, "search for proton gain"
  opt :proton_loss, "search for proton loss"
  opt :water_loss, "*all* mods are also considered with water loss"
end

opts = parser.parse(ARGV)

if ARGV.size < 2
  parser.educate
  exit
end

CHARGED_MODS = [:lithium, :ammonium, :proton_gain, :proton_loss]

unless CHARGED_MODS.any? {|key| opts[key] }
  puts "*" * 78
  puts "ArgumentError: need at least one charged mod!"
  puts "*" * 78
  parser.educate
  exit
end

(lipidmaps, *actual_mzs) = ARGV

$VERBOSE = opts[:verbose]

MSLM = MS::Lipid::Modification

mods = {
  proton_gain: MSLM.new(:proton),
  water_loss: MSLM.new(:water, :loss => true),
  lithium: MSLM.new(:lithium),
  ammonium: MSLM.new(:ammonium),
  proton_loss: MS::Lipid::Modification.new(:proton, :loss => true, :charge => -1)
}

lipids = MS::LipidMaps.parse_file(lipidmaps)

ions = []
lipids.each do |lipid| 
  CHARGED_MODS.each do |key|
    if opts[key]
      ions << MS::Lipid::Ion.new(lipid, [mods[key]])
      if opts[:water_loss]
        ions << MS::Lipid::Ion.new(lipid, [mods[key], mods[:water_loss]]) 
      end
    end
  end
end

actual_mzs.map(&:to_f).each do |actual_mz|
  puts "*" * 70
  puts "ACTUAL M/Z: #{actual_mz}"
  closest = ions.sort_by {|ion| (ion.mz - actual_mz).abs }
  closest[0,opts[:top_n]].each do |ion|
    p ion
  end
end
