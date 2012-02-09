#!/usr/bin/env ruby

require 'ms/mzml'
require 'ms/lipid/search'
require 'ms/lipid/search/query'
require 'ms/lipid_maps'

if ARGV.size != 2
  puts "need: lipimaps mzML "
  exit
end

(lipidmaps, file) = ARGV

spectra = []
MS::Mzml.foreach(file) do |spectrum|
  spectra << spectrum if spectrum.mzs.size > 700  # <<<<<<------ kludge for ms_level == 1
end

$VERBOSE = 5

master_spectrum = MS::Spectrum.merge(spectra, :bin_width => 1, :bin_unit => :ppm)

puts "LOOKING AT MASTER:"
p master_spectrum.mzs.size
p master_spectrum.mzs
p master_spectrum.mzs.size
p master_spectrum.intensities
abort 'DONE MAKING MASTER SPECTRUM'

proton = MS::Lipid::Modification.new(:proton)
h2o_loss = MS::Lipid::Modification.new(:water, :loss => true)

lipids = MS::LipidMaps.parse_file(lipidmaps)

queries = lipids.map do |lipid| 
  [[proton], [proton, h2o_loss]].map do |mods|
    MS::Lipid::Search::Query.new(lipid, mods)
  end
end.flatten(1)

searcher = MS::Lipid::Search.new(queries, :ppm => false)

hit_groups = searcher.search(spectra[3].mzs[200,20])

hit_groups.map(&:first).sort_by(&:pvalue).each do |hit|
  puts "="*80
  p hit
end
