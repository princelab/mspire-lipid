#!/usr/bin/env ruby

require 'ms/mzml'
require 'ms/lipid/search'
require 'ms/lipid/search/query'
require 'ms/lipid_maps'

class Sample
  attr_accessor :file
  attr_accessor :spectrum
  def initialize(file)
    @file = file
    @spectrum = merge_ms1_spectra(file)
  end

  # returns a single spectrum object
  def merge_ms1_spectra(file)
    spectra = []
    MS::Mzml.foreach(file) do |spectrum|
      spectra << spectrum if spectrum.mzs.size > 1000  # <<<<<<------ kludge for ms_level == 1
    end
    MS::Spectrum.merge(spectra, :bin_width => 3, :bin_unit => :ppm)
  end
end

if ARGV.size < 2
  puts "usage: #{File.basename(__FILE__)} lipidmaps.tsv <file>.mzML ..."
  exit
end

(lipidmaps, *files) = ARGV

$VERBOSE = 5

proton = MS::Lipid::Modification.new(:proton)
h2o_loss = MS::Lipid::Modification.new(:water, :loss => true)

lipids = MS::LipidMaps.parse_file(lipidmaps)

queries = lipids.map do |lipid| 
  [[proton], [proton, h2o_loss]].map do |mods|
    MS::Lipid::Search::Query.new(lipid, mods)
  end
end.flatten(1)

searcher = MS::Lipid::Search.new(queries, :ppm => false)

files.each do |file|
  puts "/\\"*80
  puts "FILE: #{file}"
  sample = Sample.new(file)

  hit_groups = searcher.search(sample.spectrum)

  hit_groups.map(&:first).sort_by(&:pvalue).each do |hit|
    puts "="*80
    p hit
  end
end
