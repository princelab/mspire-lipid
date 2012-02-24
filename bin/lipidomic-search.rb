#!/usr/bin/env ruby

require 'trollop'
require 'ms/mzml'
require 'ms/lipid/search'
require 'ms/lipid/ion'
require 'ms/lipid/search/query'
require 'ms/lipid_maps'

# for html output: (just make the id clickable)
LIPIDMAPS_SEARCH = "http://www.lipidmaps.org/data/LMSDRecord.php?LMID="

DEFAULTS = { 
  :bin_width => 5, 
  :bin_unit => :ppm,
  :search_unit => :ppm,
}

class Sample
  attr_accessor :file
  attr_accessor :spectrum
  def initialize(file, merge_opts={})
    @file = file
    @spectrum = merge_ms1_spectra(file, DEFAULTS.merge(merge_opts))
  end

  # returns a single spectrum object
  def merge_ms1_spectra(file, opts)
    spectra = []
    warn "using number of peaks as proxy for ms level right now"
    MS::Mzml.foreach(file) do |spectrum|
      spectra << spectrum if spectrum.mzs.size > 1000  # <<<<<<------ kludge for ms_level == 1
    end
    spectra.each {|spectrum| spectrum.sort! }

    MS::Spectrum.merge(spectra, opts)
  end
end

ext = ".lipidID.tsv"

parser = Trollop::Parser.new do
  banner "usage: #{File.basename(__FILE__)} [OPTIONS] <lipidmaps>.tsv <file>.mzML ..."
  text "output: <file>#{ext} ..."
  text ""
  text "note that sometimes you get an error from R like this:"
  text "(`eval': voidEval failed: Packet[cmd=2130771970,len=<nil>, con='<nil>', status=error...)"
  text "just re-run it and it will work"
  text ""
  opt :bin_width, "width of the bins for merging", :default => DEFAULTS[:bin_width]
  opt :bin_unit, "units for binning (ppm or amu)", :default => DEFAULTS[:bin_unit].to_s
  opt :search_unit, "unit for searching nearest hit (ppm or amu)", :default => DEFAULTS[:search_unit].to_s
  opt :top_n_peaks, "the number of highest intensity peaks to query the DB with", :default => 1000
  opt :display_n, "the number of best hits to display", :default => 20
  opt :verbose, "talk about it"
end

opts = parser.parse(ARGV)
opts[:bin_unit] = opts[:bin_unit].to_sym
opts[:search_unit] = opts[:search_unit].to_sym

if ARGV.size < 2
  parser.educate
  exit
end

(lipidmaps, *files) = ARGV

$VERBOSE = opts[:verbose]

proton = MS::Lipid::Modification.new(:proton)
h2o_loss = MS::Lipid::Modification.new(:water, :loss => true)

lipids = MS::LipidMaps.parse_file(lipidmaps)

ions = lipids.map do |lipid| 
  [[proton], [proton, h2o_loss]].map do |mods|
    MS::Lipid::Ion.new(lipid, mods)
  end
end.flatten(1)


searcher = MS::Lipid::Search.new(ions, :ppm => (opts[:search_unit] == :ppm))

files.each do |file|
  base = file.chomp(File.extname(file))
  puts "processing file: #{file}" if $VERBOSE
  sample = Sample.new(file, opts)

  num_points = sample.spectrum.mzs.size
  puts "#{num_points} merged peaks in #{file}" if $VERBOSE

  highest_points = sample.spectrum.points.sort_by(&:last).reverse[0,opts[:top_n_peaks]].sort

  sample.spectrum = MS::Spectrum.from_points( highest_points )

  queries = sample.spectrum.mzs.each_with_index.map {|mz,index| MS::Lipid::Search::Query.new(mz, index) }
  hit_groups = searcher.search(queries, :return_order => :sorted)

  hit_info = [:qvalue, :pvalue, :observed_mz, :theoretical_mz, :delta, :ppm]
  second_hit_info = [:ppm]

  output = base + ext
  puts "writing to #{output}" if $VERBOSE
  File.open(output, 'w') do |out|
    out.puts (hit_info + %w(2nd_hit_ppm first_isobar_name num_isobars isobars)).join("\t")
    hit_groups[0,opts[:display_n]].each_with_index do |hit_group,i|
      ar = []
      tophit = hit_group.first
      ar.push *hit_info.map {|mthd| tophit.send(mthd) }
      ar.push *second_hit_info.map {|mthd| hit_group[1].send(mthd) }
      common_name = tophit.db_isobar_group.first.lipid.common_name
      common_name = tophit.db_isobar_group.first.lipid.systematic_name if common_name == "-"
      ar.push common_name
      ar.push tophit.db_isobar_group.size
      ions = tophit.db_isobar_group.map do |ion|
        [ion.lipid.lm_id, ion.modifications.map do |mod| 
          (mod.gain? ? '+' : '-') + "(#{mod.charged_formula})"
        end.join
        ].join(":")
      end.join(' ')
      ar.push ions
      out.puts ar.join("\t")
    end
  end
end
