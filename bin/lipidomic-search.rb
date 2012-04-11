#!/usr/bin/env ruby

require 'trollop'
require 'ms/mzml'
require 'ms/lipid/search'
require 'ms/lipid/ion'
require 'ms/lipid/search/query'
require 'ms/lipid_maps'
require 'ms/error_rate/qvalue'

# for html output: (just make the id clickable)
LIPIDMAPS_SEARCH = "http://www.lipidmaps.org/data/LMSDRecord.php?LMID="

DECOY_MODULATOR = 0.8319

DEFAULTS = { 
  :bin_width => 5, 
  :bin_unit => :ppm,
  :search_unit => :ppm,
}

def LipidPoint < Array
  attr_accessor :sample
end

class Sample
  attr_accessor :file
  attr_accessor :spectrum
  def initialize(file, merge_opts={})
    @file = file
    @spectrum = merge_ms1_spectra(file, DEFAULTS.merge(merge_opts))
  end

  # returns a single spectrum object
  def self.merge_ms1_spectra(files, opts)
    files.map do |file|
      MS::Mzml.foreach(file).select {|spec| spec.ms_level == 1 }.map(&:sort!)
    end
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
  text ""
  text "modifications (at least 1 charged mod is required):"
  opt :lithium, "search for lithium adducts"
  opt :ammonium, "search for ammonium adducts"
  opt :proton_gain, "search for proton gain"
  opt :proton_loss, "search for proton loss"
  opt :water_loss, "*all* mods are also considered with water loss"
  opt :decoy, "search with an equal number of decoy modifications"
  opt :verbose, "talk about it"
end

opts = parser.parse(ARGV)
opts[:bin_unit] = opts[:bin_unit].to_sym
opts[:search_unit] = opts[:search_unit].to_sym

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

(lipidmaps, *files) = ARGV

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


searcher = MS::Lipid::Search.new(ions, :ppm => (opts[:search_unit] == :ppm))

if opts[:decoy]
  # assumes a mod group that is either the mod or a mod and water loss
  decoy_ions = ions.map do |ion|
    # modify the first mod and leave the second untouched (if any)
    mod_group = ion.modifications
    fake_mod = mod_group.first.dup
    fake_mod.massdiff *= DECOY_MODULATOR
    fake_mod.formula = "FAKE#{mod_group.first.formula}(#{fake_mod.massdiff})"
    fake_mod.name = "fake_#{mod_group.first.name}".to_sym
    new_mod_group = [fake_mod, *mod_group[1..-1]]
    MS::Lipid::Ion.new(ion.lipid, new_mod_group)
  end
  decoy_searcher = MS::Lipid::Search.new(decoy_ions, :ppm => (opts[:search_unit] == :ppm)) 
end

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
  if opts[:decoy]
    decoy_hit_groups = decoy_searcher.search(queries, :return_order => :sorted)
    hit_group_qvalue_pairs = MS::ErrorRate::Qvalue.target_decoy_qvalues(hit_groups, decoy_hit_groups, :monotonic => true, &:pvalue)
    hit_group_qvalue_pairs.each do |hit_group, qval|
      hit_group.first.decoy_qvalue = qval
    end
  end

  # all info is relative to the hit_group
  info = {
    decoy_qvalue: :decoy_qvalue.to_proc,
    qvalue:  :qvalue.to_proc,
    pvalue:  :pvalue.to_proc,
    observed_mz:  :observed_mz.to_proc,
    theoretical_mz:  :theoretical_mz.to_proc,
    delta:  :delta.to_proc,
    ppm:  :ppm.to_proc,
    hit2_ppm: proc {|hg| hg[1].ppm },
    first_isobar_name: proc {|hg| (lipid=hg.first.db_isobar_group.first.lipid).common_name || lipid.systematic_name },
    num_isobars: proc {|hg| hg.first.db_isobar_group.size },
    ions: proc {|hg|
      hg.first.db_isobar_group.map do |ion|
        [ion.lipid.lm_id, ion.modifications.map do |mod| 
          (mod.gain? ? '+' : '-') + "(#{mod.charged_formula})"
        end.join
        ].join(":")
      end.join(' ')
    }
  }

  output = base + ext
  puts "writing to #{output}" if $VERBOSE
  File.open(output, 'w') do |out|
    out.puts info.keys.join("\t")
    hit_groups[0,opts[:display_n]].each do |hit_group|
      out.puts info.values.map {|prc| prc.call(hit_group) }.join("\t")
    end
  end

  if opts[:decoy]
    decoy_output = base + '.decoy' + ext
    File.open(decoy_output, 'w') do |dout|
      decoy_info = info.dup
      [:qvalue, :decoy_qvalue].each {|key| decoy_info.delete(key) }
      dout.puts decoy_info.keys.join("\t")
      decoy_hit_groups[0,opts[:display_n]].each do |hit_group|
        dout.puts decoy_info.values.map {|prc| prc.call(hit_group) }.join("\t")
      end
    end
  end
end
