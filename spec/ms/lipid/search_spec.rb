require 'spec_helper'

require 'ms/lipid_maps'

describe MS::Lipid::Search do
  before do
    @lipids = MS::LipidMaps.parse_file(TESTFILES + '/lipidmaps_short.tsv')
  end

  it 'creates a probability distribution'

end
