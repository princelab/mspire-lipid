require 'spec_helper'

require 'ms/lipidomics'

describe "simple DB search" do

  before do
    @tfile = TESTFILES + '/lipidmaps_short.tsv'
  end
  it 'parses lipid maps files' do
    LipidMaps
  end


  it 'creates a prob distribution' do
    MS::Lipid::Search.simple_search( MS::LipidMaps.parse_file(@tfile) )
  end
end
