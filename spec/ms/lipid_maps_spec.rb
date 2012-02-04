require 'spec_helper'

require 'ms/lipid_maps'

describe MS::LipidMaps do

  before do
    @tfile = TESTFILES + '/lipidmaps_short.tsv'
  end

  it 'parses lipid maps files' do
    lipids = MS::LipidMaps.parse_file(@tfile)
    lipids.size.should == 30  # one is rejected for no formula
    ll = lipids.last
    ll.sub_class.should == 'Isoflavonoids [PK1205]'
    ll.lm_id.should == "LMPK12050388"
  end
end

