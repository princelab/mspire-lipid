require 'spec_helper'

require 'mspire/lipid_maps'
require 'mspire/molecular_formula'
HAVE_RUBABEL = 
begin
  require 'rubabel' ; true
rescue ; false
end

describe Mspire::LipidMaps do
  describe 'parsing programmatically downloaded files' do

    before do
      @tfile = TESTFILES + '/lipidmaps_programmatic_short.tsv'
    end

    it 'parses lipid maps files' do
      lipids = Mspire::LipidMaps.parse_file(@tfile)
      lipids.size.should == 30  # one is rejected for no formula
      ll = lipids.last
      ll.sub_class.should == 'Isoflavonoids [PK1205]'
      ll.lm_id.should == "LMPK12050388"
      ll.formula.should be_a(Mspire::MolecularFormula)
      # ensures a high res mass by default
      ll.mass.should be_within(0.000000001).of(314.07903816634)
    end
  end

  describe 'parsing lipidmaps downloaded files' do

    before do
      @tfile = TESTFILES + '/lipidmaps_download.tsv'
      @tfile_sd = TESTFILES + '/lipidmaps_sd_download.tsv'
    end

    it 'parses lipid maps files' do

      [@tfile, @tfile_sd].each do |file|

        lipids = Mspire::LipidMaps.parse_file(file)
        lipids.size.should == 10
        ll = lipids.last
        ll.sub_class.should == '[pretend subclass]'
        ll.lm_id.should == 'LMFA00000014'
        ll.mass.should == 347.282416
        ll.formula.should be_a(Mspire::MolecularFormula)
        ll.formula[:c].should == 22  # <- frozen
        if file =~ /_sd_/
          ll.structure.include?('25 24  0  0  0  0  0  0  0  0999 V2000').should be_true
          if HAVE_RUBABEL
            lipids = Mspire::LipidMaps.parse_file(file, :rubabel_molecules => true)
            lipids.last.structure.should be_a(Rubabel::Molecule)
          end
        else
          ll.structure.should be_nil
        end

      end
    end
  end

end

