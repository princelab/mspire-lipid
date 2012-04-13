require 'spec_helper'

require 'mspire/lipid'

describe Mspire::Lipid do

  before do
    @data = ['LMFA00000007', 'n-decanohydroxamic acid', 'N-hydroxydecanamide', 'C10H21NO2', 187.16, 'Fatty Acyls [FA]', 'Other Fatty Acyls [FA00]']
  end

  it 'can be initialized with an array' do
    lipid = Mspire::Lipid.new(*@data)
    lipid.mass.should == @data[4]
    lipid.sub_class.should == nil
  end
end
