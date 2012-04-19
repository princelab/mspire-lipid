require 'spec_helper'

require 'mspire/sdf'

describe Mspire::SDF do

  before(:each) do
    @sdf_string = %Q{
benzene
ACD/Labs0812062058

6  6  0  0  0  0  0  0  0  0  1 V2000
1.9050   -0.7932    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
1.9050   -2.1232    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
0.7531   -0.1282    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
0.7531   -2.7882    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
-0.3987   -0.7932    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
-0.3987   -2.1232    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
2  1  1  0  0  0  0
3  1  2  0  0  0  0
4  2  2  0  0  0  0
5  3  1  0  0  0  0
6  4  1  0  0  0  0
6  5  2  0  0  0  0
M  END
$$$$
}
    @sdf_string.strip!
  end

  it 'reads SDF files' do
    sdf = Mspire::SDF.new(@sdf_string)
    sdf.header.should == ['benzene', 'ACD/Labs0812062058', '']
    sdf.atoms.size.should == 6
    sdf.atoms.first.should be_an(Mspire::SDF::Atom)
    sdf.bonds.size.should == 6
    sdf.bonds.first.should be_an(Mspire::SDF::Bond)
    p sdf.atoms
    p sdf.bonds
  end

end
