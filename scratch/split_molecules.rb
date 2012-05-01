#!/usr/bin/env ruby

require 'yaml'
require 'openbabel'
include OpenBabel

if ARGV.size == 0
  puts "smiles and smart"
  exit
end

(smiles_string, smart_search) = ARGV

smi_to_svg = OBConversion.new
smi_to_svg.set_in_and_out_formats("smi", "svg")

obmol = OBMol.new

smi_to_svg.read_string(obmol, smiles_string)
obmol.add_hydrogens

p obmol.get_formula
svg_string = smi_to_svg.write_string(obmol).rstrip
File.write("before.svg", svg_string)

bond_triplets = []
(1..obmol.num_bonds).each do |i|
  bond = obmol.get_bond(i)
  if bond
    bond_triplets << [bond.get_begin_atom, bond.get_end_atom, bond]
  end
end

(1..obmol.num_atoms).each do |i|
  atom = obmol.get_atom(i)
  if atom.matches_smarts("C(=O)CC")
    bond_triplet = bond_triplets.find {|triplet| triplet.include?(atom) }
    obmol.delete_bond(bond_triplet.last)
    atom.begin_nbr_atom 
    atom.next_nbr_atom
    break
  end
end

molecules = obmol.separate
molecules.each_with_index do |mol,i|
  puts "mol #{i}"
  p mol.get_formula
  p mol.add_hydrogens
  p mol.get_formula
  svg_string = smi_to_svg.write_string(mol).rstrip
  File.write("after_#{i}.svg", svg_string)
end



#puts (obmol.methods - Object.new.methods).map(&:to_s).to_yaml
#p obmol.exact_mass




=begin
$obConversion->ReadFile($obMol, $filename);

for (1..$obMol->NumAtoms()) {
    $atom = $obMol->GetAtom($_);
    # look to see if this atom is a thiophene sulfur atom
    if ($atom->MatchesSMARTS("[#16D2]([#6D3H1])[#6D3H1]")) {
	$sulfurIdx = $atom->GetIdx();
    # see if this atom is one of the carbon atoms bonded to a thiophene sulfur
    } elsif ($atom->MatchesSMARTS("[#6D3H1]([#16D2][#6D3H1])[#6]") ) {
	if ($c2Idx == 0) { $c2Idx = $atom->GetIdx(); }
	else {$c5Idx = $atom->GetIdx(); }
    }
}

# Get the actual atom objects -- indexing will change as atoms are added and deleted!
$sulfurAtom = $obMol->GetAtom($sulfurIdx);
$c2Atom = $obMol->GetAtom($c2Idx);
$c5Atom = $obMol->GetAtom($c5Idx);

$obMol->DeleteAtom($sulfurAtom);

$obMol->DeleteHydrogens($c2Atom);
$obMol->DeleteHydrogens($c5Atom);

$c2Atom->SetAtomicNum(1);
$c5Atom->SetAtomicNum(1);

$obConversion->WriteFile($obMol, "$filename.mol");

=end
