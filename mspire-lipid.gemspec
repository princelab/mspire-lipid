# coding: utf-8
lib = File.expand_path('../lib', __FILE__)
$LOAD_PATH.unshift(lib) unless $LOAD_PATH.include?(lib)
require 'mspire/lipid/version'

Gem::Specification.new do |spec|
  spec.name          = "mspire-lipid"
  spec.version       = Mspire::Lipid::VERSION
  spec.authors       = ["John T. Prince"]
  spec.email         = ["jtprince@gmail.com"]
  spec.summary       = %q{mass spectrometry based lipidomics - especially shotgun lipidomics}
  spec.description   = %q{mass spectrometry based lipidomics - especially shotgun lipidomics.}
  spec.homepage      = ""
  spec.license       = "MIT"

  spec.files         = `git ls-files -z`.split("\x0")
  spec.executables   = spec.files.grep(%r{^bin/}) { |f| File.basename(f) }
  spec.test_files    = spec.files.grep(%r{^(test|spec|features)/})
  spec.require_paths = ["lib"]

  [
    ["mspire", "~> 0.10.8.0"],
    ["rubabel", ">= 0.1.6"],
  ].each do |args|
    spec.add_dependency(*args)
  end

  [
    ["bundler", "~> 1.6.2"],
    ["rake"],
    ["rspec", "~> 2.14.1"], 
    ["rdoc", "~> 4.1.1"], 
    ["simplecov", "~> 0.8.2"],
  ].each do |args|
    spec.add_development_dependency(*args)
  end
end
