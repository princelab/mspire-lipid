# encoding: utf-8

require 'rubygems'
require 'rake'

require 'jeweler'
Jeweler::Tasks.new do |gem|
  # gem is a Gem::Specification... see http://docs.rubygems.org/read/chapter/20 for more options
  gem.name = "mspire-lipidomics"
  gem.homepage = "http://github.com/princelab/mspire-lipidomics"
  gem.license = "MIT"
  gem.summary = %Q{mass spectrometry based lipidomics - especially shotgun lipidomics}
  gem.description = %Q{does lipidomics}
  gem.email = "jtprince@gmail.com"
  gem.authors = ["John T. Prince"]
  gem.add_dependency "mspire", "~> 0.8.4"
  gem.add_development_dependency "rubabel", ">= 0.1.6"
  gem.add_development_dependency "rspec", "~> 2.3.0"
  gem.add_development_dependency "jeweler", "~> 1.6.4"
  gem.add_development_dependency "rcov", ">= 0"
end
Jeweler::RubygemsDotOrgTasks.new

require 'rspec/core'
require 'rspec/core/rake_task'
RSpec::Core::RakeTask.new(:spec) do |spec|
  spec.pattern = FileList['spec/**/*_spec.rb']
end

RSpec::Core::RakeTask.new(:rcov) do |spec|
  spec.pattern = 'spec/**/*_spec.rb'
  spec.rcov = true
end

task :default => :spec

require 'rdoc/task'
Rake::RDocTask.new do |rdoc|
  version = File.exist?('VERSION') ? File.read('VERSION') : ""

  rdoc.rdoc_dir = 'rdoc'
  rdoc.title = "ms-lipidomics #{version}"
  rdoc.rdoc_files.include('README*')
  rdoc.rdoc_files.include('lib/**/*.rb')
end
