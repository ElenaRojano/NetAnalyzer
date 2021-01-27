# coding: utf-8
lib = File.expand_path('../lib', __FILE__)
$LOAD_PATH.unshift(lib) unless $LOAD_PATH.include?(lib)
require 'NetAnalyzer/version'

Gem::Specification.new do |spec|
  spec.name          = "NetAnalyzer"
  spec.version       = NetAnalyzer::VERSION
  spec.authors       = ["Elena Rojano, Pedro Seoane"]
  spec.email         = ["elenarojano@uma.es, seoanezonjic@hotmail.com"]

  spec.summary       = %q{Network analysis tool that calculate and validate different association indices.}
  spec.description   = %q{NetAnalyzer is a useful network analysis tool developed in Ruby that can 1) analyse any type of unweighted network, regardless of the number of layers, 2) calculate the relationship between different layers, using various association indices (Jaccard, Simpson, PCC, geometric, cosine and hypergeometric) and 3) validate the results}
  spec.homepage      = "https://github.com/ElenaRojano/NetAnalyzer"
  spec.license       = "MIT"

  spec.files         = `git ls-files -z`.split("\x0").reject { |f| f.match(%r{^(test|spec|features)/}) }
  spec.bindir        = "bin"
  spec.executables   = spec.files.grep(%r{^bin/}) { |f| File.basename(f) }
  spec.require_paths = ["lib"]

  spec.add_development_dependency "rake"
  spec.add_development_dependency "rspec"
  spec.add_dependency "numo-linalg"
  spec.add_dependency "numo-narray"
  spec.add_dependency "pp"
  spec.add_dependency "pycall"
  spec.add_dependency "npy"
  spec.add_dependency "bigdecimal"
  spec.add_dependency "gv"
  spec.add_dependency "semtools"
end
