require "bundler/gem_tasks"
#require "rspec/core/rake_task"

#RSpec::Core::RakeTask.new(:spec)

#task :default => :spec

require 'rake/testtask'

Rake::TestTask.new do |t|
  t.libs << 'test'
  t.pattern = "test/*_test.rb"
end