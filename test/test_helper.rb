require 'stringio'
require 'test/unit'

ROOT_PATH=File.dirname(__FILE__)
$: << File.expand_path(File.join(ROOT_PATH, "../lib/"))
$: << File.expand_path(File.join(ROOT_PATH, "../lib/NetAnalyzer/"))

require File.join(ROOT_PATH, '../lib/NetAnalyzer/network')