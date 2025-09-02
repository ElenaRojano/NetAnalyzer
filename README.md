# NetAnalyzer

DEPRECATED PROJECT. MIGRATED TO [python semtools](https://github.com/seoanezonjic/NetAnalyzer)
 
NetAnalyzer is a network analysis tool that can be used to calculate the associations between nodes in unweighted n-partite networks [1]. The calculation of the association between nodes is based on similarity indices (Jaccard, Simpson, geometric and cosine), statistic-based (Pearson correlation coefficient, CSI and  hypergeometric) [2] and a special metric designed only for tripartite networks (here called as 'transference' method [3]). The user can choose the association index method according to the network to analyse. The tool gives a table of results, with all the associations between nodes and the association value calculated.
 
If you use this tool, please cite us: [1] E. Rojano, P. Seoane, A. Bueno, J. R. Perkins & J. A. G. Ranea. Revealing the Relationship Between Human Genome Regions and Pathological Phenotypes Through Network Analysis. Lecture Notes in Computer Science, Vol 10208, 197-207 (2017).

[2] Fuxman-Bass et al. Using networks to measure similarity between genes: association index selection. Nature Methods, 10(12):1169-76. 2013.

[3] Alaimo et al. ncPred: ncRNA-Disease Association Prediction through Tripartite Network-Based Inference. Frontiers in Bioengineering and Biotechnology, 2:71, 2014.

## Installation

Linux & MacOS:

Please, check before your Ruby compiler (it has to be clang to install nmatrix) 

```
ruby -rrbconfig -e 'puts RbConfig::MAKEFILE_CONFIG["CC"]'
```

If not, install RVM (https://rvm.io/), and then: 

```
rvm reinstall 2.4.1 --with-gcc=clang --with-cxx=clang++
git clone https://github.com/SciRuby/nmatrix.git
cd nmatrix
gem install bundler
bundle install
bundle exec rake compile
bundle exec rake spec
```

Once nmatrix gem is installed:

````ruby
gem install 'NetAnalyzer'
```

## Usage

The program NetAnalyzer.rb can analyse an unweighted network to calculate the association index between different nodes. 

An example of use can be the following:

    $ NetAnalyzer.rb NetAnalyzer.rb -i network.txt -l 'hpo,HP:;patients,[0-9]' -m hypergeometric -u 'hpo;patients' -a 'associations_file.txt'

Where:

```
-i: Input file with the network to analyse. It must have two columns (separated by default by tabs) that represents the nodes that are related (NodeA\tNodeB). Please if you have doubts about the format, check the example providen.
-l: Layers construction. Please consider that, depending on the n-partite network you provide, NetAnalyzer will transform it into a bipartite one to perform the analysis (excepting if the association method used is 'transference'). The layers must contain a identifier of the node, and a character or pattern to identify. In this example, the bipartite network has HPO terms (with 'HP:' string in each of them) and patients that have these HPO terms (they are given as numerical patient IDs). Both layers must be separated by ';'.
-m: Association method. There are 8 different association methods to choose: 'jaccard', 'cosine', 'pcc', 'csi', 'hypergeometric', 'simpson', 'geometric' and 'transference'.
-u: Set which layer will be the one that establish connections between nodes in the other layer. In this case, we will get with patient is associated to other patient because the HPO they share.
-a: Associations output file name. Here you can find the associations between nodes in the network and the calculated association value, according to the chosen method.
``` 

Optional flags:

```
-s: Split character. Change if the layers of the network are not separated by tabs.
-o: Output file name.

```

## Development

After checking out the repo, run `bin/setup` to install dependencies. Then, run `rake spec` to run the tests. You can also run `bin/console` for an interactive prompt that will allow you to experiment.

To install this gem onto your local machine, run `bundle exec rake install`. To release a new version, update the version number in `version.rb`, and then run `bundle exec rake release`, which will create a git tag for the version, push git commits and tags, and push the `.gem` file to [rubygems.org](https://rubygems.org).

## Contributing

Bug reports and pull requests are welcome on GitHub at https://github.com/ElenaRojano/NetAnalyzer.


## License

The gem is available as open source under the terms of the [MIT License](http://opensource.org/licenses/MIT).

