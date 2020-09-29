use Bio::TreeIO;
use Bio::Tree::TreeFunctionsI;
use FileUtil;

### usage looks like this:
### find readData/ -name \*subs.best.dnd | perl treeDists.pl - > suspect.trees
### output is input to extractValidSubSeqs.pl

my $list = shift;
my $lfh = FileUtil::openFileHandle($list);
while (<$lfh>) {
  chomp;
  my $file = $_;
  $n++;
  
  print STDERR "$n\n" unless $n % 100;
  my $treeio = Bio::TreeIO->new(-format => 'newick', -file => $file);
  while ( $tree = $treeio->next_tree ) {
    my @taxa = $tree->get_leaf_nodes;
    for (my $j = 1; $j < @taxa; $j++) {
      my $b = $taxa[$j]->id;
      push @dists,$tree->distance(-nodes => [$taxa[0], $taxa[$j]]);
    }
  }
  my $n = @dists;
  for (my $i = 0; $i < $n; $i++) { push @dists,0;}
  my $l = join ",",@dists;
  print $n,"\n",$l,"\n";
  
}
$lfh->close();
