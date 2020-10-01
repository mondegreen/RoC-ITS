use Bio::TreeIO;
use Bio::Tree::TreeFunctionsI;
use FileUtil;
use BioStats;

my $list = shift;
my $cutoff = shift;
$cutoff = 0.05 unless defined($cutoff);

my $lfh = FileUtil::openFileHandle($list);
while (<$lfh>) {
  chomp;
  my $file = $_;
  $n++;
  
  my @dists;
  my $cnt;
  print STDERR "$n\n" unless $n % 100;
  my $treeio = Bio::TreeIO->new(-format => 'newick', -file => $file);
  while ( $tree = $treeio->next_tree ) {
    my @taxa = $tree->get_leaf_nodes;
    $cnt = @taxa;
    for (my $j = 1; $j < @taxa; $j++) {
      my $b = $taxa[$j]->id;
      push @dists,$tree->distance(-nodes => [$taxa[0], $taxa[$j]]);
    }
    for (my $j = 2; $j < @taxa; $j++) {
      my $b = $taxa[$j]->id;
      push @dists,$tree->distance(-nodes => [$taxa[1], $taxa[$j]]);
    }
  }
  my $mean = BioStats::getMean(\@dists);
  my $std = BioStats::getMeanStandardDeviation(\@dists,$mean);
  my $x = @taxa;
  my $l = join " ",$file,$mean,$std,$cnt,"xx",@dists;
  print $l,"\n" if $std > $cutoff;
}
$lfh->close();
