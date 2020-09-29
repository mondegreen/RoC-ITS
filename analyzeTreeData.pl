use FileUtil;
use Math::CDF;

my $f = shift;
my $format = shift;

$format = 0 unless defined($format);
my %results;
my %successTotal;
my %totalCnts;
my %lenTotal;

my $fh = FileUtil::openFileHandle($f);
while (<$fh>) {
  my $l = $_;
  chomp;
  s/.*?\"//;
  s/\"//;
  s/_(\d+)-(\d+)_/ $1 $2 /;
  my @d = split /\s+/;
  if (@l) {
    if ($d[0] == $l[0] && $d[1] eq $l[1] && $d[5] == $l[5]) {
      push @set,\@d;
    } else {
      analyzeSet(\@set);
      undef @set;
      push @set,\@d;
    }
  } else {
    push @set,\@d;
  }
  @l = @d;
}
$fh->close();
analyzeSet(\@set);

if ($format) {
} else {
  foreach my $i (keys %results) {
    foreach my $j (sort keys %{$results{$i}}) {
      my $p = $results{$i}{$j}/$totalCnts{$i}{$j};
      my $s = $successTotal{$i}{$j}/$totalCnts{$i}{$j};
      my $t = $lenTotal{$i}{$j}/$totalCnts{$i}{$j};
      print "$i\t$j\t$t\t$s\t$p\n" if ($t+1) >= ($j*3) && $p < 0.001;
    }
  }
}

sub analyzeSet {
  my ($s) = @_;

  my $size;
  my $name;
  my @set = @{$s};
  my @pattern;
  foreach my $d (@set) {
    my @d = @{$d};
    my $l = join "\t",@d;
    if ($d[6] == 0) {
      push @pattern,0;
      $size = $d[0];
      $name = $d[1];
    } else {
      push @pattern,$d[7];
    }
  }
  my %setCnts;
  for (my $i = 1; $i < @pattern; $i++) {
    my %grp;
    my $max = @pattern;
    $max = ($i+$size) < $max ? ($i+$size):$max;
    for (my $j = $i; $j < $max; $j++) {
      $grp{$pattern[$j]}++;
      push @grp,$pattern[$j];
    }
    my $x = keys(%grp);
    if (keys(%grp) == $size) {
      my $sub = join "",sort(keys(%grp));
#      print "hhh $size $i $x $sub\n";
      $setCnts{$sub}++;
    }
  }
  my $best;
  foreach my $i (sort {$setCnts{$b}<=>$setCnts{$a}; } keys %setCnts) {
    $best = $i;
  }
  my $len = @set - 1;
  my $successes = $setCnts{$best} + $size - 1;
  if (defined($best)) {
    my $p = 1-Math::CDF::pbinom($successes,$len,1/$size);
    print "$name\t$size\t$len\t$best\t$successes\t$setCnts{$best}\t$p\n" if $format == 1;
    $results{$name}{$size}+=$p;
    $successTotal{$name}{$size}+=$successes;
  } else {
    print "$name\t$size\t$len\t-\t0\t0\t1\n" if $format == 1;
    $results{$name}{$size}+=1;
  }
  $lenTotal{$name}{$size}+=$len;
  $totalCnts{$name}{$size}++;
}
