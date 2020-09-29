use FileUtil;

### Takes clustering/tree data from existing sets of sequences and identifies outliers that suggest that there are multiple inserts in a RoC-ITS read
### If the data is consistent with multiple inserts, the data is broken up into independent read sets.

my $treeData = shift;
my $treeDir = shift;
my $suspectTrees = shift;
my $validSubTrees = shift;
my $outDir = shift;

my $vfh = FileUtil::openFileHandle($validSubTrees);
while (<$vfh>) {
  chomp;
  my @d = split /\t/;
  $valid{$d[0]} = $d[1];
}
$vfh->close();

my $sfh = FileUtil::openFileHandle($suspectTrees);
while (<$sfh>) {
  chomp;
  my @d = split /\s+/;
  $d[0] =~ s/.*\///;
  $d[0] =~ s/\.subs.*//;
  if (exists($valid{$d[0]})) {
  } else {
    $bad{$d[0]} = 1;
  }
}
$sfh->close();

my $tfh = FileUtil::openFileHandle($treeData);
while (<$tfh>) {
  chomp;
  s/^.*?\"//;
  s/\"//;
  s/_/:/;
  s/_/ /;
  my @d = split /\s+/;
  $d[1] =~ /(.*):/;
  my $id = $1;
  next unless exists($valid{$id}) && $valid{$id} == $d[0];
  if (@l) {
    if ($id eq $lid && $d[3] == $l[3]) {
      push @set,\@d;
    } else {
      binSet(\@set);
      undef @set;
      push @set,\@d;
    }
  } else {
    push @set,\@d;
  }
  $lid = $id;
  @l = @d;
}
$tfh->close();
binSet(\@set);

### by targeting best.fas files we can look at those datasets with enough data to make trees
my $fafh = FileUtil::openFileHandle("find $treeDir -name \*.subs.best.fas|");
while (<$fafh>) {
  chomp;
  my $basename = $_;
  $basename =~ s/subs.*/subs.fa/;
  $basename =~ /.*\/(.*)\.subs.fa$/;
  my $id = $1;
  $id =~ /^(..)/;
  my $code = $1;
  my $cmd = "mkdir -p $outDir/$code/";
  system($cmd);
  next if exists($bad{$id});
  if (exists($valid{$id})) {
    my %map;
    my $fhin = FileUtil::openFileHandle($basename);
    while (<$fhin>) {
      if (/^>(\S+)/) {
        my $read = $1;
        my $def = $_;
        my $seq = <$fhin>;
        $map{$read} = $def.$seq;
      }
    }
    $fhin->close();
    for (my $i = 1; $i <= $valid{$id}; $i++) {
      my $fhout = new FileHandle(">$outDir/$code/$id.$i.subs.fa");
      foreach my $j (keys %map) {
        if ($clusters{$id}{$j}{$i} > 5) {
          print $fhout $map{$j};
        }
      }
      $fhout->close();
    }
  } else {
    $cmd = "cp $basename $outDir/$code/";
    system($cmd);
  }
}
$fafh->close();

sub binSet {
  my ($s) = @_;

  my $cnt;
  my $name;
  my @set = @{$s};
  foreach my $d (@set) {
    my @d = @{$d};
    if ($d[4] == 0) {
      $cnt = 1;
      $d[1] =~ /(.*):/;
      $name = $1;
    }
    if (exists($assign{$d[5]})) {
    } else {
      $assign{$d[5]} = $cnt;
      $cnt++;
    }
  }
  foreach my $d (@set) {
    my @d = @{$d};
    $clusters{$name}{$d[1]}{$assign{$d[5]}}++;
  }
}

