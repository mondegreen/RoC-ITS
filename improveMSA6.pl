use FileUtil;

my $msa = shift;
my $word = shift;
my $minCov = shift;
$word = 12 unless $word;
my $gggs = "n" x 22;

my $maxParallel = 40;

my $mfh = FileUtil::openFileHandle("perl reformatFasta.pl $msa 0 99999999|");
while (<$mfh>) {
  if (/^>(\S+)/) {
    my $id = $1;
    $_ = <$mfh>;
    chomp; 
    $seq{$id} = $gggs.uc($_).$gggs;
    my @seq = split //,$seq{$id};
    $len = @seq;
    push @seqs1,\@seq;
    push @ids,$id;
#    print $seq{$id},"\n";
  }
}
$mfh->close();

my @good;
my $cnt = @ids;
for (my $i = 0; $i < $len; $i++) {
  my $bases = 0;
  for (my $j = 0; $j < $cnt; $j++) {
    $dist{${$seqs1[$j]}[$i]}++;
    if (${$seqs1[$j]}[$i] eq "-") {
    } else {
      $bases++;
    }
  }
  if ($bases < $minCov) {
    for (my $j = 0; $j < $cnt; $j++) {
      ${$seqs1[$j]}[$i] = "";
    }
  }
}
for (my $j = 0; $j < $cnt; $j++) {
  my $tmp = join "",@{$seqs1[$j]};
  my @seq = split //,$tmp;
  $len = @seq;
  push @seqs,\@seq;
}

for (my $i = 0; $i < $len; $i++) {
  my %dist;
  for (my $j = 0; $j < $cnt; $j++) {
    $dist{${$seqs[$j]}[$i]}++;
  }
  if ($dist{"-"}/$cnt >= 0.1) {
    push @good,"-";
  } else {
    my $best;
    foreach my $j (sort {$dist{$b}<=>$dist{$a}; } keys %dist) {
      $best = $j;
      last;
    }
    if ($dist{$best}/$cnt < 0.9) {
      push @good,"-";
    } else {
      push @good,"G";
    }
  }
}

my $g = join "",@good;
#print $g,"\n";
if ($g =~ /^(-+)/) {
  my $w = "p" x length($1);
  $g =~ s/^(-+)/$w/;
}
my @gg = split //,$g;
my @regionBegins;
my @regionEnds;
my $start = 0;
my $endFound = 0;
while (!$endFound) {
  my $i = $start;
  for ($i = $start; $i < @gg; $i++) {
    if ($gg[$i] eq "-") {
      $begin = $i - $word;
      my $tot = 0;
      my $pos = $i;
      while ($tot < $word) {
        if ($gg[$pos] eq "G") {
          $tot++;
        } else {
          $tot = 0;
        }
        $pos++;
      }
      $end = $pos;
      push @regionBegins,$begin;
      push @regionEnds,$end;
      $start = $end - $word;
      last;
    }
  }
  if ($i < @gg) {
  } else {
    $endFound = 1;
  }
}

for (my $ii = 0; $ii < @regionBegins; $ii++) {
  my $fho = new FileHandle(">temp.$$.$ii.$regionBegins[$ii].$regionEnds[$ii].fa");
  for (my $i = 0; $i < $cnt; $i++) {
    my @set;
    for (my $j = $regionBegins[$ii]; $j < $regionEnds[$ii]; $j++) {
      my $chr = ${$seqs[$i]}[$j];
      if ($j < ($regionEnds[$ii]-$word)) {
        ${$seqs[$i]}[$j] = "";
        push @set,$chr;
      } else {
        ### subsequences that have gaps in the right anchor region are by definition rare but need to
        ### be replaced by N's. This is important as later non-dashes in the right anchor are removed
        ### to when the new msa is inserted.
        if ($chr eq "-") {
          push @set,"n";
        } else {
          push @set,$chr;
        }
      }
    }
    my $sub = join "",@set;
    $sub =~ s/-+//g;
    print $fho ">$ids[$i] $regionBegins[$ii] $regionEnds[$ii]\n";
    print $fho $sub,"\n";
  }
  $fho->close();
  ### requires probcons installed and in path
  my $probCmd = "probcons -c 0 temp.$$.$ii.$regionBegins[$ii].$regionEnds[$ii].fa |perl reformatFasta.pl - 0 9999999 > temp.$$.$ii.$regionBegins[$ii].out &";
  push @cmds,$probCmd;
  print STDERR $probCmd,"\n";
}
my $ccnt = 0;
my $fho = new FileHandle(">temp.$$.cmd");
foreach my $i (@cmds) {
  print $fho $i,"\n";
  $ccnt++;
  if ($cnt == $maxParallel) {
    print $fho "wait\n";
    $ccnt = 0;
  }
}
print $fho "wait\n";
$fho->close();

my $cmd = "tcsh -x temp.$$.cmd > /dev/null";
print STDERR $cmd,"\n";
my $run = `$cmd`;

for (my $ii = 0; $ii < @regionBegins; $ii++) {
  my @result;
  my $rfh = FileUtil::openFileHandle("temp.$$.$ii.$regionBegins[$ii].out");
  while (<$rfh>) {
    next if /^>/;
    chomp;
    push @result,$_;
  }
  $rfh->close();
  for (my $i = 0; $i < $cnt; $i++) {
    for (my $j = $regionBegins[$ii]; $j < $regionEnds[$ii]; $j++) {
      my @res = split //,$result[$i];
      my $e = @res;
      my $removed = 0;
      for (my $k = $e - 1; $k >= 0; $k--) {
        if ($res[$k] =~ /\w/o) {
          $res[$k] = "";
          $removed++;
        }
        last if $removed == $word;
      }
      my $res = join "",@res;
      ${$seqs[$i]}[$j] = lc($res);
      last;
    }
  }
}

for (my $i = 0; $i < $cnt; $i++) {
  my $l = join "",@{$seqs[$i]};
  $l =~ s/n/-/g;
  print ">$ids[$i]\n$l\n";
}
$cmd = "rm temp.$$*.fa";
system($cmd);
$cmd = "rm temp.$$*.out temp.$$.cmd";
system($cmd);
