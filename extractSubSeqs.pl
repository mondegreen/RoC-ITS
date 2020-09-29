use FileUtil;

### 

my $f = shift;        ### hmmer output file
my $fasta = shift;    ### nanopore fasta file
my $outDir = shift;   ### directory where the subseqs will be put

my $fafh = FileUtil::openFileHandle($fasta);
while (<$fafh>) {
  if (/^>(\S+)/) {
    $id = $1;
    $_ = <$fafh>;
    chomp;
    $seq{$id} = $_;
  }
}
$fafh->close();

my $fh = FileUtil::openFileHandle($f);
while (<$fh>) {
  if (/^>>\s*((\w\w)\S+)/) {
    my $begin = 0;
    ($id,$index) = ($1,$2);
    my $cmd = "mkdir -p $outDir/$index/";
    system($cmd);
    my $fho = new FileHandle(">$outDir/$index/$id.subs.fa");
    while (<$fh>) {
      last unless /\S+/;
      my @d = split /\s+/;
      next unless $d[2] eq "!";
      my $end = $d[11];
      ### this is bad - hard coded length of hmm used to identify the splint region
      my $hmmEnd = 386-$d[8];
      if ($d[8] < 300) {
        $hmmEnd = 0;
      }
      my $delta = $end-$begin+$hmmEnd;
      if ($delta > 1500 && $delta < 3500) {
        my $sub = substr($seq{$id},$begin,$end-$begin+$hmmEnd);
        my $bb = sprintf "%9.9d",$begin;
        my $ee = sprintf "%9.9d",$end;
        print $fho ">$id:$bb-$ee /begin=$begin /end=$end /length=$delta\n$sub\n";
      }

      $begin = $end;
    }
    $fho->close();
  }
}
$fh->close();
