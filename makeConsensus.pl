use FileUtil;

my $fasta = shift;
my $name = $fasta;
$name =~ s/.*\///;
$name =~ s/.fasta_aln$//;
my $ffh = FileUtil::openFileHandle("perl reformatFasta.pl $fasta 0 999999999 1|");
while (<$ffh>) {
	if (/^>(\S+)/) {
		$id = $1;
		$cnt++;
	} else {
		chomp;
		my $ll = length($_);
		my ($dash) = $_ =~ tr/-/-/;
		$totalChar += $ll;
		$totalDash += $dash;
		my @s = split //;
		my $l = @s;
		$len = $l if $l > $len;
		push @seq,\@s;
		push @ids,$id;
	}
}
$ffh->close();

my @con;
for (my $j = 0; $j < $len; $j++) {
	my %b;
	for (my $i = 0; $i < @ids; $i++) {
		$b{${$seq[$i]}[$j]}++;
	}
	foreach my $i (sort { $b{$b} <=> $b{$a};} keys %b) {
		push @con,$i;
		last;
	}
}
my $con = join "",@con;
$con =~ s/-//g;
$con =~ s/ //g;
my $perDash = int($totalDash/$totalChar*1000)/10;
print ">$name /count=$cnt /dash=$perDash\n";
print $con,"\n";
