$err = 2;
$length = 100;
$fragment = 180;
$dev = 0.2;
$count = 1000000;

$out = "reads_sim_acc_len" . $length . "bp_frag" . $fragment . "bp_err" . $err . "perc";
$r1 = $out . ".R1.fastq";
$r2 = $out . ".R2.fastq";
$key = $out . ".key.fasta";

open(OUT1, ">" . $r1);
open(OUT2, ">" . $r2);
open(KEY, ">" . $key);


my @bases;
$bases[0] = "A";
$bases[1] = "C";
$bases[2] = "G";
$bases[3] = "T";

my %revcomp;
$revcomp{"A"} = "T";
$revcomp{"C"} = "G";
$revcomp{"G"} = "C";
$revcomp{"T"} = "A";

$qual = "C" x $length;

for($i=0;$i<$count;$i++){

	$d = $dev * $fragment;
	$r = rand($d+1);
	$d = $r -$d;
	$act_frag = int($fragment + $d);

	#print $act_frag . "\t" . $d . "\n";

	$frag_seq = "";
	for($j=0;$j<$act_frag;$j++){
		$r = int(rand(4));
		$base = $bases[$r];
		$frag_seq .= $base;
	}

	#$seq1 = substr($frag_seq,0,$length);
	$seq1 = "";
	for($k=0;$k<$length;$k++){
		$base = substr($frag_seq,$k,1);

		$odd = int(rand(1000));
		if($odd <= (10*$err)){
			$seq1 .= $revcomp{$base};
		}else{
			$seq1 .= $base;
		}
	}

	$seq2 = substr($frag_seq,length($frag_seq) - $length, $length);
	$seq2_a = reverse($seq2);
	$rev_seq2 = "";
	for($k=0;$k<length($seq2);$k++){
		$base = $revcomp{substr($seq2_a,$k,1)};

                $odd = int(rand(1000));
                if($odd <= (10*$err)){
                        $rev_seq2 .= $revcomp{$base};
                }else{
			$rev_seq2 .= $base;
		}
	}

	$header = "reads_" .$i . "_" . $act_frag . "\n";
	print KEY ">" . $header;
	print KEY $frag_seq . "\n"; 
	
	print OUT1 "@" . $header;
	print OUT1  $seq1 . "\n";
	print OUT1 "+\n";
	print OUT1 $qual . "\n";

	print OUT2 "@" . $header;
	print OUT2  $rev_seq2 . "\n";
        print OUT2 "+\n";
        print OUT2 $qual . "\n";
	


}

close(OUT1);
close(OUT2);
close(KEY);
