$key = $ARGV[1];
my %keys;
open(FILE, $key);
while($buff = <FILE>){
	chomp $buff;
	if(substr($buff,0,1) eq ">"){
		$header = $buff;
		$header =~ m/>(.+)/g;
		$header = $1;
		
		$seq =<FILE>;
		chomp $seq;
		
		$keys{$header} = $seq;
		$total++;

	}
}

close(FILE);


$res = $ARGV[0];
open(FILE, $res);
while($buff = <FILE>){
       chomp $buff;
       if(substr($buff,0,1) eq "@"){
                $header = $buff;
                $header =~ m/@(.+)/g;
                $header = $1;

                $seq =<FILE>;
                chomp $seq;
	
		$compseq = $keys{$header};
		
		if(length($seq) == length($compseq)){
			$match++;
		}		
	
		if($seq eq $compseq){
			$acc++;
		}else{
			#print $seq . "\n";
			#print $compseq . "\n";
		}
		
		$ans++;
        }
}

close(FILE);

print STDERR "TOTAL QUERY: " . $total . "\n";
print STDERR "TOTAL SOLVED: " . $ans . "\n";
print STDERR "MATCH LENGTH: " . $match . "\n";
print STDERR "CORRECT SEQUENCE: " . $acc . "\n";

