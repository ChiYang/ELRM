#!/usr/bin/perl
#Script Apriori.pl
#Author: Chi Yang
#Date: 2014/05/25
#License: MIT licence
$minSup = $ARGV[0];
$DataFile = $ARGV[1];
$Output = $ARGV[2];
#Example: perl Apriori_frequently_cooccurred_single_base_features.pl 0.2 ../example/similarMotifs/PurR.fasta ../example/Frequent_SingleBase_Sets/PurR.apriori_output
if(@ARGV != 3){
	print STDERR "Three arguments have to be set.\n";
	print STDERR "Usage: perl Apriori.pl [minimum support] [input fasta file] [output file]\n";
	exit;
}elsif($minSup <= 0 or $minSup >= 1){
	print STDERR "Error: The range for the minimum support should be larger than 0 and less than 1\n";
	exit;
}elsif(!-e $DataFile){
	print STDERR "The input fasta file '$DataFile' does not exist.\n";
	exit;
}

my($startDate, $startTime) = &gettime();
print "StartTime: $startDate $startTime\n";
open(OUTPUT, ">$Output");
close(OUTPUT);
@T = (); #Each sequence is considered as an Transactions
$k = 1; #Frequent K-item set(s); start from 1
@C = ();
%C = ();
@F = ();
%F = ();

%sup = ();

$totalNumber = 0;
$totalScore = 0;


print "Reading the data...$DataFile\n";
open(FILE, $DataFile);
while(<FILE>){
	chomp;
	#@transaction = split(/\t/);
	
	if($_ =~ /^>(.*)$/){
		next;
	}
	@nucleotides = split('');
	@{$T[$i]} = ();
	for($j = 0; $j < @nucleotides; $j ++){
		$nuc = $nucleotides[$j];
		$pos = $nuc.(1+$j);
		#$pos = (1+$j).$nuc;
		push(@{$T[$i]}, $pos);
	}
	$avgScore = 1;	
	for($j = 0; $j < @nucleotides; $j ++){
		$nuc = $nucleotides[$j];
		$pos = $nuc.(1+$j);
		#$pos = (1+$j).$nuc;
		$C[$k-1]{$pos} += 1;
	}
	$totalScore += $avgScore;
	$i ++;
	$totalNumber++;
	print "Sequence number: $totalNumber\r";
}
close(FILE);
print "\nReading Finished\n";
print "Execute Apriori algorithm to search frequently co-occurred single-base features...\n";

foreach $pos (keys(%{$C[$k-1]})){
	$count = $C[$k-1]{$pos};
	if($count/$totalScore >= $minSup){
		$F[$k-1]{$pos} = $count;
	}
}
do{
	$k = $k + 1;
	print "Searching for frequent $k itemsets,\t";
	apriori_gen();
	print "Candidate itemsets: ".scalar(keys(%{$C[$k-1]}))."\n";
	if(scalar(keys(%{$C[$k-1]})) > 0){
		for($i = 0; $i < $totalNumber; $i ++){
			print "$i / ".$totalNumber."\r";
			MySubsets(\@{$T[$i]});
		}
		print "\n";
		$nextrun = 0;
		foreach $pos (keys(%{$C[$k-1]})){
			$count = $C[$k-1]{$pos};
			if($count/$totalScore >= $minSup){
				$F[$k-1]{$pos} = $count;
				$nextrun ++;
			}
		}
		print "Having ".$nextrun." frequent $k\-item itemsets.\n";
	}else{
		$nextrun = 0;
		print "Finished\n"
	}
}while($nextrun > 0);

print "Output to files: $Output\n";
open(OUTPUT, ">$Output");
for($i =0; $i < @F; $i++){
	foreach $p (keys(%{$F[$i]})){
		$support = $F[$i]{$p}/$totalScore;
		print OUTPUT "$p\t$F[$i]{$p}\t$support\n";
		$count ++;
	}
}
close(OUTPUT);
my($endDate, $endTime) = &gettime();

print "EndTime: $endDate $endTime\n";


#########

sub apriori_gen{
	#Genreate %{C[$k-1]} from %{$F[$k-2]}  #$k = 0 => 1 item itemset
	@temp = sort {$a <=> $b}(keys(%{$F[$k-2]}));
	#print join(",", @temp)."\n";
	for($i = 0; $i < @temp; $i++){
		$a = $temp[$i];
		@elements = split(',',$a);
		$lastA = pop(@elements);
		$testA = join(',', @elements);
		for($j = $i+1; $j < @temp; $j++){
			$b = $temp[$j];
			@elements2 = split(',',$b);
			$lastB = pop(@elements2);
			$canMerge = 1;
			$testB = join(',', @elements2);
			if($testA ne $testB){
				$canMerge = 0;
			}
			if($canMerge == 1){
				if($lastA =~ /(\d+)/){
					$posA = $1;
				}
				if($lastB =~ /(\d+)/){
					$posB = $1;
				}
				if($posA > $posB){
					push(@elements2, $lastB);
					push(@elements2, $lastA);
					
				}else{
					push(@elements2, $lastA);
					push(@elements2, $lastB);
				}
				$newCombination = join(",", sort{$a <=> $b}(@elements2));
				#Maybe the same position can be further pruned
				$C[$k-1]{$newCombination} = 0;
			}
		}
	}
}	

sub MySubsets{
	@pos = @{$_[0]};
	$transactionLength = scalar(@pos);
	@subsets = ();
	@{$subsets[0]{"All"}} = @pos;
	for($level = 1; $level < $k; $level ++){
		foreach $key (sort{$a <=> $b}(keys(%{$subsets[$level-1]}))) {
			$temp = \@{$subsets[$level-1]{$key}};
			while(scalar(@{$temp}) >= $k-$level+1){
				if($key eq "All"){
					#no prefix
					@prefix = ();
				}else{
					@prefix = split(',',$key);
				}
				$pre = shift(@{$temp});
				push(@prefix, $pre);
				$prefixLength = @prefix;
				if($prefixLength == $k - 1){
					$prefix = join(',', @prefix);
					foreach $t (@{$temp}){
						$out = $prefix.','.$t;
						if(defined($C[$k-1]{$out})){
							$C[$k-1]{$out} ++;
						}
					}
				}else{
					$prefix = join(",", @prefix);
					@{$subsets[$level]{$prefix}} = @{$temp};
				}
			}
		}
		#level = prefix number - 1
	}
	#return @output;
}
sub gettime{
	my ($sec, $min, $hour, $day, $mon, $year) = localtime(time);
	$now_date=join("-",($year+1900,$mon+1,$day));
	$now_time=join(":",($hour,$min,$sec));
	return ($now_date, $now_time);
}

