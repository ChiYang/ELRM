#!/usr/bin/perl
#Script FindRules.pl
#Author: Chi Yang
#Date: 2014/05/25
#License: MIT licence

$minConf = $ARGV[0]; #Minimum confidence; 0.6
$MCS = $ARGV[1]; #Minimum correlation strength; 0.1
$Input = $ARGV[2]; #The output file from Apriori
$output = $ARGV[3]; #Output for association rules
$output2 = $ARGV[4]; #Output for association features
#Example: perl FindRules.pl 0.6 0.1 ./example/Frequent_SingleBase_Sets/PurR.apriori_output ../example/AssociationRules/PurR.associationRules ../example/AssociationRules/PurR.associationFeatures
if(@ARGV != 5){
	print STDERR "Five arguments have to be set.\n";
	print STDERR "Usage: perl FindRules.pl [minimum confidence] [MCS] [Input file] [Output file1] [Output file2]\n";
	exit;
}elsif($minConf <= 0 or $minConf >= 1){
	print STDERR "Error: The range for the minimum confidence should be > 0 and < 1\n";
	exit;
}elsif($MCS <= 0 or $MCS >= 1){
	print STDERR "Error: The range for the minimum correlation strength should be > 0 and < 1\n";
	exit;
}elsif(!-e $Input){
	print STDERR "Error: The input file '$Input' does not exist.\n";
	exit;
}

%PAR = ();
%NAR = ();
%sup = ();
@X = ();
open(FILE, $Input);
while(<FILE>){
	chomp;
	@array = split(/\t/);
	$items = $array[0];
#	print $items."=>";
	@items = split(/\,/, $items);
	@items = sort sortFunc (@items);
	
	$items = join(",",@items);
#	print $items."\n";
	$sup{$items} = $array[2];
	$itemsetLength = scalar(@items);
	if($itemsetLength > 1){
		push(@X, $items);
	}
}
close(FILE);



%Group = ();

#@{$X[$i]} = sort{$a <=> $b}(@{$X[$i]});

foreach $itemsets (@X){	
	@itemsets = split(/\,/, $itemsets);
	#print join(",",@itemsets)."=>";
	%A2B = ();
	$setLength = @itemsets;
	for($j = 1; $j <= ($setLength/2); $j++){
		@comb = ordered_combinations([@itemsets], $j);
		foreach $c (@comb){
			@a = split(/\,/,$c);
			@a = &unique(@a);
			if($j == scalar(@a)){
				$key = join(',',sort sortFunc @a);
				@B = ();
				foreach $item (@itemsets){
					$inA = 0;
					foreach $a (@a){
						if($item eq $a){
							$inA = 1;
						}
					}
					if($inA == 0){
						push(@B, $item);
					}
				}
				$value = join(',',sort sortFunc (@B));
				$A2B{$key}= $value;
			}
		}
	}
	foreach $setA (keys(%A2B)){
		$setB = $A2B{$setA};
		$union = join(",",sort sortFunc (split(',', $setA.",".$setB)));
		if(!defined($sup{$union})){
			print $union."\n";
			exit;
		}
		$interest = $sup{$union} - ($sup{$setA}*$sup{$setB});
		$supAB = $sup{$union};
		#$sup_A = 1 - $sup{$setA};
		#$sup_B = 1 - $sup{$setB};
		#$sup_A_B = 1 - $sup{$setA} - $sup{$setB} + $supAB;
		#$supA_B = $sup{$setA} - $supAB;
		#$sup_AB = $sup{$setB} - $supAB;
		

		$correlation = $interest/(sqrt($sup{$setA}*(1-$sup{$setA})*$sup{$setB}*(1-$sup{$setB})));
		#The VARCC measure was implemented in the followings
		if(abs($correlation) > $MCS){ #Check if the correlation is larger than MCS
			#print "$setA=>$setB: \|$wsup{$union} - $wsup{$setA} * $wsup{$setB} \| = $interest\n";
			if($correlation > 0){
				#Check if the rule for A=>B is valid
				$confA2B = $sup{$union}/$sup{$setA};
				if($confA2B > $minConf){
					$PAR{$setA."=>".$setB} = $correlation."\t".$confA2B;
					$Group{$union} ++;
				}
				#Check if the rule for ~A=>~B is valid
				$conf_A2_B = (1-$sup{$setA}-$sup{$setB}+$sup{$union})/(1-$sup{$setA});
				if($conf_A2_B > $minConf){
					$notA = "~(".$setA.")";
					$notB = "~(".$setB.")";
					$NAR{$notA."=>".$notB} = $correlation."\t".$conf_A2_B;
					$Group{$union} ++;
				}
			}else{
				#Check if the rule for A=>~B is valid
				$confA2_B = 1 - ($sup{$union}/$sup{$setA});
				if($confA2_B > $minConf){
					$notB = "~(".$setB.")";
					$NAR{$setA."=>".$notB} = $correlation."\t".$confA2_B;
					$Group{$union} ++;
				}
				#Check if the rule for ~A=>B is valid
				$conf_A2B = ($sup{$setB} - $sup{$union})/(1-$sup{$setA});
				if($conf_A2B > $minConf){
					$notA = "~(".$setA.")";
					$NAR{$notA."=>".$setB} = $correlation."\t".$conf_A2B;
					$Group{$union} ++;
				}
			}
		}
	}
}
open(OUTPUT, ">$output");
	print OUTPUT "Association rule\tCorrelation Strength\tConfidence\n";
	foreach $key (keys(%PAR)){
		print OUTPUT $key."\t".$PAR{$key}."\n";
	}
	foreach $key (keys(%NAR)){
		print OUTPUT $key."\t".$NAR{$key}."\n";
	}
close(OUTPUT);

open(OUTPUT,">$output2");
@OutputKeys = sort sortFuncForOutput keys(%Group);
foreach $key (@OutputKeys){
	print OUTPUT "$key\n";
}
close(OUTPUT);

#@temp = ordered_combinations([qw(11a 12b c d)], 2);
#print join("\n", @temp)."\n";


sub ordered_combinations
{
  my ($data, $k) = @_;

  return if $k < 1;

  my $results = $data;

  while (--$k) {
    my @new;
    for my $letter (@$data) {
      push @new, map { $letter .','. $_ } @$results;
    } # end for $letter in @$data

    $results = \@new;
  } # end while --$k is not 0

  return @$results;
} # end ordered_combinations


sub unique{
	%output = ();
	foreach $x (@_){
		$output{$x}++;
	}
	#my @unique = grep { not $seen{$_} ++ } @_;
	return keys(%output);
}

sub sortFunc{
	%Temp = ();
	$Temp{'A'}=0;
	$Temp{'C'}=1;
	$Temp{'G'}=2;
	$Temp{'T'}=3;
	if($a =~ /([ATGC])(\d+)/){
		$nucA = $1;
		$posA = $2;
	}
	if($b =~ /([ATGC])(\d+)/){
		$nucB = $1;
		$posB = $2;
	}
	if($posA > $posB){
		return 1;
	}elsif($posA < $posB){
		return -1;
	}else{
		if($Temp{$nucA} > $Temp{$nucB}){
			return 1;
		}else{
			return -1;
		}
	}
}
sub sortFuncForOutput{
	%Temp = ();
	$Temp{'A'}=0;
	$Temp{'C'}=1;
	$Temp{'G'}=2;
	$Temp{'T'}=3;
	@a = split(/\,/, $a);
	@b = split(/\,/, $b);
	if(scalar(@a) > scalar(@b)){
		return 1;
	}elsif(scalar(@a) < scalar(@b)) {
		return -1;
	}
	for($_i = 0; $_i < @a; $_i++){
		$_a = $a[$_i];
		$_b = $b[$_i];
		if($_a =~ /([ATGC])(\d+)/){
			$nucA = $1;
			$posA = $2;
		}
		if($_b =~ /([ATGC])(\d+)/){
			$nucB = $1;
			$posB = $2;
		}
		if($posA > $posB){
			return 1;
		}elsif($posA < $posB){
			return -1;
		}else{
			if($Temp{$nucA} > $Temp{$nucB}){
				return 1;
			}elsif($Temp{$nucA} < $Temp{$nucB}){
				return -1;
			}else{
				next;
			}
		}
	}
	
}
