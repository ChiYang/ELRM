#!/usr/bin/perl
#Script FindRules.pl
#Author: Chi Yang
#Date: 2015/07/15
#License: MIT licence
$modelFilePath = $ARGV[0]; #The ELRM models
$outputFilePath1 = $ARGV[1]; #Association links
$outputFilePath2 = $ARGV[2]; #Single Base Features


%AssociationFeatures = ();
%SingleBaseFeatures = ();
%CodedSingleBaseFeatures = ();
@ACGT = ("A","C","G","T");
open(FILE, $modelFilePath) || die ($!."\n");
while(<FILE>){
	chomp;
	($term, $coef) = split(/\t/);
	if($term ne "(Intercept)" && $term ne "s"){
		#$coef = sprintf("%.8f",$coef);
		@terms = split(/\:/, $term);
		@assocNucs = ();
		if(scalar(@terms) > 1){
			$AssociationFeatures{$term} = $coef;
		}elsif(scalar(@terms) == 1){
			if($term =~ /^Z(\d+)$/){
				$pos = int(($1-1) / 3)+1;
				$temp = ($1-1) % 3; #
				$nuc = $ACGT[$temp];
				$decodedTerm = $nuc.$pos;
				$SingleBaseFeatures{$decodedTerm} = $coef;
				#$CodedSingleBaseFeatures{"Z".$1} = $coef;
				if($coef != 0){
					push(@selectedCodedSF, "Z".$1);
				}
				if($nuc eq "G"){ #1,2,0, 1, 2, 0
					$decodedTerm = "T".$pos;
					$SingleBaseFeatures{$decodedTerm} = -1*($SingleBaseFeatures{"A".$pos} + $SingleBaseFeatures{"C".$pos} + $SingleBaseFeatures{"G".$pos});
				}
			}
		}
	}elsif($term eq "(Intercept)"){
		#$coef = sprintf("%.8f",$coef);
		$intercept = $coef;
		#print "$term\n";
	}elsif($term eq "s"){
		$theS = $coef;
	}
}
close(FILE);



open(OUTPUT, ">$outputFilePath1");
print OUTPUT "AssociationFeatures\tSingleBaseFeature1\tSingleBaseFeature2\tCoefficients\n";
foreach $af (sort{sortFuncForOutput}(%AssociationFeatures)){
	if($AssociationFeatures{$af} != 0){
		@sf = split(/\:/, $af);
		$length  =@sf;
		if($length > 2){
			for($i = 0; $i < $length-1; $i++){
				print OUTPUT $af."\t". $sf[$i]."\t".$sf[$i+1]."\t".$AssociationFeatures{$af}."\n";
			}

		}
		print OUTPUT $af."\t".$sf[0]."\t".$sf[$length-1]."\t".$AssociationFeatures{$af}."\n";
	}
}
close(OUTPUT);

exit;
open(OUTPUT, ">$outputFilePath2");
print OUTPUT "SingleBaseFeatures\tCoefficients\n";
foreach $sf (sort {sortFuncForOutput} (%SingleBaseFeatures)){
	if($SingleBaseFeatures{$sf} != 0){
		print OUTPUT "$sf\t$SingleBaseFeatures{$sf}\n";
	}
}
close(OUTPUT);

close(OUTPUT);


sub sortFuncForOutput{
	%Temp = ();
	$Temp{'A'}=0;
	$Temp{'C'}=1;
	$Temp{'G'}=2;
	$Temp{'T'}=3;
	@a = split(/\:/, $a);
	@b = split(/\:/, $b);
	if(scalar(@a) > scalar(@b)){
		return 1;
	}elsif(scalar(@a) < scalar(@b)) {
		return -1;
	}else{
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
}