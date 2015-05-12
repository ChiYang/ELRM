#!/usr/bin/perl
#Script SequenceCoding.pl
#Author: Chi Yang
#Date: 2014/06/02
#License: MIT licence
$motifLength = $ARGV[0];
$inputFastaFile = $ARGV[1];
$interactionFile = $ARGV[2]; #Here the interaction file is the file containing association features
$outputCodeFile = $ARGV[3]; 
$isTraining = $ARGV[4];

#Coding sequences for training
#perl SequenceCoding.pl 16 ../example/TrainSeq/PurR_training.fasta ../example/AssociationRules/PurR.associationFeatures ../example/TrainSeq/PurR_training.dummyCodes 1

#Coding sequences for testing
#perl SequenceCoding.pl 16 ../example/TestSeq/PurR_test.fasta ../example/AssociationRules/PurR.associationFeatures ../example/TestSeq/PurR_testing.dummyCodes 0


if(@ARGV != 5){
	print STDERR "Five arguments have to be set.\n";
	print STDERR "Usage: perl SequenceCoding.pl [motif length] [Input file] [AssociationFeature file] [Output file1] [Is used for training]\n";
	exit;
}

if(!$motifLength){
	print STDERR "Error: Motif length is not specified.\n";
	exit 1;
}

if(!-e $inputFastaFile ){
	print STDERR "Error: Input fasta file is not existed.\n";
	exit 1;
}
if(!-e $interactionFile){
	print STDERR "Error: The association feature file is not existed";
	exit 1;
}


$Z = 1;
@terms = ();
for($i = 0; $i < $motifLength;$i++){
	$term1 = "Z".$Z;
	$Z = $Z + 1;
	$term2 = "Z".$Z;
	$Z = $Z + 1;
	$term3 = "Z".$Z;
	$Z = $Z + 1;
	push(@terms, $term1, $term2, $term3);
}
@asTerms = ();
open(FILE, "$interactionFile");
while(<FILE>){
	chomp;
	@array = split(/\t/);
	@t = split(/\,/, $array[0]);
	push(@asTerms, join(":", @t));
}
close(FILE);
open(OUTPUT, ">$outputCodeFile");
if(!$isTraining){
	print OUTPUT "Title\tStart\tEnd\tStrand\tSequence\t".join("\t", @terms);
	if(@asTerms > 0){
		print OUTPUT "\t". join("\t", @asTerms)."\n"
	}else{
		print OUTPUT "\n";
	}
}else{
	print OUTPUT "y\t".join("\t", @terms);
	if(@asTerms > 0){
		print OUTPUT "\t". join("\t", @asTerms)."\n";
	}else{
		print OUTPUT "\n";
	}
}

#Reading sequences and output the results
@testCodes = ();
%testFeatuers = ();
@testHeads = ();
$counter =0;
$title = "";
open(FILE, $inputFastaFile) || die("No such file!\n");
while(<FILE>){
	chomp;
	$_ =~ s/\r//g;
	if($_ =~ /^>(.*)/){
		if($title ne ""){
			#print $title."=>".length($sequence)."\n";
			for($i = 0; $i < (length($sequence) - $motifLength+ 1); $i++){
				#print "Processing $i / ".length($sequence)."\r";
				$subS = substr($sequence, $i, $motifLength);
				$rsubS = reverse($subS);
				$rsubS =~ tr/[ATCG]/[TAGC]/;
				if(!$isTraining){
					@testCode = &dummyCoding($subS);
					%testFeatures = &extractFeatures($subS);
					$Output = "$title\t".($i+1)."\t".($i+$motifLength)."\t+\t$subS";
					$Output .= "\t".join("\t", @testCode);
					@associationValues = ();
					for($k = 0; $k < @asTerms;$k ++){
						$value = 1;
						@t = split(/\:/, $asTerms[$k]);
						for($l = 0 ; $l < @t; $l++){
							if(defined($testFeatures{$t[$l]})){
								$value = $value * 1;
							}else{
								$value = 0;
							}
						}
						push(@associationValues, $value);
					}
					if(@associationValues > 0){
						$Output .= "\t".join("\t", @associationValues)."\n";
					}else{
						$Output .= "\n";
					}
					@rtestCode = &dummyCoding($rsubS);
					%rtestFeatures = &extractFeatures($rsubS);
					$Output .= "$title\t".($i+1)."\t".($i+$motifLength)."\t-\t$rsubS";
					$Output .= "\t".join("\t", @rtestCode);
					@associationValues = ();
					$value = 1;
					for($k = 0; $k < @asTerms;$k ++){
						$value = 1;
						@t = split(/\:/, $asTerms[$k]);
						for($l = 0 ; $l < @t; $l++){
							if(defined($rtestFeatures{$t[$l]})){
								$value = $value *1;
							}else{
								$value = 0;
							}
						}
						push(@associationValues, $value);
					}
					if(@associationValues > 0){
						$Output .= "\t".join("\t", @associationValues)."\n";
					}else{
						$Output .= "\n";
					}
					print OUTPUT $Output;
				}else{
					@testCode = &dummyCoding($subS);
					%testFeatures = &extractFeatures($subS);
					if($title =~ /positive/i){
						$Output = "1\t".join("\t", @testCode);
					}else{
						$Output = "0\t".join("\t", @testCode);
					}
					@associationValues = ();
					for($k = 0; $k < @asTerms;$k ++){
						$value = 1;
						@t = split(/\:/, $asTerms[$k]);
						for($l = 0 ; $l < @t; $l++){
							if(defined($testFeatures{$t[$l]})){
								$value = $value * 1;
							}else{
								$value = 0;
							}
						}
						push(@associationValues, $value);
					}
					if(@associationValues > 0){
						$Output .= "\t".join("\t", @associationValues)."\n";
					}else{
						$Output .= "\n";
					}
					print OUTPUT $Output;
				}
			}
		}
		$temp = $1;
		if($temp =~ /(.*?)\s/){
			$title = $1;
		}else{
			$title = $temp;
		}
		$sequence = "";
		next;
	}
	$sequence .= $_;
}
close(FILE);

if($title){
	#print $title."=>".length($sequence)."\n";
	for($i = 0; $i < (length($sequence) - $motifLength+ 1); $i++){
		#print "Processing $i / ".length($sequence)."\r";
		$subS = substr($sequence, $i, $motifLength);
		$rsubS = reverse($subS);
		$rsubS =~ tr/[ATCG]/[TAGC]/;
		if(!$isTraining){
			@testCode = &dummyCoding($subS);
			%testFeatures = &extractFeatures($subS);
			$Output = "$title\t".($i+1)."\t".($i+$motifLength)."\t+\t$subS";
			$Output .= "\t".join("\t", @testCode);
			@associationValues = ();
			for($k = 0; $k < @asTerms;$k ++){
				$value = 1;
				@t = split(/\:/, $asTerms[$k]);
				for($l = 0 ; $l < @t; $l++){
					if(defined($testFeatures{$t[$l]})){
						$value = $value * 1;
					}else{
						$value = 0;
					}
				}
				push(@associationValues, $value);
			}
			if(@associationValues > 0){
				$Output .= "\t".join("\t", @associationValues)."\n";
			}else{
				$Output .= "\n";
			}
			@rtestCode = &dummyCoding($rsubS);
			%rtestFeatures = &extractFeatures($rsubS);
			$Output .= "$title\t".($i+1)."\t".($i+$motifLength)."\t-\t$rsubS";
			$Output .= "\t".join("\t", @rtestCode);
			@associationValues = ();
			$value = 1;
			for($k = 0; $k < @asTerms;$k ++){
				$value = 1;
				@t = split(/\:/, $asTerms[$k]);
				for($l = 0 ; $l < @t; $l++){
					if(defined($rtestFeatures{$t[$l]})){
						$value = $value *1;
					}else{
						$value = 0;
					}
				}
				push(@associationValues, $value);
			}
			if(@associationValues > 0){
				$Output .= "\t".join("\t", @associationValues)."\n";
			}else{
				$Output .= "\n";
			}
			
			print OUTPUT $Output;
		}else{
			@testCode = &dummyCoding($subS);
			%testFeatures = &extractFeatures($subS);
			if($title =~ /positive/i){
				$Output = "1\t".join("\t", @testCode);
			}else{
				$Output = "0\t".join("\t", @testCode);
			}
			@associationValues = ();
			for($k = 0; $k < @asTerms;$k ++){
				$value = 1;
				@t = split(/\:/, $asTerms[$k]);
				for($l = 0 ; $l < @t; $l++){
					if(defined($testFeatures{$t[$l]})){
						$value = $value * 1;
					}else{
						$value = 0;
					}
				}
				push(@associationValues, $value);
			}
			if(@associationValues > 0){
				$Output .= "\t".join("\t", @associationValues)."\n";
			}else{
				$Output .= "\n";
			}
			print OUTPUT $Output;
		}
	}
}
sub dummyCoding{
	my $s, $seq = @_[0];
	my @s = split("", $seq);
	my @codes = ();
	for($_x = 0; $_x < @s; $_x ++){
		if($s[$_x] eq "A"){
			push(@codes, 1,0,0);
		}elsif($s[$_x] eq "C"){
			push(@codes, 0,1,0);
		}elsif($s[$_x] eq "G"){
			push(@codes, 0,0,1);
		}elsif($s[$_x] eq "T"){
			push(@codes, -1,-1,-1);
		}
	}
	return @codes;
}
sub extractFeatures{
	my $s, $seq = @_[0];
	my @s = split("", $seq);
	my %features = ();
	for($_x = 0; $_x < @s; $_x ++){
		$nuc = $s[$_x];
		$pos = $_x + 1;
		$f = $nuc.$pos;
		$features{$f}++;
	}
	return %features;
}
