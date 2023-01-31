# This script parses the output file of BLAT:
	# When one specific read reaches more than one region aligned ($TGAPbases!=0), the length of each adjacent gap should be compared individually with the length of the corresponding intron depending on the gene transcript if canonical splicing

open FILE1, "<$ARGV[0]";	# $ARGV[0] = ../patients.nsubs.rna_seq.list
while (<FILE1>){
	chomp;
	$LINE1=$_;
	@PARAMETERS1=split(/\t/,$LINE1);
	$REQUEST=$PARAMETERS1[0];
	
	$HASH{$REQUEST}=1;
}
close FILE1;

open FILE2, "<$ARGV[1]";	# $ARGV[1] = ../DATABASE/exons.GRCh37.p13.txt
while (<FILE2>){
	chomp;
	$LINE2=$_;

	if ($LINE2=~/Gene name/){
		goto END;
	}

	@PARAMETERS2=split(/\t/,$LINE2);
	$SYMBOL=$PARAMETERS2[0];
	$TRANSCRIPT=$PARAMETERS2[2];
	$START=$PARAMETERS2[4];
	$END=$PARAMETERS2[5];
	$NUMBER=$PARAMETERS2[6];
	
	$EXON{$SYMBOL}{$TRANSCRIPT}{$NUMBER}[0]=$START;
	$EXON{$SYMBOL}{$TRANSCRIPT}{$NUMBER}[1]=$END;
	
	END:
}
close FILE2;

foreach $KEY_REQUEST (sort {$a<=>$b} keys %HASH){

	$REFERENCE="/path/to/" . "$KEY_REQUEST" . "/" . "donor.splice.sites/nucleotide.substitutions/candidates.mes";
	open FILE3, "<$REFERENCE";
	LINE:
	while (<FILE3>){
		chomp;
		
		if ($_=~/CHR/){
			goto LINE;
		}
		
		@HEADER=split(/\t/,$_);
		$CHR=$HEADER[0];
		$POS=$HEADER[1];
		$OBS=$HEADER[3];
		$NAME=$HEADER[4];
		$STRAND=$HEADER[5];
		
		if (exists $MUT{$NAME}{$POS}{$OBS}){	# Several 9-mer sequences for the same mutation could be often an indication of different cryptic splicing events
			goto LINE;
		}
		$MUT{$NAME}{$POS}{$OBS}=1;
		$PATH="/path/to/" . "$KEY_REQUEST" . "/" . "donor.splice.sites/nucleotide.substitutions" . "/" . "$NAME" . "." . "$POS" . "." . "$OBS";
		
		foreach $KEY1_TRANSCRIPT (sort {$a<=>$b} keys %{$EXON{$NAME}}){
			foreach $KEY1_NUMBER (sort {$a<=>$b} keys %{$EXON{$NAME}{$KEY1_TRANSCRIPT}}){
				$FLAG=0;
				
				$START_REC1A=$EXON{$NAME}{$KEY1_TRANSCRIPT}{$KEY1_NUMBER}[0];
				$END_REC1A=$EXON{$NAME}{$KEY1_TRANSCRIPT}{$KEY1_NUMBER}[1];
				$NUMBER_exonPOST=$KEY1_NUMBER+1;
				if (exists $EXON{$NAME}{$KEY1_TRANSCRIPT}{$NUMBER_exonPOST}){	# The mutation could be located in the last exon of the gene so there would be no later intron
					$START_REC1B=$EXON{$NAME}{$KEY1_TRANSCRIPT}{$NUMBER_exonPOST}[0];
					$END_REC1B=$EXON{$NAME}{$KEY1_TRANSCRIPT}{$NUMBER_exonPOST}[1];
					if ($START_REC1A<=$POS && $POS<=$END_REC1A){	# $START_REC1A and $END_REC1A would define the exon limits where mutation of interest is located
						if ($STRAND==1){
							$START_intronPOST=$END_REC1A+1;
							$END_intronPOST=$START_REC1B;
							
							$FLAG++;
						}
						else{
							$START_intronPOST=$END_REC1B+1;
							$END_intronPOST=$START_REC1A;
							
							$FLAG++;
						}
						
						if ($FLAG!=0){
							$LENGTH_INTRONpost=$END_intronPOST-$START_intronPOST;
							$INTRON{$LENGTH_INTRONpost}=1;
							
							$IN="YES";
						}
					}
					else{
						if ($STRAND==1){
							if ($END_REC1A<$POS && $POS<$START_REC1B){
								$START_intron=$END_REC1A+1;
								$END_intron=$START_REC1B;
								
								$FLAG++;
							}
						}
						else{
							if ($END_REC1B<$POS && $POS<$START_REC1A){
								$START_intron=$END_REC1B+1;
								$END_intron=$START_REC1A;
								
								$FLAG++;
							}
						}
						
						if ($FLAG!=0){
							$LENGTH_INTRON=$END_intron-$START_intron;
							$INTRON{$LENGTH_INTRON}=1;
							
							$IN="YES";
						}
					}
				}	
			}
		}
		
		$PRINT1="$PATH" . "/" . "reads.cryptic.splicing.v5.txt";
		open FILE_print1, ">$PRINT1";
		$PRINT2="$PATH" . "/" . "reads.alternative.splicing.v5.txt";
		open FILE_print2, ">$PRINT2";
		$PRINT3="$PATH" . "/" . "reads.unknown.v5.txt";
		open FILE_print3, ">$PRINT3";
		
		$FILE="$PATH" . "/" . "reads.blat.txt";
		
		$COUNT=0;
		open FILE4, "<$FILE";
		while (<FILE4>){
			chomp;
			$LINE4=$_;
			
			if ($COUNT<5){
				print FILE_print1 "$LINE4\n";
				print FILE_print2 "$LINE4\n";
				print FILE_print3 "$LINE4\n";
				$COUNT++;
				goto END;
			}
			
			if ($IN ne YES){
				print FILE_print3 "$LINE4\n";
				goto END;
			}
			
			@ALIGNMENT=split(/\t/,$LINE4);
			$MATCH=$ALIGNMENT[0];
			$TGAPbases=$ALIGNMENT[7];
			$CHRM=$ALIGNMENT[13];
			$Tstart=$ALIGNMENT[15];
			$Tend=$ALIGNMENT[16];
			# The length of the next arrays are always the same. When the value is equal to or greater than 2 (first position of an array is numbered as 0), inner splicing events are represented due to the presence of small exons
			$inBASES_pre=substr($ALIGNMENT[18],0,-1);
			@inBASES=split(/,/,$inBASES_pre);
			$inPOS_pre=substr($ALIGNMENT[20],0,-1);
			@inPOS=split(/,/,$inPOS_pre);
		
			if ($CHRM!=$CHR || $Tstart>$POS || $POS>$Tend){
				goto END;
			}
		
			if ($TGAPbases!=0){
				$TICKET=0;
				if ($#inBASES>=2){	# It is the same to say "$#inPOS>=2". The read reaches 3 or more regions (then there are minimum 2 inner gaps)
					# The lengths of all inner gaps will be stored. On the other hand, starting and ending coordinates will be stored for the region aligned where mutation is located
					for $A(0..$#inBASES){
						# The alignment starting positions would be equal to "$inPOS[$A]+1" whereas the lengths would result from "$inBASES[$A]-1". With this in mind:
						# $lastPOS=($inBASES[$A]-1)+($inPOS[$A]+1);						
						$lastPOS=$inBASES[$A]+$inPOS[$A];
						
						if ($A<$#inBASES){
							$B=$A+1;
							# As previously noted, "$inPOS[$B]+1" corresponds to the real start position of an region aligned so we have to carry out a correction (-1) to determine the length of the adyacent gap 
							# $inLENGTH=($inPOS[$B]+1)-$lastPOS-1;
							$inLENGTH=$inPOS[$B]-$lastPOS;
							if (($inPOS[$A]+1)<=$POS && $POS<=$lastPOS){
								if ($inLENGTH==0){
									print FILE_print3 "$LINE4\n";
									goto END;
								}
								# The first region aligned has a value for the variable $A equal to 0
								$exonCOORD{$A}=1;
								# The first intron has a value for the variable $A equal to 0								
								$gapLENGTH{$A}[0]=$inLENGTH;	# Lengths are storaged for the different inner gaps
							
								$TICKET++;
							}
							if (($lastPOS+1)<=$POS && $POS<=$inPOS[$B]){
								$gapLENGTH{$A}[0]=$inLENGTH;
							}
						}
						else{
							if (($inPOS[$A]+1)<=$POS && $POS<=$lastPOS){							
								$exonCOORD{$A}=1;
								
								$TICKET++;
							}
						}
					}
					if ($TICKET==0){	# The mutation is necessarily located in a gap
						foreach $KEY_A (sort {$a<=>$b} keys %gapLENGTH){
							$gapLENGTH_mut=$gapLENGTH{$KEY_A}[0];
							foreach $KEY_LENGTH_INTRON (sort {$a<=>$b} keys %INTRON){
								if ($KEY_LENGTH_INTRON==$gapLENGTH_mut){
									print FILE_print2 "$LINE4\n";
									goto END;
								}
							}
						}
					}
					else{
						foreach $KEY_A (sort {$a<=>$b} keys %exonCOORD){
							if ($KEY_A==0){	# The mutation is located in the first region aligned
								if ($STRAND==1){
									$gapLENGTH_post=$gapLENGTH{$KEY_A}[0];
									foreach $KEY_LENGTH_INTRON (sort {$a<=>$b} keys %INTRON){
										if ($KEY_LENGTH_INTRON==$gapLENGTH_post){
											print FILE_print2 "$LINE4\n";
											goto END;
										}
									}
								}
								else{
									print FILE_print3 "$LINE4\n";
									goto END;
								}
							}
							elsif ($KEY_A==$#inBASES){	# The mutation is located in the last region aligned
								if ($STRAND==-1){
									$KEY_B=$KEY_A-1;
									$gapLENGTH_post=$gapLENGTH{$KEY_B}[0];
									foreach $KEY_LENGTH_INTRON (sort {$a<=>$b} keys %INTRON){
										if ($KEY_LENGTH_INTRON==$gapLENGTH_post){
											print FILE_print2 "$LINE4\n";
											goto END;
										}
									}
								}
								else{
									print FILE_print3 "$LINE4\n";
									goto END;
								}
							}
							else{
								if ($STRAND==1){	
									$gapLENGTH_post=$gapLENGTH{$KEY_A}[0];
									foreach $KEY_LENGTH_INTRON (sort {$a<=>$b} keys %INTRON){
										if ($KEY_LENGTH_INTRON==$gapLENGTH_post){
											print FILE_print2 "$LINE4\n";
											goto END;
										}
									}
								}
								else{
									$KEY_B=$KEY_A-1;
									$gapLENGTH_post=$gapLENGTH{$KEY_B}[0];
									foreach $KEY_LENGTH_INTRON (sort {$a<=>$b} keys %INTRON){
										if ($KEY_LENGTH_INTRON==$gapLENGTH_post){
											print FILE_print2 "$LINE4\n";
											goto END;
										}
									}
								}
							}
						}
					}
				}		
				else{	# $#inbases (or $#inPOS) will inevitable be equal to 1. The read output of BLAT reaches 2 regions (there is 1 inner gap)
					# Starting and ending coordinates will be stored for the region aligned where mutation is located
					for $C(0..$#inBASES){
						# The alignment starting positions would be equal to "$inPOS[$C]+1" whereas the lengths would result from "$inBASES[$C]-1". With this in mind:
						# $lastPOS=($inBASES[$C]-1)+($inPOS[$C]+1);						
						$lastPOS=$inBASES[$C]+$inPOS[$C];
						
						if (($inPOS[$C]+1)<=$POS && $POS<=$lastPOS){
							# The first region aligned has a value for the variable $C equal to 0
							$exonCOORD{$C}=1;
							#$KEY="ZERO";
							#$exonCOORD{$KEY}=1;
							$TICKET++;
						}
					
					}
					if ($TICKET==0){	# The mutation is necessarily located in the gap
						foreach $KEY_LENGTH_INTRON (sort {$a<=>$b} keys %INTRON){
							if ($KEY_LENGTH_INTRON==$TGAPbases){
								print FILE_print2 "$LINE4\n";
								goto END;
							}
						}
					}			
					else{
						foreach $KEY_C (sort {$a<=>$b} keys %exonCOORD){
							if ($KEY_C==0){	# The mutation is located in the first region aligned
								if ($STRAND==1){	
									foreach $KEY_LENGTH_INTRON (sort {$a<=>$b} keys %INTRON){
										if ($KEY_LENGTH_INTRON==$TGAPbases){
											print FILE_print2 "$LINE4\n";
											goto END;
										}
									}
								}
								else{
									print FILE_print3 "$LINE4\n";
									goto END;
								}
							}
							else{ # $KEY_C=1; The mutation is necessarily located in the second and last region aligned
								if ($STRAND==-1){
									foreach $KEY_LENGTH_INTRON (sort {$a<=>$b} keys %INTRON){
										if ($KEY_LENGTH_INTRON==$TGAPbases){
											print FILE_print2 "$LINE4\n";
											goto END;
										}
									}
								}
								else{
									print FILE_print3 "$LINE4\n";
									goto END;
								}
							}
						}
					}	
				}
				print FILE_print1 "$LINE4\n";
			}
		
			if ($TGAPbases==0 && $MATCH<76){	# Reads to quantify the possible cryptic splicing event
				print FILE_print3 "$LINE4\n";
			}
			
			END:
			
			%exonCOORD=();
			%gapLENGTH=();
		}
		close FILE4;
		
		%INTRON=();
		
		close FILE_print1;
		close FILE_print2;
		close FILE_print3;
	}
	close FILE3;
	
	%MUT=();

}