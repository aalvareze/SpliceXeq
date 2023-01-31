# Distinction in splice sites
	# The script labels reads with recurrent splice sites and non-homologous reads in this way

# Labeling guidelines:
	# NOTES,	
		# CSS eq Cryptic Splice Site(s)
	# - Reads with one or two CSS possibly related to the mutation will be labeled as 'css-rep:$VALUE', with $VALUE being equal to the genomic position(s) of the recurrent CSS
	# - Reads with no identifiable CSS but candidates for exposing exon skipping events will be labeled as 'ESkip:$VALUE', storing $VALUE the genomic positions of the two reiterative canonical splice sites involved in the process
	# - Reads aligned without the introduction of any gaps but candidates for extending further than canonical splice sites will be labeled as 'IRet:$VALE', being $VALUE equal to the genomic position(s) of the canonical splice site(s)
	# - NA refers to 'Not Available' and will be applied in cases of non-recurrence
	
open FILE1, "<$ARGV[0]";	# $ARGV[0] = ../DATABASE/exons.GRCh37.p13.txt
while (<FILE1>){
	chomp;
	
	if ($_=~/Gene name/){
		next;
	}

	@LINE1=split(/\t/,$_);
	$SYMBOL=$LINE1[0];
	$TRANSCRIPT=$LINE1[2];
	$START=$LINE1[4];
	$END=$LINE1[5];
	$NUMBER=$LINE1[6];
	
	# To be fully optimised, the information accessability will depend on the type and form of data being stored
	
	$EXON1{$SYMBOL}{$START}=1;
	$EXON1{$SYMBOL}{$END}=1;
	
	$EXON2{$SYMBOL}{$TRANSCRIPT}{$NUMBER}[0]=$START;
	$EXON2{$SYMBOL}{$TRANSCRIPT}{$NUMBER}[1]=$END;
}
close FILE1;

open FILE2, "<$ARGV[1]";        # $ARGV[1] = ../patients.nsubs.rna_seq.list
while (<FILE2>){
	chomp;
	@LINE2=split(/\t/,$_);
	$REQUEST=$LINE2[0];
	
	$REFERENCE="/path/to/" . "$REQUEST" . "/" . "acceptor.splice.sites/nucleotide.substitutions/candidates.mes";
	open FILE3, "<$REFERENCE";
	LINE:
	while (<FILE3>){
		chomp;
		
		if ($_=~/CHR/){
			goto LINE;
		}
		
		@HEADER1=split(/\t/,$_);
		$POS=$HEADER1[1];
		$OBS=$HEADER1[3];
		$NAME=$HEADER1[4];
		$STRAND=$HEADER1[5];
		
		if (exists $MUT{$NAME}{$POS}{$OBS}){	# Several 23-mer sequences for the same mutation could be often an indication of different cryptic splicing events
			goto LINE;
		}
		$MUT{$NAME}{$POS}{$OBS}=1;
		
		$PATH_PRINT_PRE="/path/to/" . "$REQUEST" . "/" . "acceptor.splice.sites/nucleotide.substitutions" . "/" . "$NAME" . "." . "$POS" . "." . "$OBS" . "/" . "reads.cryptic.splicing.v5.nodup.nomin_optB_PRE.txt";
		open FILE_PRINT, ">$PATH_PRINT_PRE";
		
		$PATH_REF="/path/to/" . "$REQUEST" . "/" . "acceptor.splice.sites/nucleotide.substitutions" . "/" . "$NAME" . "." . "$POS" . "." . "$OBS" . "/" . "reads.cryptic.splicing.v5.nodup.txt";		
		$COUNT=0;
		open FILE_REF, "<$PATH_REF";
		while (<FILE_REF>){
			chomp;

			if ($COUNT<5){
				print FILE_PRINT "$_\n";
				$COUNT++;
				next;
			}
			
			@HEADER2=split(/\t/,$_);
			$MATCH=$HEADER2[0];
			$TGAPbases=$HEADER2[7];
			$ID=$HEADER2[9];
			$inBASES_pre=substr($HEADER2[18],0,-1);
			@inBASES=split(/,/,$inBASES_pre);
			$inPOS_pre=substr($HEADER2[20],0,-1);
			@inPOS=split(/,/,$inPOS_pre);
			
			if ($TGAPbases!=0){
				
				$REITER_POS1=0;
				$REITER_POS2=0;
				for $A(0..$#inBASES){
					if (($inPOS[$A]+1)<=$POS && $POS<=($inBASES[$A]+$inPOS[$A])){
						if ($A==0){
							$CHECK_POS2=$inBASES[$A]+$inPOS[$A];
							if (not exists $EXON1{$NAME}{$CHECK_POS2}){
								$REITER_POS2=$CHECK_POS2;
							}
						}
						elsif ($A==$#inBASES){
							$CHECK_POS1=$inPOS[$A]+1;
							if (not exists $EXON1{$NAME}{$CHECK_POS1}){
								$REITER_POS1=$CHECK_POS1;
							}
						}
						else{
							$CHECK_POS1=$inPOS[$A]+1;
							if (not exists $EXON1{$NAME}{$CHECK_POS1}){
								$REITER_POS1=$CHECK_POS1;
							}
							$CHECK_POS2=$inBASES[$A]+$inPOS[$A];
							if (not exists $EXON1{$NAME}{$CHECK_POS2}){
								$REITER_POS2=$CHECK_POS2;
							}
						}
					}
					else{	# The mutation could be located in a gap
						if ($A>0){
							$B=$A-1;
							if ((($inBASES[$B]+$inPOS[$B])+1)<=$POS && $POS<=$inPOS[$A]){
								$CHECK_POS1=$inPOS[$A]+1;
								if (not exists $EXON1{$NAME}{$CHECK_POS1}){
									$REITER_POS1=$CHECK_POS1;
								}
								$CHECK_POS2=$inBASES[$B]+$inPOS[$B];
								if (not exists $EXON1{$NAME}{$CHECK_POS2}){
									$REITER_POS2=$CHECK_POS2;
								}
							}
						}
					}
				}
				if ($REITER_POS1!=0){
					if (exists $HASH1{$REITER_POS1}){
						# If a read initially stored in %HASH1 comes to have a reiterative CSS, it will be moved to %HASH2 together with the homologous reads. Therefore, %HASH2 will only store reads with identical CSS
						$HASH2{$REITER_POS1}[0]=2;	# Number two comes from the addition of one read first stored in %HASH1 plus one new read with the same CSS
						$ID_REC1=$HASH1{$REITER_POS1}[1];
						$HASH2{$REITER_POS1}[1]=$ID_REC1;
						delete $HASH1{$REITER_POS1};
						
						$HASH2{$REITER_POS1}[2]=$ID;
					}
					elsif (exists $HASH2{$REITER_POS1}){
						$HASH2{$REITER_POS1}[0]=$HASH2{$REITER_POS1}[0]+1;
						$LENGTH1=@{$HASH2{$REITER_POS1}};
						$HASH2{$REITER_POS1}[$LENGTH1]=$ID;
					}
					else{
						$HASH1{$REITER_POS1}[0]=1;
						$HASH1{$REITER_POS1}[1]=$ID;
					}
				}
				if ($REITER_POS2!=0){
					if (exists $HASH1{$REITER_POS2}){
						# If a read initially stored in %HASH1 comes to have a reiterative CSS, it will be moved to %HASH2 together with the homologous reads. Therefore, %HASH2 will only store reads with identical CSS
						$HASH2{$REITER_POS2}[0]=2;	# Number two comes from the addition of one read first stored in %HASH1 plus one new read with the same CSS
						$ID_REC2=$HASH1{$REITER_POS2}[1];
						$HASH2{$REITER_POS2}[1]=$ID_REC2;
						delete $HASH1{$REITER_POS2};
						
						$HASH2{$REITER_POS2}[2]=$ID;
					}
					elsif (exists $HASH2{$REITER_POS2}){
						$HASH2{$REITER_POS2}[0]=$HASH2{$REITER_POS2}[0]+1;
						$LENGTH2=@{$HASH2{$REITER_POS2}};
						$HASH2{$REITER_POS2}[$LENGTH2]=$ID;
					}
					else{
						$HASH1{$REITER_POS2}[0]=1;
						$HASH1{$REITER_POS2}[1]=$ID;
					}
				}
				
				if ($REITER_POS1==0 && $REITER_POS2==0){	# Recovering of reads with no identifiable CSS but candidates for laying bare exon skipping events
					$CHECK_POS_RIGHT=0;
					$CHECK_POS_LEFT=0;
					for $C(0..$#inBASES){
						if (($inPOS[$C]+1)<=$POS && $POS<=($inBASES[$C]+$inPOS[$C])){
							if ($C==0){
								$CHECK_POS4=$inBASES[$C]+$inPOS[$C];
								$E=$C+1;
								$CHECK_POS3=$inPOS[$E]+1;
								$CHECK_POS_RIGHT="$CHECK_POS4" . "-" . "$CHECK_POS3";
							}
							elsif ($C==$#inBASES){
								$D=$C-1;
								$CHECK_POS4=$inBASES[$D]+$inPOS[$D];
								$CHECK_POS3=$inPOS[$C]+1;
								$CHECK_POS_LEFT="$CHECK_POS4" . "-" . "$CHECK_POS3";
							}
							else{
								$D=$C-1;
								$E=$C+1;
								foreach $KEY_TRANSCRIPT (sort {$a<=>$b} keys %{$EXON2{$NAME}}){
									for $F($D..$E){	
										foreach $KEY_NUMBER (sort {$a<=>$b} keys %{$EXON2{$NAME}{$KEY_TRANSCRIPT}}){
											$START_REC=$EXON2{$NAME}{$KEY_TRANSCRIPT}{$KEY_NUMBER}[0];
											$END_REC=$EXON2{$NAME}{$KEY_TRANSCRIPT}{$KEY_NUMBER}[1];
											if ($START_REC<=($inPOS[$F]+1) && ($inBASES[$F]+$inPOS[$F])<=$END_REC){
												if ($F==$D){
													$RANK1=$KEY_NUMBER;
													goto END;
												}
												elsif ($F==$C){
													$RANK2=$KEY_NUMBER;
													goto END;
												}
												else{
													$RANK3=$KEY_NUMBER;
													goto END;
												}
											}
										}
										END:
									}
									if (abs($RANK2-$RANK1)>1){
										$CHECK_POS4=$inBASES[$D]+$inPOS[$D];
										$CHECK_POS3=$inPOS[$C]+1;
										$CHECK_POS_LEFT="$CHECK_POS4" . "-" . "$CHECK_POS3";
									}
									if (abs($RANK2-$RANK3)>1){
										$CHECK_POS4=$inBASES[$C]+$inPOS[$C];
										$CHECK_POS3=$inPOS[$E]+1;
										$CHECK_POS_RIGHT="$CHECK_POS4" . "-" . "$CHECK_POS3";
									}
								}
							}
						}
						else{	# The mutation could be located in a gap
							if ($C>0){
								$D=$C-1;
								if ((($inBASES[$D]+$inPOS[$D])+1)<=$POS && $POS<=$inPOS[$C]){
									$CHECK_POS4=$inBASES[$D]+$inPOS[$D];
									$CHECK_POS3=$inPOS[$C]+1;
									$CHECK_POS_LEFT="$CHECK_POS4" . "-" . "$CHECK_POS3";
								}
							}
						}
					}
					if ($CHECK_POS_RIGHT=~/-/){
						if (exists $HASH1{$CHECK_POS_RIGHT}){
							# If a read initially stored in %HASH1 comes to have a reiterative ESkip event, it will be moved to %HASH2 together with the homologous reads. Therefore, %HASH2 will only store reads with identical ESkip events
							$HASH2{$CHECK_POS_RIGHT}[0]=2;	# Number two comes from the addition of one read first stored in %HASH1 plus one new read with the same ESkip event
							$ID_REC1=$HASH1{$CHECK_POS_RIGHT}[1];
							$HASH2{$CHECK_POS_RIGHT}[1]=$ID_REC1;
							delete $HASH1{$CHECK_POS_RIGHT};
							
							$HASH2{$CHECK_POS_RIGHT}[2]=$ID;
						}
						elsif (exists $HASH2{$CHECK_POS_RIGHT}){
							$HASH2{$CHECK_POS_RIGHT}[0]=$HASH2{$CHECK_POS_RIGHT}[0]+1;
							$LENGTH1=@{$HASH2{$CHECK_POS_RIGHT}};
							$HASH2{$CHECK_POS_RIGHT}[$LENGTH1]=$ID;
						}
						else{
							$HASH1{$CHECK_POS_RIGHT}[0]=1;
							$HASH1{$CHECK_POS_RIGHT}[1]=$ID;
						}
					}
					if ($CHECK_POS_LEFT=~/-/){
						if (exists $HASH1{$CHECK_POS_LEFT}){
							# If a read initially stored in %HASH1 comes to have a reiterative ESkip event, it will be moved to %HASH2 together with the homologous reads. Therefore, %HASH2 will only store reads with identical ESkip events
							$HASH2{$CHECK_POS_LEFT}[0]=2;	# Number two comes from the addition of one read first stored in %HASH1 plus one new read with the same ESkip event
							$ID_REC2=$HASH1{$CHECK_POS_LEFT}[1];
							$HASH2{$CHECK_POS_LEFT}[1]=$ID_REC2;
							delete $HASH1{$CHECK_POS_LEFT};
							
							$HASH2{$CHECK_POS_LEFT}[2]=$ID;
						}
						elsif (exists $HASH2{$CHECK_POS_LEFT}){
							$HASH2{$CHECK_POS_LEFT}[0]=$HASH2{$CHECK_POS_LEFT}[0]+1;
							$LENGTH2=@{$HASH2{$CHECK_POS_LEFT}};
							$HASH2{$CHECK_POS_LEFT}[$LENGTH2]=$ID;
						}
						else{
							$HASH1{$CHECK_POS_LEFT}[0]=1;
							$HASH1{$CHECK_POS_LEFT}[1]=$ID;
						}
					}
				}
				
			}
			else{	# Recovering of reads which overhang canonical splice sites
				$REITER_POS1=0;
				$REITER_POS2=0;
				$MATCH_HALF=int(($MATCH/2)+0.5);
				$MATCH_THIRD=int(($MATCH/3)+0.5);
				
				$CHECK_POS1=$inPOS[0]+1;
				$CHECK_POS2=$inBASES[$#inBASES]+$inPOS[$#inPOS];
				foreach $KEY_TRANSCRIPT (sort {$a<=>$b} keys %{$EXON2{$NAME}}){
					foreach $KEY_NUMBER (sort {$a<=>$b} keys %{$EXON2{$NAME}{$KEY_TRANSCRIPT}}){
						$START_REC=$EXON2{$NAME}{$KEY_TRANSCRIPT}{$KEY_NUMBER}[0];
						$END_REC=$EXON2{$NAME}{$KEY_TRANSCRIPT}{$KEY_NUMBER}[1];
						if ($START_REC<=$CHECK_POS1 && $CHECK_POS1<=$END_REC && $END_REC<$CHECK_POS2){
							$REITER_POS2="#" . "$END_REC";
						}
						if ($START_REC<=$CHECK_POS2 && $CHECK_POS2<=$END_REC && $CHECK_POS1<$START_REC){
							$REITER_POS1="#" . "$START_REC";
						}
						if (($START_REC<=($CHECK_POS1+$MATCH_HALF) && ($CHECK_POS1+$MATCH_HALF)<=$END_REC && $END_REC<$CHECK_POS2 && $CHECK_POS1<$START_REC) || ($START_REC<=($CHECK_POS1+$MATCH_THIRD) && ($CHECK_POS1+$MATCH_THIRD)<=$END_REC && $END_REC<$CHECK_POS2 && $CHECK_POS1<$START_REC) || ($START_REC<=($CHECK_POS2-$MATCH_THIRD) && ($CHECK_POS2-$MATCH_THIRD)<=$END_REC && $END_REC<$CHECK_POS2 && $CHECK_POS1<$START_REC)){
							$REITER_POS1="#" . "$START_REC";
							$REITER_POS2="#" . "$END_REC";
						}
					}
				}
				if ($REITER_POS1=~/#/){
					if (exists $HASH1{$REITER_POS1}){
						# If a read initially stored in %HASH1 comes to have a reiterative CSS, it will be moved to %HASH2 together with the homologous reads. Therefore, %HASH2 will only store reads with identical CSS
						$HASH2{$REITER_POS1}[0]=2;	# Number two comes from the addition of one read first stored in %HASH1 plus one new read with the same CSS
						$ID_REC1=$HASH1{$REITER_POS1}[1];
						$HASH2{$REITER_POS1}[1]=$ID_REC1;
						delete $HASH1{$REITER_POS1};
						
						$HASH2{$REITER_POS1}[2]=$ID;
					}
					elsif (exists $HASH2{$REITER_POS1}){
						$HASH2{$REITER_POS1}[0]=$HASH2{$REITER_POS1}[0]+1;
						$LENGTH1=@{$HASH2{$REITER_POS1}};
						$HASH2{$REITER_POS1}[$LENGTH1]=$ID;
					}
					else{
						$HASH1{$REITER_POS1}[0]=1;
						$HASH1{$REITER_POS1}[1]=$ID;
					}
				}
				if ($REITER_POS2=~/#/){
					if (exists $HASH1{$REITER_POS2}){
						# If a read initially stored in %HASH1 comes to have a reiterative CSS, it will be moved to %HASH2 together with the homologous reads. Therefore, %HASH2 will only store reads with identical CSS
						$HASH2{$REITER_POS2}[0]=2;	# Number two comes from the addition of one read first stored in %HASH1 plus one new read with the same CSS
						$ID_REC2=$HASH1{$REITER_POS2}[1];
						$HASH2{$REITER_POS2}[1]=$ID_REC2;
						delete $HASH1{$REITER_POS2};
						
						$HASH2{$REITER_POS2}[2]=$ID;
					}
					elsif (exists $HASH2{$REITER_POS2}){
						$HASH2{$REITER_POS2}[0]=$HASH2{$REITER_POS2}[0]+1;
						$LENGTH2=@{$HASH2{$REITER_POS2}};
						$HASH2{$REITER_POS2}[$LENGTH2]=$ID;
					}
					else{
						$HASH1{$REITER_POS2}[0]=1;
						$HASH1{$REITER_POS2}[1]=$ID;
					}
				}
			}
		}
		close FILE_REF;
		
		foreach $CRYPTIC_SITE2 (sort {$a<=>$b} keys %HASH2){
			$LENGTH_REC2=$HASH2{$CRYPTIC_SITE2}[0];
			# $LENGTH_REC2_ALT=@{$HASH2{$CRYPTIC_SITE2}};	$LENGTH_REC2_ALT=$LENGTH_REC2+1;
			# if ($LENGTH_REC2>=3){
				for $G(1..$LENGTH_REC2){
					$ID_REC2=$HASH2{$CRYPTIC_SITE2}[$G];
					$CSS=0;
					$STOREHOUSE[$CSS]=$CRYPTIC_SITE2;
					foreach $CRYPTIC_SITE2_BIS (sort {$a<=>$b} keys %HASH2){
						if ($CRYPTIC_SITE2!~/-/ && $CRYPTIC_SITE2!~/#/ && $CRYPTIC_SITE2_BIS!~/-/ && $CRYPTIC_SITE2_BIS!~/#/){
							if ($CRYPTIC_SITE2!=$CRYPTIC_SITE2_BIS){
								$LENGTH_REC2_BIS=$HASH2{$CRYPTIC_SITE2_BIS}[0];
								for $H(1..$LENGTH_REC2_BIS){
									$ID_REC2_BIS=$HASH2{$CRYPTIC_SITE2_BIS}[$H];
									if ($ID_REC2 eq $ID_REC2_BIS){
										$CSS++;
										$STOREHOUSE[$CSS]=$CRYPTIC_SITE2_BIS;
										splice @{$HASH2{$CRYPTIC_SITE2_BIS}},$H,1;
										$HASH2{$CRYPTIC_SITE2_BIS}[0]=$HASH2{$CRYPTIC_SITE2_BIS}[0]-1;
										if ($HASH2{$CRYPTIC_SITE2_BIS}[0]==0){
											delete $HASH2{$CRYPTIC_SITE2_BIS};
										}
									}
								}
							}
						}
						if ($CRYPTIC_SITE2=~/-/ && $CRYPTIC_SITE2_BIS=~/-/){
							if ($CRYPTIC_SITE2 ne $CRYPTIC_SITE2_BIS){
								$LENGTH_REC2_BIS=$HASH2{$CRYPTIC_SITE2_BIS}[0];
								for $H(1..$LENGTH_REC2_BIS){
									$ID_REC2_BIS=$HASH2{$CRYPTIC_SITE2_BIS}[$H];
									if ($ID_REC2 eq $ID_REC2_BIS){
										$CSS++;
										$STOREHOUSE[$CSS]=$CRYPTIC_SITE2_BIS;
										splice @{$HASH2{$CRYPTIC_SITE2_BIS}},$H,1;
										$HASH2{$CRYPTIC_SITE2_BIS}[0]=$HASH2{$CRYPTIC_SITE2_BIS}[0]-1;
										if ($HASH2{$CRYPTIC_SITE2_BIS}[0]==0){
											delete $HASH2{$CRYPTIC_SITE2_BIS};
										}
									}
								}
							}
						}
						if ($CRYPTIC_SITE2=~/#/ && $CRYPTIC_SITE2_BIS=~/#/){
							if ($CRYPTIC_SITE2 ne $CRYPTIC_SITE2_BIS){
								$LENGTH_REC2_BIS=$HASH2{$CRYPTIC_SITE2_BIS}[0];
								for $H(1..$LENGTH_REC2_BIS){
									$ID_REC2_BIS=$HASH2{$CRYPTIC_SITE2_BIS}[$H];
									if ($ID_REC2 eq $ID_REC2_BIS){
										$CSS++;
										$STOREHOUSE[$CSS]=$CRYPTIC_SITE2_BIS;
										splice @{$HASH2{$CRYPTIC_SITE2_BIS}},$H,1;
										$HASH2{$CRYPTIC_SITE2_BIS}[0]=$HASH2{$CRYPTIC_SITE2_BIS}[0]-1;
										if ($HASH2{$CRYPTIC_SITE2_BIS}[0]==0){
											delete $HASH2{$CRYPTIC_SITE2_BIS};
										}
									}
								}
							}
						}
					}
					foreach $CRYPTIC_SITE1_BIS (sort {$a<=>$b} keys %HASH1){
						$ID_REC1_BIS=$HASH1{$CRYPTIC_SITE1_BIS}[1];
						if ($ID_REC2 eq $ID_REC1_BIS){
							delete $HASH1{$CRYPTIC_SITE1_BIS};
						}
					}
					@CANDIDATE2=split(/\t/,qx{awk '\$10=="$ID_REC2"' $PATH_REF | grep ','});
					for $I(0..$#CANDIDATE2){
						if ($I<$#CANDIDATE2){
							print FILE_PRINT "$CANDIDATE2[$I]\t";
						}
						else{
							chomp $CANDIDATE2[$I];
							print FILE_PRINT "$CANDIDATE2[$I]\t";
						}
					}
					if ($#STOREHOUSE==0){
						if ($STOREHOUSE[$#STOREHOUSE]!~/-/ && $STOREHOUSE[$#STOREHOUSE]!~/#/){	
							print FILE_PRINT "css-rep:$STOREHOUSE[$#STOREHOUSE]\n";
						}
						elsif ($STOREHOUSE[$#STOREHOUSE]=~/-/){
							print FILE_PRINT "ESkip:$STOREHOUSE[$#STOREHOUSE]\n";
						}
						else{
							$SUBSTR=substr($STOREHOUSE[$#STOREHOUSE],1);
							print FILE_PRINT "IRet:$SUBSTR\n";
						}
					}
					else{
						for $J(0..$#STOREHOUSE){
							if ($J<$#STOREHOUSE){
								if ($J==0){
									if ($STOREHOUSE[$J]!~/-/ && $STOREHOUSE[$J]!~/#/){
										print FILE_PRINT "css-rep:$STOREHOUSE[$J];";
									}
									elsif ($STOREHOUSE[$J]=~/-/){
										print FILE_PRINT "ESkip:$STOREHOUSE[$J];";
									}
									else{
										$SUBSTR=substr($STOREHOUSE[$J],1);
										print FILE_PRINT "IRet:$STOREHOUSE[$J];";
									}
								}
								else{
									print FILE_PRINT "$STOREHOUSE[$J];";
								}
							}
							else{
								print FILE_PRINT "$STOREHOUSE[$J]\n";
							}
						}
					}
					
					@STOREHOUSE=();
				}
			# }
		}
		foreach $CRYPTIC_SITE1 (sort {$a<=>$b} keys %HASH1){
			$LENGTH_REC1=$HASH1{$CRYPTIC_SITE1}[0];
			# $LENGTH_REC1_ALT=@{$HASH1{$CRYPTIC_SITE1}};	$LENGTH_REC1_ALT=$LENGTH_REC1+1;
			# if ($LENGTH_REC1>=3){
				for $K(1..$LENGTH_REC1){
					$ID_REC1=$HASH1{$CRYPTIC_SITE1}[$K];
					@CANDIDATE1=split(/\t/,qx{awk '\$10=="$ID_REC1"' $PATH_REF | grep ','});
					for $L(0..$#CANDIDATE1){
						if ($L<$#CANDIDATE1){
							print FILE_PRINT "$CANDIDATE1[$L]\t";
						}
						else{
							chomp $CANDIDATE1[$L];
							print FILE_PRINT "$CANDIDATE1[$L]\t";
						}
					}
					print FILE_PRINT "NA\n";
				}
			# }
		}
		close FILE_PRINT;
		
		$PATH_PRINT="/path/to/" . "$REQUEST" . "/" . "acceptor.splice.sites/nucleotide.substitutions" . "/" . "$NAME" . "." . "$POS" . "." . "$OBS" . "/" . "reads.cryptic.splicing.v5.nodup.nomin_optB.txt";
		# qx{cat $PATH_PRINT_PRE | sort | uniq > $PATH_PRINT};	# It does not work entirely well without sorting
		qx{awk '!x[\$0]++' $PATH_PRINT_PRE > $PATH_PRINT};	# OR qx{cat -n $PATH_PRINT_PRE| sort -uk2 | sort -nk1 | cut -f2- > $PATH_PRINT};
		qx{rm -r $PATH_PRINT_PRE};
		
		%HASH1=();
		%HASH2=();
	}
	close FILE3;
	
	%MUT=();
	
}
close FILE2;