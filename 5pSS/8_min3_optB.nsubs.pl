# Distinction in cryptic splicing sites

open FILE1, "<$ARGV[0]";        # $ARGV[0] = ../patients.nsubs.rna_seq.list
while (<FILE1>){
	chomp;
	@LINE=split(/\t/,$_);
	$REQUEST=$LINE[0];
	
	$REFERENCE="/path/to/" . "$REQUEST" . "/" . "donor.splice.sites/nucleotide.substitutions/candidates.mes";
	open FILE2, "<$REFERENCE";
	LINE:
	while (<FILE2>){
		chomp;
		
		if ($_=~/CHR/){
			goto LINE;
		}
		
		@HEADER1=split(/\t/,$_);
		$POS=$HEADER1[1];
		$OBS=$HEADER1[3];
		$NAME=$HEADER1[4];
		$STRAND=$HEADER1[5];
		
		if (exists $MUT{$NAME}{$POS}{$OBS}){    # Several 9-mer sequences for the same mutation could be often an indication of different cryptic splicing events
			goto LINE;
		}
		$MUT{$NAME}{$POS}{$OBS}=1;
		
		$PATH_PRINT_PRE="/path/to/" . "$REQUEST" . "/" . "donor.splice.sites/nucleotide.substitutions" . "/" . "$NAME" . "." . "$POS" . "." . "$OBS" . "/" . "reads.cryptic.splicing.v5.nodup.min3_optB_PRE.txt";
		open FILE_PRINT, ">$PATH_PRINT_PRE";
		
		$PATH_REF="/path/to/" . "$REQUEST" . "/" . "donor.splice.sites/nucleotide.substitutions" . "/" . "$NAME" . "." . "$POS" . "." . "$OBS" . "/" . "reads.cryptic.splicing.v5.nodup.txt";		
		$COUNT=0;
		open FILE_REF, "<$PATH_REF";
		while (<FILE_REF>){
			chomp;

			if ($COUNT<5){
				print FILE_PRINT "$_\n";
				$COUNT++;
				goto END;
			}
			
			@HEADER2=split(/\t/,$_);
			$ID=$HEADER2[9];
			$inBASES_pre=substr($HEADER2[18],0,-1);
			@inBASES=split(/,/,$inBASES_pre);
			$inPOS_pre=substr($HEADER2[20],0,-1);
			@inPOS=split(/,/,$inPOS_pre);
			
			for $A(0..$#inBASES){
				if (($inPOS[$A]+1)<=$POS && $POS<=($inBASES[$A]+$inPOS[$A])){
					if ($STRAND==1){
						$lastPOS=$inBASES[$A]+$inPOS[$A];
						if (exists $HASH{$lastPOS}){
							$HASH{$lastPOS}[0]=$HASH{$lastPOS}[0]+1;
							$LENGTH=@{$HASH{$lastPOS}};
							$HASH{$lastPOS}[$LENGTH+1]=$ID;
						}
						else{
							$HASH{$lastPOS}[0]=1;
							$HASH{$lastPOS}[1]=$ID;
						}
					}
					else{
						$firstPOS=$inPOS[$A]+1;
						if (exists $HASH{$firstPOS}){
							$HASH{$firstPOS}[0]=$HASH{$firstPOS}[0]+1;
							$LENGTH=@{$HASH{$firstPOS}};
							$HASH{$firstPOS}[$LENGTH+1]=$ID;
						}
						else{
							$HASH{$firstPOS}[0]=1;
							$HASH{$firstPOS}[1]=$ID;
						}
					}
				}
				else{	# The mutation could be located in a gap
					if ($A>0){
						$B=$A-1;
						if ((($inBASES[$B]+$inPOS[$B])+1)<=$POS && $POS<=$inPOS[$A]){
							if ($STRAND==1){
								$lastPOS=$inBASES[$B]+$inPOS[$B];
								if (exists $HASH{$lastPOS}){
									$HASH{$lastPOS}[0]=$HASH{$lastPOS}[0]+1;
									$LENGTH=@{$HASH{$lastPOS}};
									$HASH{$lastPOS}[$LENGTH+1]=$ID;
								}
								else{
									$HASH{$lastPOS}[0]=1;
									$HASH{$lastPOS}[1]=$ID;
								}
							}
							else{
								$firstPOS=$inPOS[$A]+1;
								if (exists $HASH{$firstPOS}){
									$HASH{$firstPOS}[0]=$HASH{$firstPOS}[0]+1;
									$LENGTH=@{$HASH{$firstPOS}};
									$HASH{$firstPOS}[$LENGTH+1]=$ID;
								}
								else{
									$HASH{$firstPOS}[0]=1;
									$HASH{$firstPOS}[1]=$ID;
								}
							}
						}
					}
				}
			}
			
			END:
		}
		close FILE_REF;
		
		foreach $CRYPTIC_SITE (sort {$a<=>$b} keys %HASH){
			if ($HASH{$CRYPTIC_SITE}[0]>=3){	
				$LENGTH_REC=@{$HASH{$CRYPTIC_SITE}};
				for $C(1..$LENGTH_REC){
					$ID_REC=$HASH{$CRYPTIC_SITE}[$C];
					$CANDIDATE=qx{awk '\$10=="$ID_REC"' $PATH_REF | grep ','};
					# chomp $CANDIDATE;
				    # 
					# print FILE_PRINT "$CANDIDATE\n";
					print FILE_PRINT "$CANDIDATE";
				}
			}
		}
		
		close FILE_PRINT;
		$PATH_PRINT="/path/to/" . "$REQUEST" . "/" . "donor.splice.sites/nucleotide.substitutions" . "/" . "$NAME" . "." . "$POS" . "." . "$OBS" . "/" . "reads.cryptic.splicing.v5.nodup.min3_optB.txt";
		# qx{cat $PATH_PRINT_PRE | sort | uniq > $PATH_PRINT};
		qx{cat $PATH_PRINT_PRE | uniq > $PATH_PRINT};
		qx{rm -r $PATH_PRINT_PRE};
		
		%HASH=();
	}
	close FILE2;
	
	%MUT=();
}
close FILE1;