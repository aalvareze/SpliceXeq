# Distinction in cryptic splicing events

open FILE1, "<$ARGV[0]";        # $ARGV[0] = ../patients.nsubs.rna_seq.list
while (<FILE1>){
	chomp;
	@LINE=split(/\t/,$_);
	$REQUEST=$LINE[0];
	
	$REFERENCE="/path/to/" . "$REQUEST" . "/" . "acceptor.splice.sites/nucleotide.substitutions/candidates.mes";
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
		
		if (exists $MUT{$NAME}{$POS}{$OBS}){    # Several 23-mer sequences for the same mutation could be often an indication of different cryptic splicing events
			goto LINE;
		}
		$MUT{$NAME}{$POS}{$OBS}=1;
		
		$PATH_PRINT="/path/to/" . "$REQUEST" . "/" . "acceptor.splice.sites/nucleotide.substitutions" . "/" . "$NAME" . "." . "$POS" . "." . "$OBS" . "/" . "reads.cryptic.splicing.v4.nodup.min3_optA.txt";
		open FILE_PRINT, ">$PATH_PRINT";
		
		$PATH_REF="/path/to/" . "$REQUEST" . "/" . "acceptor.splice.sites/nucleotide.substitutions" . "/" . "$NAME" . "." . "$POS" . "." . "$OBS" . "/" . "reads.cryptic.splicing.v4.nodup.txt";		
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
			
			for $A(1..$#inBASES){
				$B=$A-1;
				$START=($inBASES[$B]+$inPOS[$B])+1;
				$END=$inPOS[$A];
				$SPAN=($END-$START)+1;	# It is necessary to add a correction (+1) to carefully determine the span of the gap
				if (exists $HASH{$START}{$END}{$SPAN}){
					$HASH{$START}{$END}{$SPAN}[0]=$HASH{$START}{$END}{$SPAN}[0]+1;
					$LENGTH=@{$HASH{$START}{$END}{$SPAN}};
					$HASH{$START}{$END}{$SPAN}[$LENGTH+1]=$ID;
				}
				else{
					$HASH{$START}{$END}{$SPAN}[0]=1;
					$HASH{$START}{$END}{$SPAN}[1]=$ID;
				}
			}
			
			END:
		}
		close FILE_REF;
		
		foreach $KEY_START (sort {$a<=>$b} keys %HASH){
			foreach $KEY_END (sort {$a<=>$b} keys %{$HASH{$KEY_START}}){
				foreach $KEY_SPAN (sort {$a<=>$b} keys %{$HASH{$KEY_START}{$KEY_END}}){
					if ($HASH{$KEY_START}{$KEY_END}{$KEY_SPAN}[0]>=3){	
						$LENGTH_REC=@{$HASH{$KEY_START}{$KEY_END}{$KEY_SPAN}};
						for $B(1..$LENGTH_REC){
							$ID_REC=$HASH{$KEY_START}{$KEY_END}{$KEY_SPAN}[$B];
							$CANDIDATE=qx{awk '\$10=="$ID_REC"' $PATH_REF | grep ','};
							# chomp $CANDIDATE;
				            # 
							# print FILE_PRINT "$CANDIDATE\n";
							print FILE_PRINT "$CANDIDATE";
						}
					}
				}
			}
		}
		
		close FILE_PRINT;
		
		%HASH=();
	}
	close FILE2;
	
	%MUT=();
}
close FILE1;