open FILE_cases, "<$ARGV[0]";	# $ARGV[0] = ../cases.nsubs.rna_seq.list
while (<FILE_cases>){
	chomp;
	@PARAMETERS=split(/\t/,$_);
	$REQUEST=$PARAMETERS[0];
	
	$HASH{$REQUEST}=1;
}
close FILE_cases;

foreach $KEY_REQUEST (sort {$a<=>$b} keys %HASH){

	$PATH="/path/to/" . "$KEY_REQUEST" . "/" . "MES" . "/" . "acceptorss" . "/" . "nsubs";
	
	$MESfiles="/path/to/scripts" . "/" . "maxentscan.program.files.list";
	
	qx{cp $MESfiles/* $PATH};
	
	$output_minus_strand="$PATH" . "/" . "minusSTRAND.onlyFWD.ref";
	$output_bioperl="$PATH" . "/" . "FWDandRV-COMP.ref";
	$print_rv_comp="$PATH" . "/" . "onlyRV-COMP.ref";
	$print_rv_comp_mut="$PATH" . "/" . "onlyRV-COMP.obs";	
	
	$output_plus_strand="$PATH" . "/" . "plusSTRAND.onlyFWD.ref";
	$print_fwd="$PATH" . "/" . "onlyFWD.ref";
	$print_fwd_mut="$PATH" . "/" . "onlyFWD.obs";
	
	$output_MES_ref="$PATH" . "/" . "maxentscan.ref";
	$output_MES_mut="$PATH" . "/" . "maxentscan.obs";
	
	$COUNTmuts=0;
	$INPUT="/path/to/" . "$KEY_REQUEST" . "/" . "nsubs.maf";
	open FILE, "<$INPUT";
	LINE:
	while (<FILE>){
		chomp;

		if ($_=~/REQUEST NUMBER/){
			goto LINE;
		}
		
		@LIST=split(/\t/,$_);
		$CHR=$LIST[1];
		$POS=$LIST[2];
		$REF=$LIST[3];
		$OBS=$LIST[4];
		$NAME=$LIST[5];
		$STRAND=$LIST[6];
	
		if ($STRAND==-1){
			$END=$POS+22;
			$START=$POS-22;
	
			qx{/path/to/samtools-0.1.19/samtools faidx /path/to/hg19.fa $CHR:$START-$END > $output_minus_strand};
			qx{perl /path/to/k-mers.options.FWDandRV-COMP.pl $output_minus_strand $output_bioperl};
	
			$A=0;
			open FILE1, "<$output_bioperl";
			while (<FILE1>){
				chomp;
				$A++;
				
				open FILE_rv_comp_ref, ">$print_rv_comp";
				open FILE_rv_comp_mut, ">$print_rv_comp_mut";
				
				if ($A==4){
					@CANONseq=split(//,$_);
					
					$COUNT=0;
					for $B(0..($#CANONseq-22)){
						for $C($B..(22+$B)){
							$FRAGMENT[$C-$COUNT]=$CANONseq[$C];
						}
						$MERGER=join("",@FRAGMENT);
						$COUNT++;
						
						print FILE_rv_comp_ref ">RV-COMP (canon seq)\n$MERGER\n";

						if ($OBS eq T){
							$FRAGMENT[22-$B]="A";
						}
						if ($OBS eq C){
							$FRAGMENT[22-$B]="G";
						}
						if ($OBS eq G){
							$FRAGMENT[22-$B]="C";
						}
						if ($OBS eq A){
							$FRAGMENT[22-$B]="T";
						}

						print FILE_rv_comp_mut ">RV-COMP (mut seq)\n";
						print FILE_rv_comp_mut join ("",@FRAGMENT), "\n";
					}

					qx{perl /path/to/score3.pl $print_rv_comp >> $output_MES_ref};
					qx{perl /path/to/score3.pl $print_rv_comp_mut >> $output_MES_mut};

					for $D(0..($#CANONseq-22)){
						$END_MERGER=$END-$D;
						$START_MERGER=$END_MERGER-22;
					
						$mutPOSinSEQ=(22-$D)+1;	# Values of the variable $mutPOSinSEQ for genes located in the minus strand even working on its RV-COMP counterpart
					
						$MUTS{$COUNTmuts}{$CHR}{$POS}{$OBS}[0]=$START_MERGER;
						$MUTS{$COUNTmuts}{$CHR}{$POS}{$OBS}[1]=$END_MERGER;
						$MUTS{$COUNTmuts}{$CHR}{$POS}{$OBS}[2]=$REF;
						$MUTS{$COUNTmuts}{$CHR}{$POS}{$OBS}[3]=$NAME;
						$MUTS{$COUNTmuts}{$CHR}{$POS}{$OBS}[4]=$STRAND;
						$MUTS{$COUNTmuts}{$CHR}{$POS}{$OBS}[5]=$mutPOSinSEQ;
						$COUNTmuts++;
					}	
				}

				close FILE_rv_comp_ref;
				close FILE_rv_comp_mut;
				
			}
			close FILE1;
		}
		elsif ($STRAND==1){
			$END=$POS+22;
			$START=$POS-22;
			
			qx{/path/to/samtools-0.1.19/samtools faidx /path/to/hg19.fa $CHR:$START-$END > $output_plus_strand};
			
			$A=0;
			open FILE2, "<$output_plus_strand";
			while (<FILE2>){
				chomp;
				$A++;				

				open FILE_fwd_ref, ">$print_fwd";
				open FILE_fwd_mut, ">$print_fwd_mut";

				if ($A==2){
					@CANONseq=split(//,$_);
					
					$COUNT=0;
					for $B(0..($#CANONseq-22)){
						for $C($B..(22+$B)){
							$FRAGMENT[$C-$COUNT]=$CANONseq[$C];
						}
						$MERGER=join("",@FRAGMENT);
						$COUNT++;
						
						print FILE_fwd_ref ">FWD (canon seq)\n$MERGER\n";

						$FRAGMENT[22-$B]=$OBS;

						print FILE_fwd_mut ">FWD (mut seq)\n";
						print FILE_fwd_mut join ("",@FRAGMENT), "\n";
					}

					qx{perl /path/to/score3.pl $print_fwd >> $output_MES_ref};
					qx{perl /path/to/score3.pl $print_fwd_mut >> $output_MES_mut};

					for $D(0..($#CANONseq-22)){
						$START_MERGER=$START+$D;
						$END_MERGER=$START_MERGER+22;
					
						$mutPOSinSEQ=(22-$D)+1;	# Values of the variable $mutPOSinSEQ for genes located in the plus strand
					
						$MUTS{$COUNTmuts}{$CHR}{$POS}{$OBS}[0]=$START_MERGER;
						$MUTS{$COUNTmuts}{$CHR}{$POS}{$OBS}[1]=$END_MERGER;
						$MUTS{$COUNTmuts}{$CHR}{$POS}{$OBS}[2]=$REF;
						$MUTS{$COUNTmuts}{$CHR}{$POS}{$OBS}[3]=$NAME;
						$MUTS{$COUNTmuts}{$CHR}{$POS}{$OBS}[4]=$STRAND;
						$MUTS{$COUNTmuts}{$CHR}{$POS}{$OBS}[5]=$mutPOSinSEQ;
						$COUNTmuts++;				
					}
				}
				
				close FILE_fwd_ref;
				close FILE_fwd_mut;
				
			}
			close FILE2;	
		}
		else{
			goto LINE;
		}
	}
	close FILE;

	qx{rm -r $output_minus_strand};
	qx{rm -r $output_bioperl};
	qx{rm -r $print_rv_comp};
	qx{rm -r $print_rv_comp_mut};
	
	qx{rm -r $output_plus_strand};
	qx{rm -r $print_fwd};
	qx{rm -r $print_fwd_mut};

	$PRINT="$PATH" . "/" . "muts.mes";
	open FILE_print, ">$PRINT";
	print FILE_print "CHR\tPOS\tREF\tALT\tSYMBOL\tSTRAND\tSTART\tEND\trefSEQ\trefSCORE\taltSEQ\taltPOSinSEQ\taltSCORE\tdeltaSCORE\n";

	$COUNT_MES_ref=0;
	open FILE_MES_ref, "<$output_MES_ref";
	while (<FILE_MES_ref>){
		chomp;
		@OUTPUT1=split(/\t/,$_);
		$SEQ_ref=$OUTPUT1[0];
		$SCORE_ref=$OUTPUT1[1];
		
		$MESref{$COUNT_MES_ref}{$SEQ_ref}[0]=$SCORE_ref;
		$COUNT_MES_ref++;
	}
	close FILE_MES_ref;

	$COUNT_MES_mut=0;
	open FILE_MES_mut, "<$output_MES_mut";
	while (<FILE_MES_mut>){
		chomp;
		@OUTPUT2=split(/\t/,$_);
		$SEQ_mut=$OUTPUT2[0];
		$SCORE_mut=$OUTPUT2[1];
		
		$MESmut{$COUNT_MES_mut}{$SEQ_mut}[0]=$SCORE_mut;
		$COUNT_MES_mut++;
	}
	close FILE_MES_mut;

	foreach $KEY_COUNT (sort {$a<=>$b} keys %MUTS){
		foreach $KEY_CHR (sort {$a<=>$b} keys %{$MUTS{$KEY_COUNT}}){
			foreach $KEY_POS (sort {$a<=>$b} keys %{$MUTS{$KEY_COUNT}{$KEY_CHR}}){
				foreach $KEY_OBS (sort {$a<=>$b} keys %{$MUTS{$KEY_COUNT}{$KEY_CHR}{$KEY_POS}}){
					$START_MERGER_prima=$MUTS{$KEY_COUNT}{$KEY_CHR}{$KEY_POS}{$KEY_OBS}[0];
					$END_MERGER_prima=$MUTS{$KEY_COUNT}{$KEY_CHR}{$KEY_POS}{$KEY_OBS}[1];
					$REF_prima=$MUTS{$KEY_COUNT}{$KEY_CHR}{$KEY_POS}{$KEY_OBS}[2];
					$NAME_prima=$MUTS{$KEY_COUNT}{$KEY_CHR}{$KEY_POS}{$KEY_OBS}[3];
					$STRAND_prima=$MUTS{$KEY_COUNT}{$KEY_CHR}{$KEY_POS}{$KEY_OBS}[4];
					$mutPOSinSEQ_prima=$MUTS{$KEY_COUNT}{$KEY_CHR}{$KEY_POS}{$KEY_OBS}[5];
					
					foreach $KEY_SEQ_ref (sort {$a<=>$b} keys %{$MESref{$KEY_COUNT}}){
						$SCORE_ref_prima=$MESref{$KEY_COUNT}{$KEY_SEQ_ref}[0];
						
						foreach $KEY_SEQ_mut (sort {$a<=>$b} keys %{$MESmut{$KEY_COUNT}}){
							$SCORE_mut_prima=$MESmut{$KEY_COUNT}{$KEY_SEQ_mut}[0];
							
							$DIFFERENCE=$SCORE_mut_prima-$SCORE_ref_prima;
							
							print FILE_print "$KEY_CHR\t$KEY_POS\t$REF_prima\t$KEY_OBS\t$NAME_prima\t$STRAND_prima\t$START_MERGER_prima\t$END_MERGER_prima\t$KEY_SEQ_ref\t$SCORE_ref_prima\t$KEY_SEQ_mut\t$mutPOSinSEQ_prima\t$SCORE_mut_prima\t$DIFFERENCE\n";
						}
					}
				}
			}
		}
	}
	
	%MUTS=();
	%MESref=();
	%MESmut=();

	close FILE_print;
	
	qx{rm -r $output_MES_ref};
	qx{rm -r $output_MES_mut};
	
}