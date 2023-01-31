# PCR duplicate removal

open FILE1, "<$ARGV[0]";        # $ARGV[0] = ../patients.nsubs.rna_seq.list
while (<FILE1>){
	chomp;
	@LINE=split(/\t/,$_);
	$ID=$LINE[0];
	
	$REFERENCE="/path/to/" . "$ID" . "/" . "donor.splice.sites/nucleotide.substitutions/candidates.mes";
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
		
		if (exists $MUT{$NAME}{$POS}{$OBS}){    # Several 9-mer sequences for the same mutation could be often an indication of different cryptic splicing events
			goto LINE;
		}
		$MUT{$NAME}{$POS}{$OBS}=1;
		
		$PATH_PRINT1="/path/to/" . "$ID" . "/" . "donor.splice.sites/nucleotide.substitutions" . "/" . "$NAME" . "." . "$POS" . "." . "$OBS" . "/" . "reads.cryptic.splicing.v5.nodup.txt";
		open FILE_PRINT1, ">$PATH_PRINT1";
		
		$PATH_REF1="/path/to/" . "$ID" . "/" . "donor.splice.sites/nucleotide.substitutions" . "/" . "$NAME" . "." . "$POS" . "." . "$OBS" . "/" . "reads.cryptic.splicing.v5.txt";		
		$COUNT1=0;
		open FILE_REF1, "<$PATH_REF1";
		while (<FILE_REF1>){
			chomp;

			if ($COUNT1<5){
				print FILE_PRINT1 "$_\n";
				$COUNT1++;
				goto END;
			}
			
			@HEADER2=split(/\t/,$_);
			$inBASES1=$HEADER2[18];
			$inPOS1=$HEADER2[20];
			
			if (exists $HASH1{$inBASES1}{$inPOS1}){
				goto END;
			}
			$HASH1{$inBASES1}{$inPOS1}=1;
			
			print FILE_PRINT1 "$_\n";
			
			END:
		}
		close FILE_REF1;
		
		close FILE_PRINT1;
		
		%HASH1=();
		
		$PATH_PRINT2="/path/to/" . "$ID" . "/" . "donor.splice.sites/nucleotide.substitutions" . "/" . "$NAME" . "." . "$POS" . "." . "$OBS" . "/" . "reads.alternative.splicing.v5.nodup.txt";
		open FILE_PRINT2, ">$PATH_PRINT2";
		
		$PATH_REF2="/path/to/" . "$ID" . "/" . "donor.splice.sites/nucleotide.substitutions" . "/" . "$NAME" . "." . "$POS" . "." . "$OBS" . "/" . "reads.alternative.splicing.v5.txt";		
		$COUNT2=0;
		open FILE_REF2, "<$PATH_REF2";
		while (<FILE_REF2>){
			chomp;

			if ($COUNT2<5){
				print FILE_PRINT2 "$_\n";
				$COUNT2++;
				goto END;
			}
			
			@HEADER3=split(/\t/,$_);
			$inBASES2=$HEADER3[18];
			$inPOS2=$HEADER3[20];
			
			if (exists $HASH2{$inBASES2}{$inPOS2}){
				goto END;
			}
			$HASH2{$inBASES2}{$inPOS2}=1;
			
			print FILE_PRINT2 "$_\n";
			
			END:
		}
		close FILE_REF2;
		
		close FILE_PRINT2;
		
		%HASH2=();
		
		$PATH_PRINT3="/path/to/" . "$ID" . "/" . "donor.splice.sites/nucleotide.substitutions" . "/" . "$NAME" . "." . "$POS" . "." . "$OBS" . "/" . "reads.unknown.v5.nodup.txt";
		open FILE_PRINT3, ">$PATH_PRINT3";
		
		$PATH_REF3="/path/to/" . "$ID" . "/" . "donor.splice.sites/nucleotide.substitutions" . "/" . "$NAME" . "." . "$POS" . "." . "$OBS" . "/" . "reads.unknown.v5.txt";		
		$COUNT3=0;
		open FILE_REF3, "<$PATH_REF3";
		while (<FILE_REF3>){
			chomp;

			if ($COUNT3<5){
				print FILE_PRINT3 "$_\n";
				$COUNT3++;
				goto END;
			}
			
			@HEADER4=split(/\t/,$_);
			$inBASES3=$HEADER4[18];
			$inPOS3=$HEADER4[20];
			
			if (exists $HASH3{$inBASES3}{$inPOS3}){
				goto END;
			}
			$HASH3{$inBASES3}{$inPOS3}=1;
			
			print FILE_PRINT3 "$_\n";
			
			END:
		}
		close FILE_REF3;
		
		close FILE_PRINT3;
		
		%HASH3=();
		
	}
	close FILE2;
	
	%MUT=();
}
close FILE1;
