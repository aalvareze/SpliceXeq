# This script filters duplicates and formats reads to make them compatible with BLAT. Moreover, it distinguishes reads named equally

open FILE1, "<$ARGV[0]";	# $ARGV[0] = ../patients.nsubs.rna_seq.list
while (<FILE1>){
	chomp;
	@LINE1=split(/\t/,$_);
	$REQUEST=$LINE1[0];
	
	$HASH1{$REQUEST}=1;
}
close FILE1;

foreach $KEY_REQUEST (sort {$a<=>$b} keys %HASH1){

	$PATH_FILE="/path/to/" . "$KEY_REQUEST" . "/" . "donor.splice.sites/nucleotide.substitutions/candidates.mes";
	open FILE2, "<$PATH_FILE";
	LINE:
	while (<FILE2>){
		chomp;
		
		if ($_=~/CHR/){
			goto LINE;
		}
		
		@LINE2=split(/\t/,$_);
		$POS=$LINE2[1];
		$OBS=$LINE2[3];
		$NAME=$LINE2[4];

		if (exists $HASH2{$NAME}{$POS}{$OBS}){	# Several 9-mer sequences for the same mutation could be often an indication of different cryptic splicing events
			goto LINE;
		}
		$HASH2{$NAME}{$POS}{$OBS}=1;
		$PATH_PRINT="/path/to/" . "$KEY_REQUEST" . "/" . "donor.splice.sites/nucleotide.substitutions" . "/" . "$NAME" . "." . "$POS" . "." . "$OBS";
		
		$PRINT="$PATH_PRINT" . "/" . "reads.fa";
		open FILE_print, ">$PRINT";
		
		$FILE="$PATH_PRINT" . "/" . "reads.mutpos.sam";
		open FILE3, "<$FILE";
		while (<FILE3>){
			chomp;
			@LINE3=split(/\t/,$_);
			$ID=$LINE3[0];
			$SEQ=$LINE3[9];
			
			if (exists $HASH3{$ID}){
				$HASH3{$ID}[0]=$HASH3{$ID}[0]+1;
				$ID_BIS="$ID" . "_" . "$HASH3{$ID}[0]";
				
				print FILE_print ">$ID_BIS\n$SEQ\n";
			}
			else{
				$HASH3{$ID}[0]=1;
				$ID_BIS="$ID" . "_" . "$HASH3{$ID}[0]";
				
				print FILE_print ">$ID_BIS\n$SEQ\n";
			}
		}
		close FILE3;
		
		close FILE_print;
		
		%HASH3=();

	}
	close FILE2;
	
	%HASH2=();
	
}