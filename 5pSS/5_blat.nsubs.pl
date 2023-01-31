# BLAT

open FILE1, "<$ARGV[0]";	# $ARGV[0] = ../patients.nsubs.rna_seq.list
while (<FILE1>){
	chomp;
	@LINE1=split(/\t/,$_);
	$ID=$LINE1[0];
	
	$HASH1{$ID}=1;
}
close FILE1;

foreach $KEY_ID (sort {$a<=>$b} keys %HASH1){

	$PATH_FILE="/path/to/" . "$KEY_ID" . "/" . "donor.splice.sites/nucleotide.substitutions/candidates.mes";
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
		$PATH_PRINT="/path/to/" . "$KEY_ID" . "/" . "donor.splice.sites/nucleotide.substitutions" . "/" . "$NAME" . "." . "$POS" . "." . "$OBS";
		
		$FILE="$PATH_PRINT" . "/" . "reads.fa";
		$PRINT="$PATH_PRINT" . "/" . "reads.blat.txt";
		
		if (! -s $FILE){	# "reads.fa" may be an empty file
			open FILE_print, ">$PRINT";
			close FILE_print;
		}
		else{
			qx{gfClient -out=psl localhost 9050 "" $FILE stdout > $PRINT};
		}
	}
	close FILE2;
	
	%HASH2=();

}