# This script parses the BAM file 'ID.Aligned.sortedByCoord.out.bam' written by STAR

open FILE1, "<$ARGV[0]";	# $ARGV[0] = ../patients.nsubs.rna_seq.list
while (<FILE1>){
	chomp;
	@LINE1=split(/\t/,$_);
	$ID=$LINE1[0];
	$BAM=$LINE1[1];
	
	$HASH1{$ID}[0]=$BAM;
}
close FILE1;

$samtools = "/path/to/samtools-0.1.19/samtools";

foreach $KEY_ID (sort {$a<=>$b} keys %HASH1){

	$PATH_BAM=$HASH1{$KEY_ID}[0];
	
	$PATH_FILE="/path/to/" . "$KEY_ID" . "/" . "donor.splice.sites/nucleotide.substitutions/candidates.mes";
	open FILE2, "<$PATH_FILE";
	LINE:
	while (<FILE2>){
		chomp;
		
		if ($_=~/CHR/){
			goto LINE;
		}
		
		@LINE2=split(/\t/,$_);
		$CHR=$LINE2[0];
		$POS=$LINE2[1];
		$REF=$LINE2[2];
		$OBS=$LINE2[3];
		$NAME=$LINE2[4];

		if (exists $HASH2{$NAME}{$POS}{$OBS}){	# Several 9-mer sequences for the same mutation could be often an indication of different cryptic splicing events
			goto LINE;
		}
		
		$HASH2{$NAME}{$POS}{$OBS}=1;
		$PATH_PRINT="/path/to/" . "$KEY_ID" . "/" . "donor.splice.sites/nucleotide.substitutions" . "/" . "$NAME" . "." . "$POS" . "." . "$OBS";
		qx{mkdir $PATH_PRINT};

		qx{$samtools view $PATH_BAM $CHR:$POS-$POS >> $PATH_PRINT/reads.mutpos.sam};
	}
	close FILE2;
	
	%HASH2=();

}