# This script tries to identify reliable splicing mutations

open FILE, "<$ARGV[0]";	# $ARGV[0] = ../cases.nsubs.rna_seq.list
while (<FILE>){
	chomp;
	@LINE=split(/\t/,$_);
	$REQUEST=$LINE[0];
	
	$HASH{$REQUEST}=1;
}
close FILE;

foreach $KEY_REQUEST (sort {$a<=>$b} keys %HASH){
	
	$PRINT1="/path/to/" . "$KEY_REQUEST" . "/" . "MES" . "/" . "donorss" . "/" . "nsubs" . "/" . "candidates.mes";
	open FILE_PRINT1, ">$PRINT1";
	
	$PATH1="/path/to/" . "$KEY_REQUEST" . "/" . "MES" . "/" . "donorss" . "/" . "nsubs" . "/" . "muts.mes";
	open FILE_PATH1, "<$PATH1";
	while (<FILE_PATH1>){
		chomp;
		@LINE1=split(/\t/,$_);
		$obsSCORE1=$LINE1[12];
		$DIFFERENCE1=$LINE1[13];
		
		if ($_=~/CHR/){
			print FILE_PRINT1 "$_\n";
		}
		elsif ($obsSCORE1>0 && $DIFFERENCE1>=5){
			print FILE_PRINT1 "$_\n";
		}
		else{
			next;
		}
	}
	close FILE_PATH1;
	
	close FILE_PRINT1;
	
	$PRINT2="/path/to/" . "$KEY_REQUEST" . "/" . "MES" . "/" . "acceptorss" . "/" . "nsubs" . "/" . "candidates.mes";
	open FILE_PRINT2, ">$PRINT2";
	
	$PATH2="/path/to/" . "$KEY_REQUEST" . "/" . "MES" . "/" . "acceptorss" . "/" . "nsubs" . "/" . "muts.mes";
	open FILE_PATH2, "<$PATH2";
	while (<FILE_PATH2>){
		chomp;
		@LINE2=split(/\t/,$_);
		$obsSCORE2=$LINE2[12];
		$DIFFERENCE2=$LINE2[13];
		
		if ($_=~/CHR/){
			print FILE_PRINT2 "$_\n";
		}
		elsif ($obsSCORE2>0 && $DIFFERENCE2>=5){
			print FILE_PRINT2 "$_\n";
		}
		else{
			next;
		}
	}
	close FILE_PATH2;
	
	close FILE_PRINT2;	
	
}