use strict;
use warnings;

#>TIPP_plastid_r1_edge_1=Tul-0_1 [41 - 757] (REVERSE SENSE)
#location record in getorf output is 1-based;

open IN1,"$ARGV[0]";
while(<IN1>){
	if( />/ ){
		$_=~s/>//;
		my$name=(split /#/,$_)[1];
		$name=(split / /,$name)[0];
		my$chr=(split / /,$_)[0];
		$chr=~ s/_[^_]*$//;
		my$geneID=(split / /,$_)[0];
		my$location=(split /\[/,$_)[1];
		$location=(split /\]/,$location)[0];
		my($start,$end)=(split / /,$location)[0,2];
		my$strand="+";
		if($_=~m/(REVERSE SENSE)/){
			$strand="-";
			print "$chr\tgetorf\tgene\t$end\t$start\t.\t$strand\t.\tID=$geneID\n";
			print "$chr\tgetorf\tCDS\t$end\t$start\t.\t$strand\t.\tID=$geneID;Parent=$geneID\n";
		}
		else{
			print "$chr\tgetorf\tgene\t$start\t$end\t.\t$strand\t.\tID=$geneID\n";
			print "$chr\tgetorf\tCDS\t$start\t$end\t.\t$strand\t.\tID=$geneID;Parent=$geneID\n";
		}
	}
}
