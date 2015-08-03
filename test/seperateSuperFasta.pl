use strict;

my $file=$ARGV[0];

my $init=0;
open(IN,$file);
while(my $line=<IN>) {
	chomp($line);
	
	if($line =~ /^>/) {
		my $fileName=$line;
		
		$fileName =~ s/ //g; # remove any weird white space
		$fileName =~ s/\t/-/g; # remove any weird white space
		$fileName =~ s/^>//;
		$fileName = (split(/\|/,$fileName))[0];
		$fileName = "chr".$fileName if($fileName !~ /^chr/);
		
		close(OUT) if($init != 0);
		print "opening $fileName.fa ...\n";
		open(OUT,">".$fileName.".fa");
		
		print OUT ">".$fileName;
		
		$init=1;
	} else {
		$line=uc($line);
		print OUT "\n$line";
	}
}
