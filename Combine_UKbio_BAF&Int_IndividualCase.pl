#!/usr/local/bin/perl

# These files are read just once, so they do not need to be copied to local disk
$file1 = "/data/LNG/CORNELIS_TEMP/WILL_PRKN/parsed/prknbaf.txt";
$file2 = "/data/LNG/CORNELIS_TEMP/WILL_PRKN/parsed/prknl2ratio.txt";
$variants = "/data/LNG/CORNELIS_TEMP/WILL_PRKN/parsed/variantsbp.bim";

$linenum=0;
# read the prknbaf file into an array
print "Reading in file $file1\n";
open (INFILE1, "< $file1")
    || die "cannot open file $file1 $!*";
while (<INFILE1>) {
    # split a single line into an array
    my @tmp = split(/ /,$_);
    # assign that array to an element of the prknbaf hash
    $prknbaf{"$linenum"} = [@tmp];
    #increase the counter
    $linenum++;
}
close (INFILE1);
print "Finished reading in file $file1\n";

# read the prknl2ratio file into an array
print "Reading in file $file2\n";
$linenum = 0;
open (INFILE2, "< $file2") 
    || die "cannot open file $file2 $!*";
while (<INFILE2>) {
    # split a single line into an array
    @tmp = split(/ /,$_);
    # assign that array to an element of the prkn2ratio hash
    $prknl2ratio{"$linenum"} = [@tmp];
    #increase the counter
    $linenum++;
}
close (INFILE2);
print "Finished reading in file $file2\n";

# read the variants file into an array
print "Reading in variants file\n";
$linenum = 0;
open (INFILE3, "< $variants")
    || die "cannot open file $variants $!*";
while (<INFILE3>) {
    chop $_;
    $variants{"$linenum"} = $_;
    #increase the counter
    $linenum++;
}
close (INFILE3);
print "Finished reading in file $variants\n";

system("mkdir -p output");

foreach $pat (0..488376) {
    # divide into directories with 1000 patients each
    $olddir = $dir;
    $dir = int(($pat + 1) /1000);
    # create a new directory if necessary 
    if ($dir ne $olddir) {
	print "Creating directory ${dir}/488\n";
	system("mkdir -p output/$dir");
    }
    my $outfile = "output/${dir}/pt" . ($pat + 1);
    open (FILE, "> $outfile") 
	|| die "Cannot open file output/${dir}/pt${pat} $!*";
    foreach $line (0..1105) {
	print FILE "$variants{$line}\t$prknl2ratio{$line}[$pat]\t$prknbaf{$line}[$pat]\n";
    }
    close ($FILE);
}
	


    
