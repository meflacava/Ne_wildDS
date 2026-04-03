#!/usr/bin/perl

if ($#ARGV == 3) { #read in arguments
    $file1 = $ARGV[0];
    $file2 = $ARGV[1];
    $barcode = $ARGV[2];
    $prefix = $ARGV[3];
} else {
    die;
}

@commas = split(/\,/, $barcode); #if barcode column has comma in it, takes everything up to comma, but if no comma in barcode column, then just takes whole barcode, so treat @commas like @barcode still
$barcode_length = length($commas[0]); #length of barcode

$x=0;
while ($x <= $#commas) { #make empty fastq files
    $hash_r1{$commas[$x]} = $prefix . "_RA_" . $commas[$x] . ".fastq";
    $hash_r2{$commas[$x]} = $prefix . "_RB_" . $commas[$x] . ".fastq";
    $filename_r1 = $hash_r1{$commas[$x]};
    $filename_r2 = $hash_r2{$commas[$x]};
    open($filename_r1, ">$filename_r1") or die;
    open($filename_r2, ">$filename_r2") or die;
    $x++;
}


open(FILE1, "<$file1") or die;
open(FILE2, "<$file2") or die;


while (<FILE1>) { #file 1 = multiplexed R1 reads, file 2 = multiplexed R2 reads

    $f1a = $_; #the four lines of R1 reads
    $f1b = <FILE1>;
    $f1c = <FILE1>;
    $f1d = <FILE1>;

    $f2a = <FILE2>; #the four lines of R2 reads
        $f2b = <FILE2>;
        $f2c = <FILE2>;
        $f2d = <FILE2>;
    
    $bc1 = substr($f1b,0,$barcode_length); #list of barcodes for all multiplexed reads in R1
    $bc2 = substr($f2b,0,$barcode_length); #R2 read barcodes? but doesn't R2 not have barcodes?
    
    if ($hash_r1{$bc1} ne "" && $hash_r1{$bc2} eq "")  { #if target read in R2 matches barcode of individual with file currently open (eq = equal, ne = not equal)

        $f1b_2 = substr($f1b, $barcode_length, length($f1b)); #grab read, but exclude barcode
        $f1d_2 = substr($f1d, $barcode_length, length($f1d)); #grab quality scores, except for barcode section

        $out1 = $hash_r1{$bc1};
        $out2 = $hash_r2{$bc1};

        print $out1 $f1a . $f1b_2 . $f1c . $f1d_2; #grab lines 1,3 of read as is, but grab shortened lines 2,4 with barcode removed from read, and put into individual's fastq file
        print $out2 $f2a . $f2b . $f2c . $f2d; #grab same read and put it in individual's second fastq file

    } elsif ($hash_r1{$bc1} eq "" && $hash_r1{$bc2} ne "")  { #if target read in R1 matches barcode of individual with file currently open (MAYBE only relevant to dual indexed reads? Maybe never does anything for my reads since I have a single barcode?)

        $f2b_2 = substr($f2b, $barcode_length, length($f2b));
        $f2d_2 = substr($f2d, $barcode_length, length($f2d));

        $out1 = $hash_r1{$bc2};
        $out2 = $hash_r2{$bc2};

        print $out1 $f2a . $f2b_2 . $f2c . $f2d_2;
        print $out2 $f1a . $f1b . $f1c . $f1d;

    } elsif ($hash_r1{$bc1} ne "" && $hash_r1{$bc2} ne "")  { #if neither barcode from R1/R2 matches individual barcode? but then list of double barcodes in SLURM file are barcodes from real individuals that had other reads assigned to their fastqs, so I don't think I understand what double barcodes means, and this output doesn't indicate which read is the problem, so I cannot troubleshoot - try adding line here to print read info after double barcode warning?

        print "Double Barcode!\t$bc1\t$bc2\n";
        print "$f1a\n";
        print "$f1b\n";
        print "$f2a\n";
        print "$f2b\n";

    }

}
close FILE1; close FILE2;



$x=0;
while ($x <= $#commas) {
        close($hash_r1{$commas[$x]});
    close($hash_r2{$commas[$x]});
        $x++;
}
