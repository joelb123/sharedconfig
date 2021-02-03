#!/usr/bin/perl

#Program takes an input file (in FASTA format, with the title field ten 
#characters or less) and converts the file into a format that rind can
#understand -- writes file "data".

$file = shift @ARGV;

$open = open (FILE, "$file");

unless ($open) {
    die "\nTarget file not found.  Specify FASTA file (label <= 10 chars) as argument.\n\n";
    }

rindFormat ($file);

sub rindFormat ($) {
    
    my $file = shift;
    my $input = "";
    my $format = "";
    open (FILE, "$file");

    while ($input = <FILE>) {
	if ($input =~ /^>/) {
	    chomp $input;
	    $input =~ s/^>//;
	    $input =~ s/ $//;
	    if (length($input) == 10) {
		$input = "$input   ";
	        }
	    if (length($input) == 9) {
		$input = "$input    ";
	        }
	    if (length($input) == 8) {
		$input = "$input     ";
	        }
	    if (length($input) == 7) {
		$input = "$input      ";
	        }
	    if (length($input) == 6) {
		$input = "$input       ";
	        }
	    if (length($input) == 5) {
		$input = "$input        ";
	        }
	    if (length($input) == 4) {
		$input = "$input         ";
	        }
	    if (length($input) == 3) {
		$input = "$input          ";
	        }
	    if (length($input) == 2) {
		$input = "$input           ";
	        }
	    if (length($input) == 1) {
		$input = "$input            ";
	        }
	    $format = "$format\n".$input;
	    }
	else {
	    chomp $input;
	    $input =~ s/^[\.\*]//;
	    $input =~ s/[\*\.]$//;
	    $input =~ s/\*/X/g;
	    $input =~ tr/./-/;
	    $format = $format.$input;
	    }
    }
    
    close (FILE);
    open (DATA, ">data");
    print DATA $format;
    close (DATA);
}
