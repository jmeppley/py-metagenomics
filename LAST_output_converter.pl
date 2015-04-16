#!/usr/bin/perl

use strict;
use Getopt::Long;

my $suffix;
my $e_cut;
my $bit_cut;
my $top_pct_cut;
my $nhits_cut;
my $lambda = 0.267;
my $k = 0.0410;
my $H = 0.140;
my $format = 'gene';
my $optpipe;
my $db_fasta;
my $opthelp;
my $quiet;
my $optabout;

GetOptions(
	   "s=s" => \$suffix,
	   "b=f" => \$bit_cut,
       "p=f" => \$top_pct_cut,
       "n=i" => \$nhits_cut,
	   "l=f" => \$lambda,
	   "k=f" => \$k,
	   "h=f" => \$H,
	   "f=s" => \$format,
	   "d=s" => \$db_fasta,
       "pipe" => \$optpipe,
       "help"=> \$opthelp,
       "q" => \$quiet,
       "about"=> \$optabout
);

	my $description = "This takes files in last format and coverts to m8 format.\n";
	my $usage =  "Usage $0 [list of lastex_output_files or --pipe to use STDIN/STDOUT] -b bitscore -s file suffix (default alltophits) -f (option m8 format type: blast, liz, or gene(default)) -d (optional FASTA database with descriptions)\n";


if ($opthelp) {
  print $usage;
  exit(0)
}

if ($optabout) {
  print $description;
  exit(0)
}

if (not defined $suffix) { 
    if ($format eq 'gene') {
        $suffix = "m8";
    } else {
        $suffix = "m8l"; 
    }
}

my %hit_count;
my %read_score_cut;

if (defined($optpipe)) {
    unshift @ARGV, 'STDIN';
}
if ( @ARGV < 1 ) {
	print STDERR $description;
	print STDERR $usage;
}

my %acc_info;
if (defined $db_fasta) {
	my $break = $/; 
	$/ = ">";

	my $count;

	print STDERR "Reading database file\n" unless defined $quiet;
	open (FASTA, "$db_fasta");
	while (my $entry = <FASTA>) {
		chomp $entry;
		next unless $entry;
		$count++;
		$entry =~ /^(.*?) (.*?)\n/;
		my $acc = $1;
		my $desc = $2;
		$acc_info{$acc} = $desc;
#		print STDERR "$acc\t$desc\n";
#		exit if ($count > 10);
	}
	close FASTA;
	print STDERR "Database read, $count entries found\n\n";
	$/ = "$break";
}

$/ = "\na ";

foreach my $input (@ARGV) {
    my ($LAST,$OUT);
    if ($input eq 'STDIN') {
        $LAST=\*STDIN;
        $OUT=\*STDOUT;
    } else {
        open (LAST, $input);
        my $filename = $input;
        $filename =~ s/.*\///;
        open (OUT, ">$filename\_$suffix");
        $LAST=\*LAST;
        $OUT=\*OUT;
    }

	my $m = 0;
	my $total = 0;
	my $reads = 0;
	my $last;
	my $recognized;
	my $pass;

	while (my $line = <$LAST>) {
		chomp $line;
		if ($line =~ m/^#/) {
			if ($line =~ m/Reference sequences=(\d+).*letters=(\d+)/) {
				$m = $2;
			}
			next;
		}
		$total++;
        # record should contain 3 lines
		my @splits = split(/\n/, $line);
        # First line is score
		next unless $splits[0] =~ m/score=(\d+)/;
		my $score = $1;
		my $bit = ($lambda * $score - log($k)) / log(2);		
        
        # check the bit score cutoff
		next if ($bit < $bit_cut);

        # read is 3rd line
		my @query_info = split(/\s+/, $splits[2]);		
		my $query = $query_info[1];

        # check some cutoffs
        if (defined $nhits_cut) {
            next if $hit_count{$query}>=$nhits_cut;
            $hit_count{$query}++
        }
        if (defined $top_pct_cut) {
            if (defined $read_score_cut{$query}) {
                next if $score < $read_score_cut{$query};
            } else {
                $read_score_cut{$query}=$score*(1.0-($top_pct_cut/100.0))
            }
        }
        
        # continue parsing read info
		my $q_start = $query_info[2] + 1;
		my $q_end = $query_info[2] + $query_info[3];
		if ($query_info[4] =~ m/-/) {
			$q_start = $query_info[5] - $query_info[2];
			$q_end = $q_start - $query_info[3]+1;
		}
		my $a_frac = $query_info[3] / $query_info[5];
		my $a_len = length($query_info[6]);

        # parse hit (aka subject) info
		my @subject_info = split(/\s+/, $splits[1]);
		my $subject = $subject_info[1];
#		my $a_len = $subject_info[3];
		my $s_start = $subject_info[2] + 1;
		my $s_end = $subject_info[2] + $subject_info[3];

        # compare sequences to get pct identities
		my $identities;
		my @s_align = split(//, $subject_info[6]);
		my @q_align = split(//, $query_info[6]);
        my $mismatch = 0;
        my $gaps = 0;
		for (my $i = 0; exists $s_align[$i]; $i++) {
			if ( $s_align[$i] =~ m/\Q$q_align[$i]\E/i ) {
				$identities++;
			} elsif ( $q_align[$i] =~ m/[-\/]/ ) {
                $gaps++;
            } else {
                $mismatch++;
            }
		}
		my $align_ID = 100*$identities/$a_len;
			
		my $n = int($query_info[5]/3);
		my $logKMN = int(log($k * $m * $n) / $H);
		my $mp = $m - $logKMN;
		my $np = $n - $logKMN;
		if ($np < 1/$k) {
			$np = int(1/$k);
		}
		my $eval = $mp * $np * 2**(-1*$bit);

        # count number of reads passing cutoffs
		$pass++;

		unless ($query eq $last) {
			$reads++;
			$last = $query;
		}
	
		my $desc = "NA";
		if (defined $db_fasta) {
			if (exists $acc_info{$subject}) {
				$desc = $acc_info{$subject};
				$recognized++;
			}
		}

		my $bit = sprintf("%.1f", $bit);
		my $eval = sprintf("%.1g", $eval);
        if ($format eq 'gene') {
            print $OUT "$query\tNA\t$subject\t$desc\t$align_ID\t$a_len\t$q_start\t$q_end\t$s_start\t$s_end\t$bit\t$eval\t$a_frac\n";
        } elsif ($format eq 'blast') {
            print $OUT "$query\t$subject\t$align_ID\t$a_len\t$mismatch\t$gaps\t$q_start\t$q_end\t$s_start\t$s_end\t$eval\t$bit\n";
        } elsif ($format eq 'liz') {
            print $OUT "$query\t$subject\t$desc\t$align_ID\t$a_len\t$q_start\t$q_end\t$s_start\t$s_end\t$bit\t$eval\t$a_frac\n";
        }
#		last if ($total >= 200);
	}
    if ($input ne 'STDIN') {
        close OUT;
        close LAST;
    }
    unless (defined $quiet) {
        print STDERR "$total results found, $pass hits passed cutoffs\n";
        if (defined $db_fasta) {
            print STDERR "$recognized hits recognized\n";
        }
    }
}


