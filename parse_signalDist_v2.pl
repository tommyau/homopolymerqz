#!perl
# parse_signalDist_v2.pl
#   Parse pyrosequencing signal intensity values from Amplicon Variant Analyzer projects
#
# Chun Hang AU (chau@hksh.com)
#   Hong Kong Sanatorium and Hospital
#
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Tie::IxHash;
use IO::File;
use File::Basename;

# init and take options
my %opts = (sample => "",
	    refseq => "",
	    debug  => 0,
	    header => 0,);
my (%define, %projectdef);
tie(%define, "Tie::IxHash");

GetOptions (
	    "input=s" => \$opts{input},
	    "projectdef=s" => \$opts{projectdef},
	    "define=s" => \%define,
	    "debug" => \$opts{debug},
	    "header" => \$opts{header},
);
# check options
die "Error: please define input signalDist file name" if !defined $opts{input} || $opts{input} eq "" || ! -f $opts{input};


# defined projectdef?
{
    my $input_basename = basename($opts{input});
    my ($sample, $refseq) = ($input_basename =~ m/^([^_]+)_vs_([^_]+)\.signalDist\.txt$/);
    if (defined $opts{projectdef}) {
	die "Error: please define project definition file name" if $opts{projectdef} eq "" || !-f $opts{projectdef};
	my $projectdef_fh = IO::File->new($opts{projectdef}) || die "Error in opening file $opts{projectdef}";
	LINE: while (<$projectdef_fh>) {
	    chomp;
	    my @fields = split("\t");
	    if ($fields[0] eq "Sample") {
		$projectdef{$fields[0]}{$fields[1]} = $fields[3];
	    } elsif ($fields[0] eq "RefSeq") {
		$projectdef{$fields[0]}{$fields[1]} = $fields[3];
	    }
	$opts{sample} = $projectdef{Sample}{$sample};
	$opts{refseq} = $projectdef{RefSeq}{$refseq};
	}
    } else {
	$opts{sample} = $sample;
	$opts{refseq} = $refseq;
    }
}

my @current_position_info;
my @current_genotype_info;

# print header if instructed
if ($opts{header}) {
    print join("\t", qw(sample refseq pos gapcol direction direction_readcnt gaps ns column numgroups genotype genotype_readcnt flowvalue flowvalue_readcnt), keys(%define))."\n";
}

# extract info from input filename
{
    my $input_basename = basename($opts{input});
}

# process every line from IO::File
my $input_fh = IO::File->new($opts{input}) || die "Error in opening file $opts{input}";
LINE:while (<$input_fh>) {
    chomp;
    my $line = $_;
    if ($line =~ /^>(\d+(?:\.\d+)?) ([fr])\S+ (\d+) reads gaps=(\d+)  ns=(\d+)  column=(\d+)  numgroups=(\d+)$/) {
	# new position
	&debuginfo(join("\t", $1, $2, $3, $4, $5, $6, $7));
	my @tmp_current_position_info = ($2, $3, $4, $5, $6, $7);
	my $rawrefpos = $1;
	my ($refpos, $gapcol);
	if ($rawrefpos =~ /^(\d+)\.(\d+)$/) {
	    $refpos=$1;
	    $gapcol=$2;
	} elsif ($rawrefpos =~ /^\d+$/) {
	    $refpos=$rawrefpos;
	    $gapcol=0;
	} else {
	    die "Error: unexpected position: $rawrefpos";
	}
	@current_position_info = ($refpos, $gapcol, @tmp_current_position_info);
	#&debuginfo(Dumper(\@current_position_info));
	&debuginfo(join("\t",@current_position_info));
    } elsif ($line =~ /^\@([ACGT]) signals, (\d+) reads?$/) {
	# new genotype of current position
	&debuginfo(join("\t", $1, $2));
	@current_genotype_info = ($1, $2);
    } elsif ($line =~ /^\s+(\d+\.\d) \[ *(\d+)\]:\**$/) {
	# signal distribution of current genotype of current position
	&debuginfo(join("\t", $1, $2));
	next LINE if $2 == 0;
	print join("\t", $opts{sample}, $opts{refseq}, @current_position_info, @current_genotype_info, $1, $2, values(%define))."\n";
    } elsif ($line eq "" || $line eq "### SignalDist v1.0" || $line =~ /^\s+\.{3}$/) {
	# skip empty lines
	next LINE;
    } else {
	# unexpected lines
	die "Error: unexpected line: $line";
    }
    
}


sub debuginfo {
    if ($opts{debug}) {print "$_\n" for @_}
}