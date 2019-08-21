#!/usr/bin/perl

use strict;

my $file = $ARGV[0];;

my $path = "./";

my %config = ();
my $cc = \%config;

gethash();

sub gethash {

    open(F, "< ".$path.$file);

    my $hhl = \%config;

    while (my $line = <F>) {
        chomp $line;
        my @spl = split(/\t/, $line);
        my $id = $spl[2];
        my $seek = $spl[3];
        my $seq = $spl[9];
        my $len = length($seq);
        push(@{$hhl->{$id}->{$seek}}, $len);
    }
}

my %out = ();
my $ao = \%out;

foreach my $id (sort keys %{$cc}) {
    foreach my $pos (sort keys %{$cc->{$id}}) {
        foreach my $len (sort @{$cc->{$id}->{$pos}}) {
            for (my $i = $pos; $i <= ($pos + $len); $i++) {
                if (!defined($ao->{$id}->{$i})) {
                    $ao->{$id}->{$i} = 1;
                }
                else {
                    $ao->{$id}->{$i}++;
                }
            }
        }
    }
}

foreach my $id (sort keys %{$ao}) {
    foreach my $pos (sort keys %{$ao->{$id}}) {
        print "$id\t$pos\t$ao->{$id}->{$pos}\n";
    }
}

close(F);

