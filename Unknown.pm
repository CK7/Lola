#!/usr/bin/perl

package Unknown;

use assembler::Alignment;
use strict;

our @ISA = qw(Alignment);

#######################################################################################################################
sub new {
	my $class = shift;

	# UNKNOWN mol-32-1605-031898      mol-32-1605-084167      13488   9037    Major
        (@_ == 5) || die "\nError: Container::new expected 9 arguments but recieved only " . scalar(@_) . ": @_\n\n";
	my ($seq1, $seq2, $len1, $len2, $status) = @_;

        my $self = $class->SUPER::new($seq1, $seq2, -1, -1, -1, -1, $len1, $len2, $status);

	bless $self, $class;
	return $self;
}

#######################################################################################################################
sub type {
        my( $self ) = @_;

	return 'UNKNOWN';
}

1;
