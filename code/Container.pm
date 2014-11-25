#!/usr/bin/perl

package Container;

use assembler::Alignment;
use strict;

our @ISA = qw(Alignment);

#######################################################################################################################
sub new {
	my $class = shift;

        (@_ == 9) || die "\nError: Container::new expected 9 arguments but recieved only " . scalar(@_) . ": @_\n\n";
	my ($container, $contained, $s1, $e1, $s2, $e2, $len1, $len2, $status) = @_;

        my $self = $class->SUPER::new($container, $contained, $s1, $e1, $s2, $e2, $len1, $len2, $status);
	bless $self, $class;
	return $self;
}

#######################################################################################################################
sub type {
        my( $self ) = @_;

	return 'CONTAINER';
}

#######################################################################################################################
sub container {
        my( $self ) = @_;

	return $self->seq1;
}

#######################################################################################################################
sub contained {
        my( $self ) = @_;

	return $self->seq2;
}

1;
